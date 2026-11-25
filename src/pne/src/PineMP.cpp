// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "pne/PineMP.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <string>
#include <thread>
#include <unordered_map>

#include "ppl/IOPlacer.h"
#include "utl/Logger.h"
#include "sta/Network.hh"
#include "db_sta/dbSta.hh"

#include "pne/BStarTree.h"
#include "pne/CostEvaluator.h"
#include "pne/WeightScheduler.h"
#include "pne/SimulatedAnnealing.h"
#include "pne/SoftMacro.h"

namespace pne {

namespace {

bool isMovableMacro(odb::dbInst* inst)
{
  const odb::dbMasterType master_type = inst->getMaster()->getType();
  return inst->isBlock() && !inst->isFixed() && !master_type.isPad()
         && !master_type.isCover() && master_type != odb::dbMasterType::RING;
}

// Indexed by BStarTree::Anchor.
constexpr std::array<const char*, BStarTree::kNumCorners> kCornerNames
    = {"bottom-left", "bottom-right", "top-left", "top-right"};

}  // namespace

PineMP::PineMP(odb::dbDatabase* db,
         utl::Logger* logger,
         ppl::IOPlacer* io_placer,
         par::PartitionMgr* partition_mgr)
  : logger_(logger), db_(db), io_placer_(io_placer),
    partition_mgr_(partition_mgr)
{
  // Initialize components
  tree_ = std::make_unique<BStarTree>();
  cost_evaluator_ = std::make_unique<CostEvaluator>(db_, nullptr, logger_);
  weight_scheduler_ = std::make_unique<WeightScheduler>();
  sa_optimizer_ = std::make_unique<SimulatedAnnealing>(logger_);
}

PineMP::~PineMP() = default;

bool PineMP::place(int num_threads)
{
  logger_->info(utl::PNE, 1,
                "PineMP: Starting pin-aware iterative macro placement");

  num_threads_ = std::max(1, num_threads);
  
  if (!initializePlacement()) {
    logger_->error(utl::PNE, 2, "Failed to initialize placement");
    return false;
  }

  // No-op success path for designs without macros.
  if (tree_->getNumNodes() == 0) {
    logger_->info(utl::PNE, 11, "PineMP: Placement completed successfully");
    return true;
  }
  
  // Run iterative co-optimization
  runIterativeOptimization();
  
  // Apply final placement to database
  applyFinalPlacement();
  
  logger_->info(utl::PNE, 3, "PineMP: Placement completed successfully");
  
  return true;
}

bool PineMP::initializePlacement()
{
  // Get block
  odb::dbChip* chip = db_->getChip();
  if (chip == nullptr) {
    logger_->error(utl::PNE, 4, "No chip found");
    return false;
  }
  
  odb::dbBlock* block = chip->getBlock();
  if (block == nullptr) {
    logger_->error(utl::PNE, 5, "No block found");
    return false;
  }

  // Clear stale state when PineMP is invoked multiple times.
  tree_->clear();

  // Compute the placement region (core area inside IO pad ring).
  computePlacementRegion();
  
  // Collect macros
  std::vector<odb::dbInst*> macros = collectMacros();
  
  if (macros.empty() && !use_soft_macros_) {
    logger_->warn(utl::PNE, 6, "No macros found for placement");
    return true;
  }
  
  logger_->info(utl::PNE, 7, "Found {} macros for placement", macros.size());

  // Before the first pack(), which already enforces part of the rule.
  configureChannelRule();
  
  // Apply halo configuration before building the B*-Tree so tree construction
  // accounts for halo padding from the start.
  applyHalos();

  // Build B*-Tree from hard macros
  buildBStarTree(macros);

  // Attach soft macros to the tree when enabled.
  if (use_soft_macros_) {
    attachSoftMacros();
  }

  // Provide the core extent so the tree can reflect the packed layout
  // toward the right / top corners during anchor selection.
  tree_->setCoreDimensions(placement_core_.dx(), placement_core_.dy());

  // Pack the tree now that halos are configured and nodes exist.
  tree_->pack();
  logger_->info(utl::PNE, 9,
                "B*-Tree: width={}, height={}, area={}",
                tree_->getWidth(), tree_->getHeight(), tree_->getArea());
  
  // Initialize cost evaluator with network
  sta::dbNetwork* network = getNetwork();
  cost_evaluator_ = std::make_unique<CostEvaluator>(db_, network, logger_);
  cost_evaluator_->setPlacementCore(placement_core_);
  cost_evaluator_->setChannelRule(channel_rule_,
                                  fixed_channel_macros_,
                                  weight_scheduler_->getChannelWeight());
  cost_evaluator_->setObjective(objective_, spread_weight_);
  cost_evaluator_->setSoftOverlapWeight(
      use_corner_subtrees_ ? soft_overlap_weight_ : 1.0);

  // Configure weight scheduler
  weight_scheduler_->setInitialInternalWeight(initial_internal_weight_);
  weight_scheduler_->setInitialIOWeight(initial_io_weight_);
  weight_scheduler_->setFinalInternalWeight(final_internal_weight_);
  weight_scheduler_->setFinalIOWeight(final_io_weight_);
  weight_scheduler_->setNumIterations(num_iterations_);
  weight_scheduler_->reset();
  
  // Configure SA
  SAConfig sa_config;
  sa_config.initial_temperature = initial_temp_;
  sa_config.cooling_rate = cooling_rate_;
  sa_config.max_iterations = max_sa_iterations_;
  // Scale the temperature schedule to the actual cost deltas of this
  // design; a fixed temperature is either greedy descent or chaos
  // depending on the design size.  (Used only by the geometric schedule;
  // Fast-SA seeds its own stage-1 temperature.)
  sa_config.auto_calibrate_temperature = true;
  sa_config.use_fast_sa = use_fast_sa_;
  sa_config.fast_sa_accept_prob = fast_sa_accept_prob_;
  sa_config.fast_sa_c = fast_sa_c_;
  sa_config.fast_sa_k = fast_sa_k_;
  sa_config.use_slack_moves = use_slack_moves_;
  sa_config.slack_move_prob = slack_move_prob_;
  sa_optimizer_->setConfig(sa_config);

  applyPlacementWithOffset();
  if (!runPplIOPlacement("initial")) {
    logger_->info(utl::PNE,
                  12,
                  "PineMP: Continuing with internal pin assignment for initial placement");
  }

  // Classify nets and cache pin geometry.  This must run after the initial
  // ppl call so the cached IO pin positions reflect the initial assignment.
  cost_evaluator_->classifyNets(macros, tree_.get());

  return true;
}

std::vector<odb::dbInst*> PineMP::collectMacros()
{
  std::vector<odb::dbInst*> macros;
  
  odb::dbBlock* block = db_->getChip()->getBlock();
  
  for (odb::dbInst* inst : block->getInsts()) {
    if (isMovableMacro(inst)) {
      macros.push_back(inst);
    }
  }
  
  return macros;
}

void PineMP::buildBStarTree(const std::vector<odb::dbInst*>& macros)
{
  // Use placement core height for initial tree construction so macros
  // are arranged to fit within the region inside the IO pad ring.
  // Pass halo_y_ so the column-depth calculation accounts for the vertical
  // halo that will be added to every macro: without this the initial packed
  // height already overflows the core before SA starts.
  const int core_height = placement_core_.dy();

  if (use_corner_subtrees_) {
    // Each corner builds its columns in a quarter of the core.
    const auto corners = assignCorners(macros);
    tree_->buildForest(corners,
                       placement_core_.dx() / 2,
                       core_height / 2,
                       halo_x_,
                       halo_y_);
    logger_->info(utl::PNE,
                  133,
                  "Corner subtrees: {} bottom-left, {} bottom-right, {} "
                  "top-left, {} top-right macros",
                  corners[0].size(),
                  corners[1].size(),
                  corners[2].size(),
                  corners[3].size());
    return;
  }

  // Build a structurally-valid initial tree:
  //   tall macros → horizontal left-child chain (no right descendants, no overflow)
  //   short macros → vertical columns (right-child chains) to the right of talls
  tree_->buildFromMacros(macros, placement_core_.dx(), core_height,
                         halo_x_, halo_y_);
  // Don't pack yet - halos will be applied first, then pack once with halos accounted for
}

std::array<std::vector<odb::dbInst*>, BStarTree::kNumCorners>
PineMP::assignCorners(const std::vector<odb::dbInst*>& macros) const
{
  // Sum up the IO pins and pads on each macro's signal nets.
  std::unordered_map<odb::dbInst*, size_t> index;
  for (size_t i = 0; i < macros.size(); ++i) {
    index[macros[i]] = i;
  }
  std::vector<double> io_x(macros.size(), 0.0);
  std::vector<double> io_y(macros.size(), 0.0);
  std::vector<int> num_io(macros.size(), 0);
  odb::dbBlock* block = db_->getChip()->getBlock();
  for (odb::dbNet* net : block->getNets()) {
    const odb::dbSigType sig = net->getSigType();
    if (sig == odb::dbSigType::POWER || sig == odb::dbSigType::GROUND
        || sig == odb::dbSigType::CLOCK) {
      continue;
    }
    std::vector<size_t> on_net;
    std::vector<odb::Point> io;
    for (odb::dbITerm* iterm : net->getITerms()) {
      odb::dbInst* inst = iterm->getInst();
      const auto it = index.find(inst);
      int x;
      int y;
      if (it != index.end()) {
        on_net.push_back(it->second);
      } else if (inst->isPad() && (inst->isPlaced() || inst->isFixed())
                 && iterm->getAvgXY(&x, &y)) {
        io.emplace_back(x, y);
      }
    }
    if (on_net.empty()) {
      continue;
    }
    for (odb::dbBTerm* bterm : net->getBTerms()) {
      for (odb::dbBPin* bpin : bterm->getBPins()) {
        for (odb::dbBox* box : bpin->getBoxes()) {
          const odb::Rect rect = box->getBox();
          io.emplace_back((rect.xMin() + rect.xMax()) / 2,
                          (rect.yMin() + rect.yMax()) / 2);
        }
      }
    }
    std::sort(on_net.begin(), on_net.end());
    on_net.erase(std::unique(on_net.begin(), on_net.end()), on_net.end());
    for (const size_t i : on_net) {
      for (const odb::Point& point : io) {
        io_x[i] += point.x();
        io_y[i] += point.y();
        ++num_io[i];
      }
    }
  }

  auto area = [&macros](size_t i) {
    odb::dbMaster* master = macros[i]->getMaster();
    return static_cast<int64_t>(master->getWidth()) * master->getHeight();
  };
  int64_t total_area = 0;
  for (size_t i = 0; i < macros.size(); ++i) {
    total_area += area(i);
  }
  std::vector<size_t> order(macros.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&area](size_t a, size_t b) {
    return area(a) > area(b);
  });

  const odb::Rect& core = placement_core_;
  std::array<int64_t, BStarTree::kNumCorners> load{};
  std::vector<int> corner_of(macros.size(), 0);
  for (const size_t i : order) {
    // Without IO pins, the least loaded corner comes first.
    std::array<double, BStarTree::kNumCorners> key{};
    for (int c = 0; c < BStarTree::kNumCorners; ++c) {
      if (num_io[i] == 0) {
        key[c] = static_cast<double>(load[c]);
        continue;
      }
      const int x = (c & 1) != 0 ? core.xMax() : core.xMin();
      const int y = (c & 2) != 0 ? core.yMax() : core.yMin();
      key[c] = -(std::abs(x - io_x[i] / num_io[i])
                 + std::abs(y - io_y[i] / num_io[i]));
    }
    std::array<int, BStarTree::kNumCorners> corners{0, 1, 2, 3};
    std::stable_sort(corners.begin(), corners.end(), [&key](int a, int b) {
      return key[a] < key[b];
    });
    int chosen = -1;
    for (const int c : corners) {
      if (load[c] == 0
          || BStarTree::kNumCorners * (load[c] + area(i)) <= total_area) {
        chosen = c;
        break;
      }
    }
    if (chosen < 0) {
      chosen = static_cast<int>(std::min_element(load.begin(), load.end())
                                - load.begin());
    }
    load[chosen] += area(i);
    corner_of[i] = chosen;
  }

  std::array<std::vector<odb::dbInst*>, BStarTree::kNumCorners> corners;
  for (size_t i = 0; i < macros.size(); ++i) {
    corners[corner_of[i]].push_back(macros[i]);
  }
  return corners;
}

void PineMP::reportSubtrees() const
{
  if (!tree_->isForest()) {
    return;
  }
  odb::dbBlock* block = db_->getChip()->getBlock();
  std::vector<int> sizes(tree_->getNumSubtrees(), 0);
  for (const auto& node : tree_->getNodes()) {
    ++sizes[node->getSubtree()];
  }
  // In corner order, whichever subtree the anchor reflects there.
  std::vector<int> order(tree_->getNumSubtrees());
  std::iota(order.begin(), order.end(), 0);
  auto rank = [this](int s) {
    return tree_->isSubtreeCentered(s)
               ? BStarTree::kNumCorners
               : static_cast<int>(tree_->getSubtreeCorner(s));
  };
  std::stable_sort(order.begin(), order.end(), [&rank](int a, int b) {
    return rank(a) < rank(b);
  });
  std::string text;
  for (const int s : order) {
    const char* name
        = tree_->isSubtreeCentered(s)
              ? "center"
              : kCornerNames[static_cast<int>(tree_->getSubtreeCorner(s))];
    text += fmt::format("{}{} {} ({:.1f} x {:.1f} um)",
                        text.empty() ? "" : ", ",
                        name,
                        sizes[s],
                        block->dbuToMicrons(tree_->getSubtreeWidth(s)),
                        block->dbuToMicrons(tree_->getSubtreeHeight(s)));
  }
  logger_->info(utl::PNE, 134, "Subtrees: {}", text);
}

void PineMP::reportHardOverlaps() const
{
  if (!tree_->isForest()) {
    return;
  }
  // Packing keeps the macros of one subtree apart, but only the overlap
  // penalty keeps different subtrees apart.
  std::vector<odb::dbInst*> macros;
  for (const auto& node : tree_->getNodes()) {
    if (node->isHardMacro()) {
      macros.push_back(node->getInst());
    }
  }
  odb::dbBlock* block = db_->getChip()->getBlock();
  for (size_t i = 0; i < macros.size(); ++i) {
    const odb::Rect a = macros[i]->getBBox()->getBox();
    for (size_t j = i + 1; j < macros.size(); ++j) {
      const odb::Rect b = macros[j]->getBBox()->getBox();
      const int dx
          = std::min(a.xMax(), b.xMax()) - std::max(a.xMin(), b.xMin());
      const int dy
          = std::min(a.yMax(), b.yMax()) - std::max(a.yMin(), b.yMin());
      if (dx > 0 && dy > 0) {
        logger_->warn(utl::PNE,
                      135,
                      "Macros {} and {} overlap by {:.3f} x {:.3f} um",
                      macros[i]->getName(),
                      macros[j]->getName(),
                      block->dbuToMicrons(dx),
                      block->dbuToMicrons(dy));
      }
    }
  }
}

void PineMP::attachSoftMacros()
{
  // (Re-)build the soft macro manager each time placement is invoked so
  // that a fresh partitioning result is picked up automatically.
  soft_macro_mgr_ = std::make_unique<SoftMacroMgr>(db_, logger_);

  // Partition the standard cells internally by recursive bisection unless
  // the user asked to consume externally generated partition_id
  // properties (e.g. from a manual triton_part run).
  partition_tree_.clear();
  if (use_internal_partitioning_) {
    if (partition_mgr_ == nullptr) {
      logger_->warn(utl::PNE,
                    92,
                    "No partitioner available; falling back to existing "
                    "partition_id properties");
    } else {
      RecursivePartitioner partitioner(db_, partition_mgr_, logger_);

      int target = partition_target_;
      int64_t max_area = 0;
      if (partition_max_area_fraction_ > 0.0) {
        // The ceiling is on cell area: a partition with cell area A
        // inflates to a soft macro footprint of A / utilization, which
        // must stay below the requested fraction of the core area.
        const double core_area
            = static_cast<double>(placement_core_.dx()) * placement_core_.dy();
        max_area = static_cast<int64_t>(
            partition_max_area_fraction_ * core_area * soft_macro_utilization_);
      }
      if (target <= 0 && max_area <= 0) {
        // Neither knob configured: match the partition count the external
        // triton_part flow used.
        target = 10;
        logger_->info(utl::PNE,
                      93,
                      "No partitioning limits configured; defaulting to {} "
                      "partitions",
                      target);
      }

      partitioner.setTargetPartitions(target);
      partitioner.setMaxPartitionArea(max_area);
      partitioner.setMinPartitionCells(partition_min_cells_);
      partitioner.setSeed(static_cast<unsigned>(partition_seed_));
      if (partitioner.partition() > 0) {
        partition_tree_ = partitioner.getTree();
      }
    }
  }

  const int n = soft_macro_mgr_->buildFromPartitions(soft_macro_utilization_,
                                                     soft_macro_aspect_ratio_);
  if (n == 0) {
    logger_->warn(
        utl::PNE,
        87,
        "No partition_id properties found in the DB. "
        "Run partitioning (internal or triton_part) before pine_mp "
        "to generate soft macros.");
    return;
  }

  logger_->info(utl::PNE, 88, "Built {} soft macros from par partitions", n);
  soft_macro_mgr_->reportStats();

  std::vector<SoftMacro*> soft_macro_ptrs;
  soft_macro_ptrs.reserve(soft_macro_mgr_->getSoftMacros().size());
  double soft_area = 0.0;
  for (auto& sm : soft_macro_mgr_->getSoftMacros()) {
    soft_macro_ptrs.push_back(&sm);
    soft_area += static_cast<double>(sm.width) * sm.height;
  }
  int region_width = placement_core_.dx();
  int region_height = placement_core_.dy();
  if (use_corner_subtrees_ && region_width > 0 && region_height > 0) {
    // Start the centered subtree as a block of the core's aspect ratio.
    const double width
        = std::sqrt(soft_area * region_width / region_height);
    const double height = width > 0.0 ? soft_area / width : 0.0;
    region_width = std::clamp(
        static_cast<int>(std::ceil(width)), 1, region_width);
    region_height = std::clamp(
        static_cast<int>(std::ceil(height)), 1, region_height);
  }
  tree_->addSoftMacros(soft_macro_ptrs, region_width,
                         region_height, halo_x_, halo_y_);
  // Don't pack yet - halos will be applied next, then pack once with halos accounted for
}

void PineMP::reportSoftMacros() const
{
  if (soft_macro_mgr_ == nullptr || !soft_macro_mgr_->hasSoftMacros()) {
    logger_->warn(utl::PNE, 89, "No soft macros available. "
                                 "Enable with set_pine_mp_soft_macros and run pine_mp.");
    return;
  }
  soft_macro_mgr_->reportStats();

  for (const auto& sm : soft_macro_mgr_->getSoftMacros()) {
    logger_->info(utl::PNE,
                  90,
                  "  {} placed at ({}, {})",
                  sm.getName(),
                  sm.x,
                  sm.y);
  }
}

void PineMP::seedStdCellPlacement(bool report)
{
  if (!seed_std_cells_ || soft_macro_mgr_ == nullptr
      || !soft_macro_mgr_->hasSoftMacros()) {
    return;
  }

  // Collapse each partition onto its soft macro center: not a legal
  // placement, a seed for ppl and gpl, which read DB locations only for
  // PLACED instances and otherwise fall back to the core center.
  int num_seeded = 0;
  int num_skipped = 0;

  for (const auto& sm : soft_macro_mgr_->getSoftMacros()) {
    const int center_x = sm.x + sm.width / 2;
    const int center_y = sm.y + sm.height / 2;

    for (odb::dbInst* inst : sm.instances) {
      // isFixed() covers LOCKED, FIRM and COVER.
      if (inst->isFixed()) {
        ++num_skipped;
        continue;
      }

      odb::dbBox* bbox = inst->getBBox();
      inst->setLocation(center_x - static_cast<int>(bbox->getDX()) / 2,
                        center_y - static_cast<int>(bbox->getDY()) / 2);
      inst->setPlacementStatus(odb::dbPlacementStatus::PLACED);
      ++num_seeded;
    }
  }

  if (report) {
    logger_->info(utl::PNE,
                  94,
                  "Seeded {} standard cells at their soft macro centers"
                  " ({} fixed/locked cells left untouched)",
                  num_seeded,
                  num_skipped);
  }
}

void PineMP::reportPartitionTree() const
{
  pne::reportPartitionTree(partition_tree_, logger_);
}

void PineMP::applyHalos()
{
  if (halo_x_ <= 0 && halo_y_ <= 0 && macro_halo_overrides_.empty()) {
    return;
  }

  // Register per-macro overrides into the tree
  odb::dbBlock* block = db_->getChip()->getBlock();
  for (const auto& [name, halo] : macro_halo_overrides_) {
    odb::dbInst* inst = block->findInst(name.c_str());
    if (inst != nullptr) {
      tree_->setMacroHalo(inst, halo);
    } else {
      logger_->warn(utl::PNE, 80,
                    "Macro '{}' not found for halo override", name);
    }
  }

  // Configure halo settings on the tree.
  if (pin_aware_halo_) {
    tree_->computePinAwareHalos(halo_x_, halo_y_);
    logger_->info(utl::PNE, 81,
                  "Configured pin-aware halos: halo_x={}, halo_y={}",
                  halo_x_, halo_y_);
  } else {
    tree_->setUniformHalo(halo_x_, halo_y_);
    logger_->info(utl::PNE, 82,
                  "Configured uniform halos: halo_x={}, halo_y={}",
                  halo_x_, halo_y_);
  }

  if (tree_->getNumNodes() > 0) {
    // Repack with halos only once nodes are present.
    tree_->pack();
    logger_->info(utl::PNE, 83,
                  "B*-Tree after halos: width={}, height={}, area={}",
                  tree_->getWidth(), tree_->getHeight(), tree_->getArea());
    // Enforce bounds: ensure tree (with halos) fits within placement core
    enforceBoundsCompliance();
  }

}

void PineMP::configureChannelRule()
{
  odb::dbBlock* block = db_->getChip()->getBlock();
  odb::dbTech* tech = db_->getTech();
  const int grid = (tech != nullptr && tech->hasManufacturingGrid())
                       ? std::max(1, tech->getManufacturingGrid())
                       : 1;
  // Default to 2x the halo; round up to the grid so shifts stay on it.
  auto resolve = [grid](int requested, int halo) {
    const int min_channel = requested >= 0 ? requested : 2 * std::max(0, halo);
    return (min_channel + grid - 1) / grid * grid;
  };

  channel_rule_ = ChannelRule();
  channel_rule_.core = {placement_core_.xMin(),
                        placement_core_.yMin(),
                        placement_core_.xMax(),
                        placement_core_.yMax()};
  channel_rule_.min_x = resolve(min_channel_x_, halo_x_);
  channel_rule_.min_y = resolve(min_channel_y_, halo_y_);
  channel_rule_.grid = grid;

  int row_pitch = 0;
  int row_origin = 0;
  if (align_to_rows_ && findRowGrid(row_pitch, row_origin)) {
    channel_rule_.row_pitch = row_pitch;
    channel_rule_.row_origin = row_origin;
    // Whole rows, so a channel between aligned macros can be exactly the
    // minimum.
    channel_rule_.min_y
        = (channel_rule_.min_y + row_pitch - 1) / row_pitch * row_pitch;
  }

  // tapcell cuts the rows this far around every macro (its default
  // distance), and a row it touches at all is cut, so the cut reaches whole
  // rows past each edge.
  constexpr double kTapcellRowHalo = 2.0;  // um
  channel_rule_.row_cut = block->micronsToDbu(kTapcellRowHalo);
  if (row_pitch > 0) {
    channel_rule_.row_cut
        = (channel_rule_.row_cut + row_pitch - 1) / row_pitch * row_pitch;
  }

  tree_->setMinChannel(channel_rule_.min_x, channel_rule_.min_y);
  tree_->setRowGrid(channel_rule_.row_pitch,
                    channel_rule_.row_origin - placement_core_.yMin());

  fixed_channel_macros_.clear();
  fixed_channel_insts_.clear();
  if (!channel_rule_.enabled()) {
    return;
  }

  // Pre-placed macros are obstacles the rule applies to but never moves.
  for (odb::dbInst* inst : block->getInsts()) {
    if (!inst->isBlock() || !inst->isFixed()) {
      continue;
    }
    const odb::Rect bbox = inst->getBBox()->getBox();
    ChannelMacro macro;
    macro.body = {bbox.xMin(), bbox.yMin(), bbox.xMax(), bbox.yMax()};
    macro.movable = false;
    fixed_channel_macros_.push_back(macro);
    fixed_channel_insts_.push_back(inst);
  }

  if (channel_rule_.min_x > 0 || channel_rule_.min_y > 0) {
    logger_->info(utl::PNE,
                  142,
                  "Minimum macro channel: x={:.3f} um, y={:.3f} um "
                  "(macros side by side face within {:.3f} um vertically)",
                  block->dbuToMicrons(channel_rule_.min_x),
                  block->dbuToMicrons(channel_rule_.min_y),
                  block->dbuToMicrons(2 * channel_rule_.row_cut));
  }
  if (channel_rule_.row_pitch > 0) {
    logger_->info(utl::PNE,
                  145,
                  "Aligning macros to the row grid: pitch {:.3f} um, "
                  "origin y={:.3f} um",
                  block->dbuToMicrons(channel_rule_.row_pitch),
                  block->dbuToMicrons(channel_rule_.row_origin));
  }
}

bool PineMP::findRowGrid(int& pitch, int& origin) const
{
  odb::dbBlock* block = db_->getChip()->getBlock();
  auto is_core_row = [](odb::dbRow* row) {
    odb::dbSite* site = row->getSite();
    return site != nullptr && site->getClass() != odb::dbSiteClass::PAD;
  };

  // The shortest row sets the pitch; taller sites must stack on it.
  pitch = 0;
  origin = std::numeric_limits<int>::max();
  for (odb::dbRow* row : block->getRows()) {
    if (!is_core_row(row)) {
      continue;
    }
    const int height = static_cast<int>(row->getSite()->getHeight());
    pitch = pitch == 0 ? height : std::min(pitch, height);
    origin = std::min(origin, row->getBBox().yMin());
  }
  if (pitch <= 0) {
    return false;
  }
  for (odb::dbRow* row : block->getRows()) {
    if (!is_core_row(row)) {
      continue;
    }
    const int y = row->getBBox().yMin();
    const int height = static_cast<int>(row->getSite()->getHeight());
    if ((y - origin) % pitch != 0 || height % pitch != 0) {
      logger_->warn(utl::PNE,
                    148,
                    "Rows are not on a single {:.3f} um grid; macros are not "
                    "aligned to rows",
                    block->dbuToMicrons(pitch));
      return false;
    }
  }
  return true;
}

void PineMP::legalizeChannels()
{
  if (!channel_rule_.enabled()) {
    return;
  }

  std::vector<ChannelMacro> macros
      = cost_evaluator_->channelMacros(tree_.get());
  std::vector<ChannelBox> obstacles;
  for (const odb::Rect& rect : cost_evaluator_->getPlacementBlockages()) {
    obstacles.push_back({rect.xMin(), rect.yMin(), rect.xMax(), rect.yMax()});
  }

  const size_t found = findNarrowChannels(macros, channel_rule_).size();
  const size_t found_off_grid = findOffGridMacros(macros, channel_rule_).size();
  const ChannelViolations remaining
      = legalizeNarrowChannels(macros, obstacles, channel_rule_);

  // channelMacros() lists the hard nodes first, in node order.
  std::vector<std::string> names;
  size_t index = 0;
  for (const auto& node : tree_->getNodes()) {
    if (!node->isHardMacro()) {
      continue;
    }
    const ChannelBox& body = macros[index++].body;
    node->setX(body.xlo - node->getHalo().left - placement_core_.xMin());
    node->setY(body.ylo - node->getHalo().bottom - placement_core_.yMin());
    names.push_back(node->getInst()->getName());
  }
  for (odb::dbInst* inst : fixed_channel_insts_) {
    names.push_back(inst->getName());
  }

  // Fixes can trade one kind of violation for the other, so "fixed" is net.
  auto fixed = [](size_t before, size_t after) {
    return before > after ? before - after : 0;
  };
  if (channel_rule_.min_x > 0 || channel_rule_.min_y > 0) {
    logger_->info(utl::PNE,
                  140,
                  "Narrow channels: {} found, {} fixed, {} remaining",
                  found,
                  fixed(found, remaining.channels.size()),
                  remaining.channels.size());
  }
  if (channel_rule_.row_pitch > 0) {
    logger_->info(utl::PNE,
                  146,
                  "Macros off the row grid: {} found, {} fixed, {} remaining",
                  found_off_grid,
                  fixed(found_off_grid, remaining.off_grid.size()),
                  remaining.off_grid.size());
  }

  odb::dbBlock* block = db_->getChip()->getBlock();
  auto um = [block](int dbu) { return block->dbuToMicrons(dbu); };
  for (const NarrowChannel& channel : remaining.channels) {
    auto side = [&](int macro, const char* edge) {
      return macro == NarrowChannel::kCoreEdge
                 ? fmt::format("the core {} edge", edge)
                 : names[macro];
    };
    const ChannelBox& region = channel.region;
    logger_->warn(utl::PNE,
                  141,
                  "Unresolved {} channel of {:.3f} um (minimum {:.3f} um) "
                  "between {} and {} at ({:.3f}, {:.3f}) - ({:.3f}, {:.3f})",
                  channel.along_x ? "vertical" : "horizontal",
                  um(channel.gap),
                  um(channel.required),
                  side(channel.low, channel.along_x ? "left" : "bottom"),
                  side(channel.high, channel.along_x ? "right" : "top"),
                  um(region.xlo),
                  um(region.ylo),
                  um(region.xhi),
                  um(region.yhi));
  }
  for (const int macro : remaining.off_grid) {
    const ChannelBox& body = macros[macro].body;
    logger_->warn(utl::PNE,
                  147,
                  "Macro {} at y={:.3f} um is {:.3f} um above the row grid",
                  names[macro],
                  um(body.ylo),
                  um(channel_rule_.rowOffset(body)));
  }
}

void PineMP::setMacroHalo(const std::string& macro_name, const Halo& halo)
{
  macro_halo_overrides_[macro_name] = halo;
}

void PineMP::enforceBoundsCompliance()
{
  // Ensure all macros with halos stay within placement core bounds.
  // If any macro extends out of bounds, clamp its position.
  
  const int core_xmin = placement_core_.xMin();
  const int core_xmax = placement_core_.xMax();
  const int core_ymin = placement_core_.yMin();
  const int core_ymax = placement_core_.yMax();
  
  int num_clamped = 0;
  
  for (const auto& node : tree_->getNodes()) {
    int x = node->getX();
    int y = node->getY();
    const int width = node->getWidth();   // includes halos
    const int height = node->getHeight(); // includes halos
    
    bool clamped = false;
    
    // Check and clamp X bounds
    if (x < core_xmin) {
      x = core_xmin;
      clamped = true;
    }
    if (x + width > core_xmax) {
      x = core_xmax - width;
      clamped = true;
    }
    
    // Check and clamp Y bounds
    if (y < core_ymin) {
      y = core_ymin;
      clamped = true;
    }
    if (y + height > core_ymax) {
      y = core_ymax - height;
      clamped = true;
    }
    
    if (clamped) {
      node->setX(x);
      node->setY(y);
      num_clamped++;
    }
  }
  
  if (num_clamped > 0) {
    logger_->warn(utl::PNE, 84,
                  "Clamped {} macros to fit within placement core bounds",
                  num_clamped);
  }
}

void PineMP::selectBestAnchor(double internal_weight,
                             double io_weight,
                             double overlap_weight,
                             double outline_weight)
{
  if (!use_corner_anchoring_) {
    return;
  }

  static constexpr std::array<BStarTree::Anchor, 4> kAnchors = {
      BStarTree::Anchor::BOTTOM_LEFT,
      BStarTree::Anchor::BOTTOM_RIGHT,
      BStarTree::Anchor::TOP_LEFT,
      BStarTree::Anchor::TOP_RIGHT};
  static constexpr std::array<const char*, 4> kAnchorNames = {
      "bottom-left", "bottom-right", "top-left", "top-right"};

  BStarTree::Anchor best_anchor = tree_->getAnchor();
  double best_cost = std::numeric_limits<double>::max();
  int best_index = 0;

  for (size_t i = 0; i < kAnchors.size(); ++i) {
    tree_->setAnchor(kAnchors[i]);
    tree_->pack();
    const double cost = cost_evaluator_->computeCost(
        tree_.get(), internal_weight, io_weight, overlap_weight, outline_weight);
    if (cost < best_cost) {
      best_cost = cost;
      best_anchor = kAnchors[i];
      best_index = static_cast<int>(i);
    }
  }

  tree_->setAnchor(best_anchor);
  tree_->pack();

  logger_->info(utl::PNE, 130,
                "Corner anchoring: selected {} (cost {:.4f})",
                kAnchorNames[best_index], best_cost);
}

void PineMP::runMultiStartSA(double internal_weight,
                             double io_weight,
                             double overlap_weight,
                             double outline_weight,
                             int iteration)
{
  const int num_chains = num_threads_;

  if (num_chains <= 1) {
    // Legacy single-chain path: the persistent optimizer keeps its RNG
    // state across iterations, bit-compatible with previous behavior.
    auto cost_function = [this, internal_weight, io_weight, overlap_weight,
                          outline_weight](BStarTree* tree) {
      return cost_evaluator_->computeCost(
          tree, internal_weight, io_weight, overlap_weight, outline_weight);
    };
    sa_optimizer_->optimize(tree_.get(), cost_function);
    return;
  }

  logger_->info(utl::PNE, 131,
                "Multi-start SA: running {} parallel chains", num_chains);

  // Clone the starting tree per chain up front, on this thread (node
  // construction reads the DB for instance orientations).
  struct Chain
  {
    std::unique_ptr<BStarTree> tree;
    double cost = 0.0;
  };
  std::vector<Chain> chains(num_chains);
  for (int i = 0; i < num_chains; ++i) {
    chains[i].tree = tree_->clone();
  }

  const SAConfig base_config = sa_optimizer_->getConfig();

  // Worker threads share the DB strictly read-only (master dimensions, pin
  // geometry): SA never writes the DB — orientations are node-local and
  // placements are committed only after the winner is adopted.  Each chain
  // gets its own tree clone, cost evaluator copy (evaluation mutates cached
  // net bounding boxes), and independently seeded annealer.  Only chain 0
  // logs, so the logger has a single writer.
  std::vector<std::thread> workers;
  workers.reserve(num_chains);
  for (int i = 0; i < num_chains; ++i) {
    workers.emplace_back([this, i, iteration, &chains, &base_config,
                          internal_weight, io_weight, overlap_weight,
                          outline_weight]() {
      CostEvaluator evaluator = *cost_evaluator_;

      SAConfig config = base_config;
      config.quiet = (i != 0);

      // Distinct, reproducible seed per (iteration, chain).
      const unsigned seed = 42u + static_cast<unsigned>(iteration) * 7919u
                            + static_cast<unsigned>(i) * 104729u;
      SimulatedAnnealing sa(logger_, seed);
      sa.setConfig(config);

      BStarTree* tree = chains[i].tree.get();
      sa.optimize(tree,
                  [&evaluator, internal_weight, io_weight, overlap_weight,
                   outline_weight](BStarTree* t) {
                    return evaluator.computeCost(t,
                                                 internal_weight,
                                                 io_weight,
                                                 overlap_weight,
                                                 outline_weight);
                  });

      tree->pack();
      chains[i].cost = evaluator.computeCost(
          tree, internal_weight, io_weight, overlap_weight, outline_weight);
    });
  }
  for (std::thread& worker : workers) {
    worker.join();
  }

  int best = 0;
  double worst_cost = chains[0].cost;
  for (int i = 1; i < num_chains; ++i) {
    if (chains[i].cost < chains[best].cost) {
      best = i;
    }
    worst_cost = std::max(worst_cost, chains[i].cost);
  }

  tree_->copyStateFrom(*chains[best].tree);
  tree_->pack();

  logger_->info(utl::PNE, 132,
                "Multi-start SA: chain {} wins with cost {:.4f} "
                "(worst chain: {:.4f})",
                best, chains[best].cost, worst_cost);
}

void PineMP::runIterativeOptimization()
{
  logger_->info(utl::PNE, 40,
                "Starting iterative co-optimization with {} iterations",
                num_iterations_);
  
  std::vector<odb::dbInst*> macros = collectMacros();

  // Keep a fixed validation objective so iterations remain comparable.
  const double validation_internal_weight = final_internal_weight_;
  const double validation_io_weight = final_io_weight_;
  const double validation_overlap_weight = weight_scheduler_->getOverlapWeight();
  const double validation_outline_weight = weight_scheduler_->getOutlineWeight();

  tree_->pack();

  // Fix the wirelength normalization baselines to the initial placement so
  // every iteration (and the validation objective) is scored on the same
  // scale.  SPREAD minimizes only the soft parts of the nets.
  const bool spread = objective_ == Objective::SPREAD;
  cost_evaluator_->computeWeightedWirelength(tree_.get(), 1.0, 1.0);
  if (spread) {
    cost_evaluator_->setWirelengthBaselines(
        cost_evaluator_->getInternalSoftWirelength(),
        cost_evaluator_->getIOSoftWirelength());
    logger_->info(utl::PNE,
                  48,
                  "Spread objective: hard macro HPWL weight {:.2f} over {} "
                  "nets",
                  spread_weight_,
                  cost_evaluator_->getNumHardNets());
  } else {
    cost_evaluator_->setWirelengthBaselines(
        cost_evaluator_->getInternalWirelength(),
        cost_evaluator_->getIOWirelength());
  }

  // Anchor the initial cluster before recording the first global best so
  // the baseline placement is already corner-aligned.
  selectBestAnchor(validation_internal_weight,
                   validation_io_weight,
                   validation_overlap_weight,
                   validation_outline_weight);

  double global_best_cost
      = cost_evaluator_->computeCost(tree_.get(),
                                     validation_internal_weight,
                                     validation_io_weight,
                                     validation_overlap_weight,
                                     validation_outline_weight);
  tree_->saveSnapshot(BStarTree::SnapshotSlot::GLOBAL);

  for (int iter = 0; iter < num_iterations_; ++iter) {
    double internal_weight = weight_scheduler_->getInternalWeight();
    double io_weight = weight_scheduler_->getIOWeight();
    double overlap_weight = weight_scheduler_->getOverlapWeight();
    double outline_weight = weight_scheduler_->getOutlineWeight();
    
    logger_->info(utl::PNE, 41,
                  "--- Iteration {} ---", iter + 1);
    logger_->info(utl::PNE, 42,
                  "Weights - Internal: {:.2f}, IO: {:.2f}, "
                  "Overlap: {:.0e}, Outline: {:.0e}",
                  internal_weight, io_weight, overlap_weight, outline_weight);
    
    // Run SA optimization: one chain, or num_threads parallel seeded
    // chains with the best result adopted.
    runMultiStartSA(internal_weight, io_weight, overlap_weight,
                    outline_weight, iter);

    // Choose the best corner anchoring for the arrangement SA converged to,
    // before committing the placement and reassigning pins.  Evaluated with
    // the fixed validation objective so the choice is consistent with the
    // global-best tracking below.
    selectBestAnchor(validation_internal_weight,
                     validation_io_weight,
                     validation_overlap_weight,
                     validation_outline_weight);

    // Apply current placement to database (for pin assignment)
    applyPlacementWithOffset();

    // Refresh wirelength statistics for the state we just committed.
    cost_evaluator_->computeWeightedWirelength(tree_.get(), 1.0, 1.0);

    // Track best placement over all outer iterations with a fixed objective.
    const double validation_cost
        = cost_evaluator_->computeCost(tree_.get(),
                                       validation_internal_weight,
                                       validation_io_weight,
                                       validation_overlap_weight,
                                       validation_outline_weight);
    if (validation_cost < global_best_cost) {
      global_best_cost = validation_cost;
      tree_->saveSnapshot(BStarTree::SnapshotSlot::GLOBAL);
    }
    
    // Report statistics
    double total_wl = cost_evaluator_->getTotalWirelength();
    double internal_wl = cost_evaluator_->getInternalWirelength();
    double io_wl = cost_evaluator_->getIOWirelength();
    
    logger_->info(utl::PNE, 43,
                  "Wirelength - Total: {:.0f}, Internal: {:.0f}, IO: {:.0f}",
                  total_wl, internal_wl, io_wl);
    if (spread) {
      logger_->info(utl::PNE,
                    49,
                    "Spread - Soft internal: {:.0f}, Soft IO: {:.0f}, Hard "
                    "extension: {:.0f}",
                    cost_evaluator_->getInternalSoftWirelength(),
                    cost_evaluator_->getIOSoftWirelength(),
                    cost_evaluator_->getHardExtension());
    }
    logger_->info(utl::PNE, 44,
                  "Floorplan - Width: {}, Height: {}, Area: {}",
                  tree_->getWidth(), tree_->getHeight(), tree_->getArea());
    reportSubtrees();
    
    // Update weights for next iteration
    if (use_adaptive_io_weighting_) {
      // Adaptively adjust weights based on IO wirelength proportion, of the
      // part of the nets the weights apply to.
      if (spread) {
        weight_scheduler_->updateWeightsBasedOnIOProportion(
            cost_evaluator_->getIOSoftWirelength(),
            cost_evaluator_->getInternalSoftWirelength());
      } else {
        weight_scheduler_->updateWeightsBasedOnIOProportion(io_wl,
                                                            internal_wl);
      }
      logger_->info(utl::PNE, 47,
                    "Adaptive weights - Internal: {:.2f}, IO: {:.2f}",
                    weight_scheduler_->getInternalWeight(),
                    weight_scheduler_->getIOWeight());
    } else {
      weight_scheduler_->nextIteration();
    }
    
    // Refresh pin assignment based on current placement.
    if (iter < num_iterations_ - 1) {
      logger_->info(utl::PNE, 45, "Updating pin assignment for next iteration");
      const double soft_io_wl = cost_evaluator_->getIOSoftWirelength();
      const double extension = cost_evaluator_->getHardExtension();
      const bool pins_moved = runPplIOPlacement("iteration");

      // Re-cache net geometry so the next SA run sees the new IO pin
      // positions.  This is the feedback path of the co-optimization.
      cost_evaluator_->classifyNets(macros, tree_.get());

      // How far ppl moved the pins toward (or away from) the macros they
      // connect to, as the cost model sees them.
      if (pins_moved) {
        cost_evaluator_->computeWeightedWirelength(tree_.get(), 1.0, 1.0);
        logger_->info(utl::PNE,
                      66,
                      "Pin placement changed the IO wirelength from {:.0f} "
                      "to {:.0f}",
                      io_wl,
                      cost_evaluator_->getIOWirelength());
        if (spread) {
          logger_->info(utl::PNE,
                        67,
                        "Pin placement changed the soft IO wirelength from "
                        "{:.0f} to {:.0f} and the hard extension from {:.0f} "
                        "to {:.0f}",
                        soft_io_wl,
                        cost_evaluator_->getIOSoftWirelength(),
                        extension,
                        cost_evaluator_->getHardExtension());
        }
      }
    }
  }

  // Apply the best state seen across all outer iterations.
  tree_->restoreSnapshot(BStarTree::SnapshotSlot::GLOBAL);
  tree_->pack();
  
  logger_->info(utl::PNE, 46, "Iterative co-optimization completed");
}

void PineMP::applyFinalPlacement()
{
  logger_->info(utl::PNE, 50, "Applying final placement to database");

  // Collect before locking below: LOCKED instances count as fixed and
  // would no longer be reported by collectMacros().
  std::vector<odb::dbInst*> macros = collectMacros();

  tree_->pack();
  legalizeChannels();
  applyPlacementWithOffset();
  reportHardOverlaps();

  // Write halos to the database
  for (const auto& node : tree_->getNodes()) {
    if (node->isSoftMacro()) {
      // Soft macros are virtual; no DB blockage is created for them here.
      // Future work: create a fence/region for the soft macro area.
      continue;
    }
    const Halo& halo = node->getHalo();
    if (halo.hasNonZero()) {
      int x, y;
      node->getInst()->getLocation(x, y);
      // Use the placed instance location (x,y) and the node-specific
      // macro dimensions + halo extents to build the blockage rectangle.
      const int left = x - halo.left;
      const int bottom = y - halo.bottom;
      const int right = x + node->getMacroWidth() + halo.right;
      const int top = y + node->getMacroHeight() + halo.top;

      odb::dbBlockage* blockage = odb::dbBlockage::create(
          node->getInst()->getBlock(), left, bottom, right, top,
          node->getInst());
      blockage->setSoft();
    }
  }

  // Mark hard macros as LOCKED so downstream tools treat them as final.
  // Soft macro positions are recorded in their SoftMacro structs.
  for (const auto& node : tree_->getNodes()) {
    if (node->isHardMacro()) {
      node->getInst()->setPlacementStatus(odb::dbPlacementStatus::LOCKED);
    }
  }

  if (!runPplIOPlacement("final")) {
    logger_->info(utl::PNE,
                  58,
                  "PineMP: Keeping internal pin assignment results for final placement");
  }

  // Ensure final reported wirelength matches the final applied state,
  // including the IO pin positions from the final ppl run.
  cost_evaluator_->classifyNets(macros, tree_.get());
  cost_evaluator_->computeWeightedWirelength(tree_.get(), 1.0, 1.0);
  
  // Report final statistics
  double total_wl = cost_evaluator_->getTotalWirelength();
  double internal_wl = cost_evaluator_->getInternalWirelength();
  double io_wl = cost_evaluator_->getIOWirelength();
  
  logger_->info(utl::PNE, 51,
                "Final Results:");
  logger_->info(utl::PNE, 52,
                "  Total Wirelength: {:.0f}", total_wl);
  logger_->info(utl::PNE, 53,
                "  Internal Wirelength: {:.0f}", internal_wl);
  logger_->info(utl::PNE, 54,
                "  IO Wirelength: {:.0f}", io_wl);
  if (objective_ == Objective::SPREAD) {
    logger_->info(utl::PNE,
                  59,
                  "  Soft Wirelength: {:.0f}, Hard Extension: {:.0f}",
                  cost_evaluator_->getInternalSoftWirelength()
                      + cost_evaluator_->getIOSoftWirelength(),
                  cost_evaluator_->getHardExtension());
  }
  logger_->info(utl::PNE, 55,
                "  Floorplan Width: {}", tree_->getWidth());
  logger_->info(utl::PNE, 56,
                "  Floorplan Height: {}", tree_->getHeight());
  logger_->info(utl::PNE, 57,
                "  Floorplan Area: {}", tree_->getArea());
  reportSubtrees();

  // The final ppl call already seeded the std cells, at these same points;
  // this seeds them when ppl did not run, and reports the counts.
  seedStdCellPlacement(true);
}

sta::dbNetwork* PineMP::getNetwork()
{
  if (sta_ != nullptr) {
    return sta_->getDbNetwork();
  }
  return nullptr;
}

bool PineMP::runPplIOPlacement(const char* stage_label)
{
  if (io_placer_ == nullptr) {
    return false;
  }

  odb::dbChip* chip = db_->getChip();
  if (chip == nullptr || chip->getBlock() == nullptr) {
    return false;
  }

  odb::dbBlock* block = chip->getBlock();
  odb::dbTech* tech = db_->getTech();
  if (tech == nullptr) {
    return false;
  }

  std::vector<odb::dbTechLayer*> horizontal_layers;
  std::vector<odb::dbTechLayer*> vertical_layers;
  for (odb::dbTechLayer* layer : tech->getLayers()) {
    if (layer->getType() != odb::dbTechLayerType::ROUTING) {
      continue;
    }

    odb::dbTrackGrid* track_grid = block->findTrackGrid(layer);
    if (track_grid == nullptr) {
      continue;
    }

    if (layer->getDirection() == odb::dbTechLayerDir::HORIZONTAL
        && track_grid->getNumGridPatternsY() > 0) {
      horizontal_layers.push_back(layer);
    } else if (layer->getDirection() == odb::dbTechLayerDir::VERTICAL
               && track_grid->getNumGridPatternsX() > 0) {
      vertical_layers.push_back(layer);
    }
  }

  auto by_level = [](odb::dbTechLayer* lhs, odb::dbTechLayer* rhs) {
    return lhs->getRoutingLevel() < rhs->getRoutingLevel();
  };
  std::sort(horizontal_layers.begin(), horizontal_layers.end(), by_level);
  std::sort(vertical_layers.begin(), vertical_layers.end(), by_level);

  if (horizontal_layers.empty() || vertical_layers.empty()) {
    logger_->warn(utl::PNE,
                  61,
                  "PineMP: Skipping PPL IO placement at {} stage (missing routing layers/tracks)",
                  stage_label);
    return false;
  }

  // ppl reads unplaced instances at the die center, so without the seed it
  // would see every std cell there.
  seedStdCellPlacement(false);
  const std::vector<std::pair<odb::dbInst*, odb::dbPlacementStatus>> hidden
      = hideHardMacros(stage_label);

  bool placed = false;
  try {
    for (odb::dbTechLayer* layer : horizontal_layers) {
      io_placer_->addHorLayer(layer);
    }
    for (odb::dbTechLayer* layer : vertical_layers) {
      io_placer_->addVerLayer(layer);
    }
    io_placer_->runHungarianMatching(false);
    placed = true;
  } catch (const std::exception& e) {
    logger_->warn(utl::PNE,
                  63,
                  "PineMP: PPL IO placement failed at {} stage: {}",
                  stage_label,
                  e.what());
  } catch (...) {
    logger_->warn(utl::PNE,
                  64,
                  "PineMP: PPL IO placement failed at {} stage",
                  stage_label);
  }

  for (const auto& [inst, status] : hidden) {
    inst->setPlacementStatus(status);
  }
  if (placed) {
    logger_->info(utl::PNE,
                  62,
                  "PineMP: Updated IO placement with PPL at {} stage",
                  stage_label);
  }
  return placed;
}

std::vector<std::pair<odb::dbInst*, odb::dbPlacementStatus>>
PineMP::hideHardMacros(const char* stage_label)
{
  std::vector<std::pair<odb::dbInst*, odb::dbPlacementStatus>> hidden;
  if (!hide_macros_from_ppl_) {
    return hidden;
  }

  const odb::Rect die = db_->getChip()->getBlock()->getDieArea();
  int num_hard = 0;
  for (const auto& node : tree_->getNodes()) {
    if (!node->isHardMacro()) {
      continue;
    }
    ++num_hard;
    // ppl keeps pins off the die edge a placed macro touches, so those stay.
    odb::dbInst* inst = node->getInst();
    if (!die.inside(inst->getBBox()->getBox())) {
      continue;
    }
    hidden.emplace_back(inst, inst->getPlacementStatus());
    inst->setPlacementStatus(odb::dbPlacementStatus::UNPLACED);
  }

  logger_->info(utl::PNE,
                65,
                "PineMP: Hiding {} of {} hard macros from PPL at {} stage",
                hidden.size(),
                num_hard,
                stage_label);
  return hidden;
}

void PineMP::computePlacementRegion()
{
  odb::dbBlock* block = db_->getChip()->getBlock();
  const odb::Rect die_area = block->getDieArea();
  const odb::Rect core_area = block->getCoreArea();

  // If rows define a core area strictly inside the die, use it directly.
  if (core_area.dx() < die_area.dx() || core_area.dy() < die_area.dy()) {
    placement_core_ = core_area;
    logger_->info(utl::PNE,
                  70,
                  "Using core area for macro placement: ({}, {}) - ({}, {})",
                  placement_core_.xMin(),
                  placement_core_.yMin(),
                  placement_core_.xMax(),
                  placement_core_.yMax());
    return;
  }

  // Fallback: compute border margins from IO pad / cover / ring geometry.
  int left_margin = 0;
  int right_margin = 0;
  int bottom_margin = 0;
  int top_margin = 0;

  for (odb::dbInst* inst : block->getInsts()) {
    const odb::dbMasterType master_type = inst->getMaster()->getType();
    if (!master_type.isPad() && !master_type.isCover()
        && master_type != odb::dbMasterType::RING) {
      continue;
    }
    if (!inst->isPlaced() && !inst->isFixed()) {
      continue;
    }

    const odb::Rect bbox = inst->getBBox()->getBox();

    // Determine the closest die edge and record the intrusion depth.
    const int dist_left = std::abs(bbox.xMin() - die_area.xMin());
    const int dist_right = std::abs(die_area.xMax() - bbox.xMax());
    const int dist_bottom = std::abs(bbox.yMin() - die_area.yMin());
    const int dist_top = std::abs(die_area.yMax() - bbox.yMax());

    const int min_dist
        = std::min({dist_left, dist_right, dist_bottom, dist_top});

    if (min_dist == dist_left) {
      left_margin = std::max(left_margin, bbox.xMax() - die_area.xMin());
    } else if (min_dist == dist_right) {
      right_margin = std::max(right_margin, die_area.xMax() - bbox.xMin());
    } else if (min_dist == dist_bottom) {
      bottom_margin = std::max(bottom_margin, bbox.yMax() - die_area.yMin());
    } else {
      top_margin = std::max(top_margin, die_area.yMax() - bbox.yMin());
    }
  }

  placement_core_.set_xlo(die_area.xMin() + left_margin);
  placement_core_.set_ylo(die_area.yMin() + bottom_margin);
  placement_core_.set_xhi(die_area.xMax() - right_margin);
  placement_core_.set_yhi(die_area.yMax() - top_margin);

  // Sanity check: core must have positive area.
  if (placement_core_.dx() <= 0 || placement_core_.dy() <= 0) {
    logger_->warn(utl::PNE,
                  71,
                  "Computed placement core has non-positive dimensions, "
                  "falling back to die area");
    placement_core_ = die_area;
  }

  if (placement_core_ != die_area) {
    logger_->info(
        utl::PNE,
        72,
        "Computed IO pad border margins: left={}, right={}, bottom={}, top={}",
        left_margin,
        right_margin,
        bottom_margin,
        top_margin);
    logger_->info(utl::PNE,
                  73,
                  "Macro placement region: ({}, {}) - ({}, {})",
                  placement_core_.xMin(),
                  placement_core_.yMin(),
                  placement_core_.xMax(),
                  placement_core_.yMax());
  }
}

void PineMP::applyPlacementWithOffset()
{
  tree_->applyPlacement();

  // Offset macros from B*-tree origin (0,0) into the placement core region.
  const int offset_x = placement_core_.xMin();
  const int offset_y = placement_core_.yMin();

  if (offset_x == 0 && offset_y == 0) {
    return;
  }

  for (const auto& node : tree_->getNodes()) {
    if (node->isSoftMacro()) {
      // No dbInst: applyPlacement() left the origin in the SoftMacro struct
      // in tree-local coordinates, so it needs the same shift.
      SoftMacro* sm = node->getSoftMacro();
      sm->x += offset_x;
      sm->y += offset_y;
      continue;
    }
    odb::dbInst* inst = node->getInst();
    if (!inst) {
      continue;
    }
    int x, y;
    inst->getLocation(x, y);
    inst->setLocation(x + offset_x, y + offset_y);
  }
}

}  // namespace pne
