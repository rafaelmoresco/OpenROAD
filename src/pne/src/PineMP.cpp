// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "pne/PineMP.h"

#include <algorithm>

#include "ppl/IOPlacer.h"
#include "utl/Logger.h"
#include "sta/Network.hh"
#include "db_sta/dbSta.hh"

#include "pne/BStarTree.h"
#include "pne/CostEvaluator.h"
#include "pne/WeightScheduler.h"
#include "pne/PinAssigner.h"
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

}  // namespace

PineMP::PineMP(odb::dbDatabase* db,
         utl::Logger* logger,
         ppl::IOPlacer* io_placer)
  : logger_(logger), db_(db), io_placer_(io_placer)
{
  // Initialize components
  tree_ = std::make_unique<BStarTree>();
  cost_evaluator_ = std::make_unique<CostEvaluator>(db_, nullptr, logger_);
  weight_scheduler_ = std::make_unique<WeightScheduler>();
  pin_assigner_ = std::make_unique<PinAssigner>(db_, logger_);
  sa_optimizer_ = std::make_unique<SimulatedAnnealing>(logger_);
}

PineMP::~PineMP() = default;

bool PineMP::place(int num_threads)
{
  logger_->info(utl::PNE, 1, 
                "PineMP: Starting pin-aware iterative macro placement");
  
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
  
  // Build B*-Tree from hard macros
  buildBStarTree(macros);

  // Attach soft macros to the tree when enabled.
  if (use_soft_macros_) {
    attachSoftMacros();
  }

  // Apply halo configuration to the tree
  applyHalos();
  
  // Initialize cost evaluator with network
  sta::dbNetwork* network = getNetwork();
  cost_evaluator_ = std::make_unique<CostEvaluator>(db_, network, logger_);
  cost_evaluator_->setPlacementCore(placement_core_);
  cost_evaluator_->classifyNets(macros);
  
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
  sa_optimizer_->setConfig(sa_config);
  
  // Initial pin assignment (uniform distribution)
  logger_->info(utl::PNE, 8, "Performing initial pin assignment");
  pin_assigner_->assignPins(macros);

  tree_->pack();
  applyPlacementWithOffset();
  if (!runPplIOPlacement("initial")) {
    logger_->info(utl::PNE,
                  12,
                  "PineMP: Continuing with internal pin assignment for initial placement");
  }
  
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

  // Build a structurally-valid initial tree:
  //   tall macros → horizontal left-child chain (no right descendants, no overflow)
  //   short macros → vertical columns (right-child chains) to the right of talls
  tree_->buildFromMacros(macros, core_height, halo_y_);
  // Initial packing
  tree_->pack();
  
  logger_->info(utl::PNE, 9, 
                "Initial B*-Tree: width={}, height={}, area={}",
                tree_->getWidth(), tree_->getHeight(), tree_->getArea());
}

void PineMP::attachSoftMacros()
{
  // (Re-)build the soft macro manager each time placement is invoked so
  // that a fresh triton_part result is picked up automatically.
  soft_macro_mgr_ = std::make_unique<SoftMacroMgr>(db_, logger_);

  const int n = soft_macro_mgr_->buildFromPartitions(soft_macro_utilization_,
                                                     soft_macro_aspect_ratio_);
  if (n == 0) {
    logger_->warn(
        utl::PNE,
        87,
        "No partition_id properties found in the DB. "
        "Run triton_part before pine_mp to generate soft macros.");
    return;
  }

  logger_->info(utl::PNE, 88, "Built {} soft macros from par partitions", n);
  soft_macro_mgr_->reportStats();

  // Append soft macro nodes to the existing B*-Tree (after hard macros).
  for (auto& sm : soft_macro_mgr_->getSoftMacros()) {
    tree_->addSoftMacro(&sm);
  }

  tree_->pack();
  logger_->info(utl::PNE,
                9,
                "B*-Tree with soft macros: width={}, height={}, area={}",
                tree_->getWidth(),
                tree_->getHeight(),
                tree_->getArea());
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
                  89,
                  "  {} placed at ({}, {})",
                  sm.getName(),
                  sm.x,
                  sm.y);
  }
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

  // Apply halos
  if (pin_aware_halo_) {
    tree_->computePinAwareHalos(halo_x_, halo_y_);
    logger_->info(utl::PNE, 81,
                  "Applied pin-aware halos: halo_x={}, halo_y={}",
                  halo_x_, halo_y_);
  } else {
    tree_->setUniformHalo(halo_x_, halo_y_);
    logger_->info(utl::PNE, 82,
                  "Applied uniform halos: halo_x={}, halo_y={}",
                  halo_x_, halo_y_);
  }

  // Repack with halos
  tree_->pack();
  logger_->info(utl::PNE, 83,
                "B*-Tree after halos: width={}, height={}, area={}",
                tree_->getWidth(), tree_->getHeight(), tree_->getArea());
}

void PineMP::setMacroHalo(const std::string& macro_name, const Halo& halo)
{
  macro_halo_overrides_[macro_name] = halo;
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
    
    // Define cost function with current weights
    auto cost_function = [&](BStarTree* tree) {
      return cost_evaluator_->computeCost(tree,
                                          internal_weight,
                                          io_weight,
                                          overlap_weight,
                                          outline_weight);
    };
    
    // Run SA optimization
    sa_optimizer_->optimize(tree_.get(), cost_function);
    
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
    logger_->info(utl::PNE, 44,
                  "Floorplan - Width: {}, Height: {}, Area: {}",
                  tree_->getWidth(), tree_->getHeight(), tree_->getArea());
    
    // Update weights for next iteration
    weight_scheduler_->nextIteration();
    
    // Refresh pin assignment based on current placement.
    if (iter < num_iterations_ - 1) {
      logger_->info(utl::PNE, 45, "Updating pin assignment for next iteration");
      if (!runPplIOPlacement("iteration")) {
        pin_assigner_->assignPins(macros);
      }
      
      // Reclassify nets after pin reassignment
      cost_evaluator_->classifyNets(macros);
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
  
  tree_->pack();
  applyPlacementWithOffset();

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
      odb::dbBlockage* blockage = odb::dbBlockage::create(node->getInst()->getBlock(),
                                        x - halo.left,
                                        y - halo.bottom,
                                        x + node->getWidth(),
                                        node->getY() + node->getHeight() + halo_y_,
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

  // Ensure final reported wirelength matches the final applied state.
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
  logger_->info(utl::PNE, 55,
                "  Floorplan Width: {}", tree_->getWidth());
  logger_->info(utl::PNE, 56,
                "  Floorplan Height: {}", tree_->getHeight());
  logger_->info(utl::PNE, 57,
                "  Floorplan Area: {}", tree_->getArea());
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

  try {
    for (odb::dbTechLayer* layer : horizontal_layers) {
      io_placer_->addHorLayer(layer);
    }
    for (odb::dbTechLayer* layer : vertical_layers) {
      io_placer_->addVerLayer(layer);
    }
    io_placer_->runHungarianMatching(false);
    logger_->info(utl::PNE,
                  62,
                  "PineMP: Updated IO placement with PPL at {} stage",
                  stage_label);
    return true;
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

  return false;
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
    odb::dbInst* inst = node->getInst();
    // Check if its Soft or Hard macro
    if (!inst) {
      continue;
    }
    int x, y;
    inst->getLocation(x, y);
    inst->setLocation(x + offset_x, y + offset_y);
  }
}

void PineMP::setPinAssignmentStrategy(const std::string& strategy)
{
  if (strategy == "uniform") {
    pin_assigner_->setStrategy(PinAssignmentStrategy::UNIFORM);
  } else if (strategy == "connectivity") {
    pin_assigner_->setStrategy(PinAssignmentStrategy::CONNECTIVITY);
  } else if (strategy == "random") {
    pin_assigner_->setStrategy(PinAssignmentStrategy::RANDOM);
  } else if (strategy == "hungarian") {
    pin_assigner_->setStrategy(PinAssignmentStrategy::HUNGARAIN);
  } else {
    logger_->warn(utl::PNE, 60,
                  "Unknown pin assignment strategy '{}', using uniform",
                  strategy);
    pin_assigner_->setStrategy(PinAssignmentStrategy::UNIFORM);
  }
}

}  // namespace pne
