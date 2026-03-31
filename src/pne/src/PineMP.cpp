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
  
  // Collect macros
  std::vector<odb::dbInst*> macros = collectMacros();
  
  if (macros.empty()) {
    logger_->warn(utl::PNE, 6, "No macros found for placement");
    return true;
  }
  
  logger_->info(utl::PNE, 7, "Found {} macros for placement", macros.size());
  
  // Build B*-Tree
  buildBStarTree(macros);
  
  // Initialize cost evaluator with network
  sta::dbNetwork* network = getNetwork();
  cost_evaluator_ = std::make_unique<CostEvaluator>(db_, network, logger_);
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
  tree_->applyPlacement();
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
  odb::dbBlock* block = db_->getChip()->getBlock();
  odb::Rect die_area = block->getDieArea();
  const int die_height = die_area.dy();

  // Build a structurally-valid initial tree:
  //   tall macros → horizontal left-child chain (no right descendants, no overflow)
  //   short macros → vertical columns (right-child chains) to the right of talls
  tree_->buildFromMacros(macros, die_height);
  // Initial packing
  tree_->pack();
  
  logger_->info(utl::PNE, 9, 
                "Initial B*-Tree: width={}, height={}, area={}",
                tree_->getWidth(), tree_->getHeight(), tree_->getArea());
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
    tree_->applyPlacement();

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
  tree_->applyPlacement();

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
