// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "pne/PineMP.h"

#include "utl/Logger.h"
#include "sta/Network.hh"
#include "dbSta/dbSta.hh"

#include "pne/BStarTree.h"
#include "pne/CostEvaluator.h"
#include "pne/WeightScheduler.h"
#include "pne/PinAssigner.h"
#include "pne/SimulatedAnnealing.h"

namespace pne {

PineMP::PineMP(odb::dbDatabase* db, utl::Logger* logger)
    : logger_(logger), db_(db)
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
  
  // Collect macros
  std::vector<odb::dbInst*> macros = collectMacros();
  
  if (macros.empty()) {
    logger_->warn(utl::PNE, 6, "No macros found for placement");
    return false;
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
  
  return true;
}

std::vector<odb::dbInst*> PineMP::collectMacros()
{
  std::vector<odb::dbInst*> macros;
  
  odb::dbBlock* block = db_->getChip()->getBlock();
  
  for (odb::dbInst* inst : block->getInsts()) {
    // Check if instance is a macro (block instance)
    if (inst->isBlock()) {
      macros.push_back(inst);
    }
  }
  
  return macros;
}

void PineMP::buildBStarTree(const std::vector<odb::dbInst*>& macros)
{
  tree_->clear();
  
  for (odb::dbInst* macro : macros) {
    tree_->addMacro(macro);
  }
  
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
    
    // Reassign pins based on current placement
    if (iter < num_iterations_ - 1) {
      logger_->info(utl::PNE, 45, "Reassigning pins for next iteration");
      pin_assigner_->assignPins(macros);
      
      // Reclassify nets after pin reassignment
      cost_evaluator_->classifyNets(macros);
    }
  }
  
  logger_->info(utl::PNE, 46, "Iterative co-optimization completed");
}

void PineMP::applyFinalPlacement()
{
  logger_->info(utl::PNE, 50, "Applying final placement to database");
  
  tree_->applyPlacement();
  
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
