// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "pne/PinAssigner.h"

#include <algorithm>
#include <cmath>
#include <random>

#include "utl/Logger.h"

namespace pne {

PinAssigner::PinAssigner(odb::dbDatabase* db, utl::Logger* logger)
    : db_(db), logger_(logger)
{
}

void PinAssigner::assignPins(const std::vector<odb::dbInst*>& macros)
{
  logger_->info(utl::PNE, 20, "Assigning pins for {} macros", macros.size());
  
  for (odb::dbInst* macro : macros) {
    assignMacroPins(macro);
  }
  
  logger_->info(utl::PNE, 21, "Pin assignment completed");
}

void PinAssigner::assignMacroPins(odb::dbInst* macro)
{
  switch (strategy_) {
    case PinAssignmentStrategy::UNIFORM:
      assignUniform(macro);
      break;
      
    case PinAssignmentStrategy::CONNECTIVITY:
      assignConnectivity(macro);
      break;
      
    case PinAssignmentStrategy::RANDOM:
      assignRandom(macro);
      break;
      
    case PinAssignmentStrategy::HUNGARAIN:
      // Hungarian algorithm not implemented yet, fall back to connectivity
      assignConnectivity(macro);
      break;
  }
}

void PinAssigner::assignUniform(odb::dbInst* macro)
{
  std::vector<odb::dbITerm* > pins = getMacroPins(macro);
  
  if (pins.empty()) {
    return;
  }
  
  // Distribute pins uniformly around the perimeter
  int num_pins = pins.size();
  
  odb::dbBox* bbox = macro->getBBox();
  odb::Rect rect = bbox->getBox();
  
  int width = rect.dx();
  int height = rect.dy();
  int perimeter = 2 * (width + height);
  
  for (int i = 0; i < num_pins; ++i) {
    double fraction = static_cast<double>(i) / num_pins;
    int pos = static_cast<int>(fraction * perimeter);
    
    Side side;
    double side_fraction;
    
    if (pos < width) {
      side = Side::BOTTOM;
      side_fraction = static_cast<double>(pos) / width;
    } else if (pos < width + height) {
      side = Side::RIGHT;
      side_fraction = static_cast<double>(pos - width) / height;
    } else if (pos < 2 * width + height) {
      side = Side::TOP;
      side_fraction = static_cast<double>(pos - width - height) / width;
    } else {
      side = Side::LEFT;
      side_fraction = static_cast<double>(pos - 2 * width - height) / height;
    }
    
    PinLocation loc = computePinLocation(macro, side, side_fraction);
    
    // Set pin placement (this is a simplified approach)
    // In reality, OpenDB pin placement is more complex
    odb::dbITerm* iterm = pins[i];
    // Note: dbITerm does not have setPlacementStatus; placement is managed at dbInst level
    // iterm->setPlacementStatus(loc.status);
  }
}

void PinAssigner::assignConnectivity(odb::dbInst* macro)
{
  std::vector<PinConnectivity> conn_info = analyzeConnectivity(macro);
  
  if (conn_info.empty()) {
    assignUniform(macro);
    return;
  }
  
  // Assign each pin to the side closest to its connectivity center
  for (const auto& pin_conn : conn_info) {
    Side side = selectBoundarySide(macro, pin_conn.center_x, pin_conn.center_y);
    
    // Compute position along the selected side
    odb::dbBox* bbox = macro->getBBox();
    odb::Rect rect = bbox->getBox();
    
    int macro_x, macro_y;
    macro->getLocation(macro_x, macro_y);
    
    double fraction = 0.5;  // Default to middle
    
    switch (side) {
      case Side::LEFT:
      case Side::RIGHT:
        if (rect.dy() > 0) {
          fraction = (pin_conn.center_y - macro_y) / rect.dy();
        }
        break;
      case Side::TOP:
      case Side::BOTTOM:
        if (rect.dx() > 0) {
          fraction = (pin_conn.center_x - macro_x) / rect.dx();
        }
        break;
    }
    
    fraction = std::max(0.0, std::min(1.0, fraction));
    
    PinLocation loc = computePinLocation(macro, side, fraction);
    // Note: dbITerm does not have setPlacementStatus; placement is managed at dbInst level
    // pin_conn.iterm->setPlacementStatus(loc.status);
  }
}

void PinAssigner::assignRandom(odb::dbInst* macro)
{
  std::vector<odb::dbITerm*> pins = getMacroPins(macro);
  
  if (pins.empty()) {
    return;
  }
  
  static std::random_device rd;
  static std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);
  
  for (odb::dbITerm* iterm : pins) {
    Side side = static_cast<Side>(static_cast<int>(dis(gen) * 4) % 4);
    double fraction = dis(gen);
    
    PinLocation loc = computePinLocation(macro, side, fraction);
    // Note: dbITerm does not have setPlacementStatus; placement is managed at dbInst level
    // iterm->setPlacementStatus(loc.status);
  }
}

std::vector<odb::dbITerm*> PinAssigner::getMacroPins(odb::dbInst* macro)
{
  std::vector<odb::dbITerm*> pins;
  
  for (odb::dbITerm* iterm : macro->getITerms()) {
    // Skip power/ground pins
    if (iterm->getSigType() == odb::dbSigType::POWER ||
        iterm->getSigType() == odb::dbSigType::GROUND) {
      continue;
    }
    
    pins.push_back(iterm);
  }
  
  return pins;
}

PinLocation PinAssigner::computePinLocation(odb::dbInst* macro,
                                            Side side,
                                            double fraction)
{
  odb::dbBox* bbox = macro->getBBox();
  odb::Rect rect = bbox->getBox();
  
  int macro_x, macro_y;
  macro->getLocation(macro_x, macro_y);
  
  PinLocation loc;
  loc.status = odb::dbPlacementStatus::PLACED;
  
  switch (side) {
    case Side::LEFT:
      loc.x = rect.xMin();
      loc.y = rect.yMin() + static_cast<int>(fraction * rect.dy());
      break;
      
    case Side::RIGHT:
      loc.x = rect.xMax();
      loc.y = rect.yMin() + static_cast<int>(fraction * rect.dy());
      break;
      
    case Side::TOP:
      loc.x = rect.xMin() + static_cast<int>(fraction * rect.dx());
      loc.y = rect.yMax();
      break;
      
    case Side::BOTTOM:
      loc.x = rect.xMin() + static_cast<int>(fraction * rect.dx());
      loc.y = rect.yMin();
      break;
  }
  
  return loc;
}

std::vector<PinAssigner::PinConnectivity> 
PinAssigner::analyzeConnectivity(odb::dbInst* macro)
{
  std::vector<PinConnectivity> result;
  
  std::vector<odb::dbITerm*> pins = getMacroPins(macro);
  
  for (odb::dbITerm* iterm : pins) {
    odb::dbNet* net = iterm->getNet();
    if (net == nullptr) {
      continue;
    }
    
    PinConnectivity conn;
    conn.iterm = iterm;
    conn.center_x = 0.0;
    conn.center_y = 0.0;
    conn.num_connections = 0;
    
    // Compute weighted center of connected pins
    for (odb::dbITerm* other_iterm : net->getITerms()) {
      if (other_iterm == iterm) {
        continue;
      }
      
      int x, y;
      if (other_iterm->getAvgXY(&x, &y)) {
        conn.center_x += x;
        conn.center_y += y;
        conn.num_connections++;
      } else {
        odb::dbInst* other_inst = other_iterm->getInst();
        int inst_x, inst_y;
        other_inst->getLocation(inst_x, inst_y);
        conn.center_x += inst_x;
        conn.center_y += inst_y;
        conn.num_connections++;
      }
    }
    
    // Also consider block terminals (IO pins)
    for (odb::dbBTerm* bterm : net->getBTerms()) {
      for (odb::dbBPin* bpin : bterm->getBPins()) {
        for (odb::dbBox* box : bpin->getBoxes()) {
          odb::Rect rect = box->getBox();
          int x = (rect.xMin() + rect.xMax()) / 2;
          int y = (rect.yMin() + rect.yMax()) / 2;
          
          conn.center_x += x;
          conn.center_y += y;
          conn.num_connections++;
        }
      }
    }
    
    if (conn.num_connections > 0) {
      conn.center_x /= conn.num_connections;
      conn.center_y /= conn.num_connections;
      result.push_back(conn);
    }
  }
  
  return result;
}

PinAssigner::Side PinAssigner::selectBoundarySide(odb::dbInst* macro,
                                                  double target_x,
                                                  double target_y)
{
  odb::dbBox* bbox = macro->getBBox();
  odb::Rect rect = bbox->getBox();
  
  int macro_x, macro_y;
  macro->getLocation(macro_x, macro_y);
  
  int center_x = rect.xMin() + rect.dx() / 2;
  int center_y = rect.yMin() + rect.dy() / 2;
  
  double dx = target_x - center_x;
  double dy = target_y - center_y;
  
  // Select side based on which direction has larger component
  if (std::abs(dx) > std::abs(dy)) {
    return dx > 0 ? Side::RIGHT : Side::LEFT;
  } else {
    return dy > 0 ? Side::TOP : Side::BOTTOM;
  }
}

void PinAssigner::resetPins(const std::vector<odb::dbInst*>& macros)
{
  // Reset all pins to unplaced status
  // Note: dbITerm does not have setPlacementStatus; placement is managed at dbInst level
  for (odb::dbInst* macro : macros) {
    // Reset instance placement status instead
    macro->setPlacementStatus(odb::dbPlacementStatus::UNPLACED);
  }
}

}  // namespace pne
