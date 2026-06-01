// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "pne/SoftMacro.h"

#include <algorithm>
#include <cmath>
#include <map>

#include "utl/Logger.h"

namespace pne {

int64_t SoftMacro::getCellArea() const
{
  int64_t total = 0;
  for (odb::dbInst* inst : instances) {
    odb::dbMaster* master = inst->getMaster();
    total += static_cast<int64_t>(master->getWidth())
             * static_cast<int64_t>(master->getHeight());
  }
  return total;
}

SoftMacroMgr::SoftMacroMgr(odb::dbDatabase* db, utl::Logger* logger)
    : db_(db), logger_(logger)
{
}

bool SoftMacroMgr::isStdCell(odb::dbInst* inst)
{
  const odb::dbMasterType master_type = inst->getMaster()->getType();
  // Exclude hard blocks, pads, covers, and ring cells.
  return !inst->isBlock() && !master_type.isPad() && !master_type.isCover()
         && master_type != odb::dbMasterType::RING;
}

std::pair<int, int> SoftMacroMgr::computeDimensions(int64_t cell_area,
                                                     double utilization,
                                                     double aspect_ratio,
                                                     int manufacturing_grid)
{
  if (cell_area <= 0 || utilization <= 0.0 || aspect_ratio <= 0.0) {
    return {manufacturing_grid, manufacturing_grid};
  }

  // soft macro area = cell_area / utilization
  // width  = sqrt(area_sm * aspect_ratio)
  // height = sqrt(area_sm / aspect_ratio)
  const double area_sm = static_cast<double>(cell_area) / utilization;
  int width
      = static_cast<int>(std::ceil(std::sqrt(area_sm * aspect_ratio)));
  int height
      = static_cast<int>(std::ceil(std::sqrt(area_sm / aspect_ratio)));

  // Snap dimensions to manufacturing grid
  auto snapToGrid = [manufacturing_grid](int value) {
    return ((value + manufacturing_grid - 1) / manufacturing_grid) * manufacturing_grid;
  };

  width = snapToGrid(std::max(1, width));
  height = snapToGrid(std::max(1, height));

  return {width, height};
}

int SoftMacroMgr::buildFromPartitions(double target_utilization,
                                      double aspect_ratio)
{
  soft_macros_.clear();

  odb::dbChip* chip = db_->getChip();
  if (chip == nullptr) {
    return 0;
  }

  odb::dbBlock* block = chip->getBlock();
  if (block == nullptr) {
    return 0;
  }

  // Group instances by partition_id.
  std::map<int, std::vector<odb::dbInst*>> groups;
  for (odb::dbInst* inst : block->getInsts()) {
    if (!isStdCell(inst)) {
      continue;
    }
    odb::dbIntProperty* prop = odb::dbIntProperty::find(inst, "partition_id");
    if (prop == nullptr) {
      continue;
    }
    groups[prop->getValue()].push_back(inst);
  }

  if (groups.empty()) {
    return 0;
  }

  soft_macros_.reserve(groups.size());
  for (auto& [part_id, insts] : groups) {
    SoftMacro sm;
    sm.partition_id = part_id;
    sm.instances = std::move(insts);

    const int64_t cell_area = sm.getCellArea();
    odb::dbTech* tech = db_->getTech();
    int manufacturing_grid = tech->getManufacturingGrid();
    const auto [w, h]
        = computeDimensions(cell_area, target_utilization, aspect_ratio, manufacturing_grid);
    sm.width = w;
    sm.height = h;

    soft_macros_.push_back(std::move(sm));
  }

  return static_cast<int>(soft_macros_.size());
}

void SoftMacroMgr::reportStats() const
{
  logger_->info(
      utl::PNE, 91, "Soft macros from partitions: {}", soft_macros_.size());
  for (const auto& sm : soft_macros_) {
    logger_->info(utl::PNE,
                  86,
                  "  Partition {}: {} cells, {}x{} DBU (area={})",
                  sm.partition_id,
                  sm.instances.size(),
                  sm.width,
                  sm.height,
                  sm.getArea());
  }
}

}  // namespace pne
