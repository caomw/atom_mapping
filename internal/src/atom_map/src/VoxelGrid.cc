/*
 * Copyright (c) 2016, The Regents of the University of California (Regents).
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *    1. Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *    3. Neither the name of the copyright holder nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS AS IS
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * Please contact the author(s) of this library if you have any questions.
 * Authors: David Fridovich-Keil   ( dfk@eecs.berkeley.edu )
 *          Erik Nelson            ( eanelson@eecs.berkeley.edu )
 */

#include <atom_map/VoxelGrid.h>

namespace atom {
  VoxelGrid::VoxelGrid(float leaf_size, float limit_min, float limit_max)
    : limit_min_(limit_min), limit_max_(limit_max),
      inverse_leaf_size_(1.0 / leaf_size) {}
  VoxelGrid::~VoxelGrid() {}

  // Filter a point cloud.
  bool VoxelGrid::Filter(const PointCloud::ConstPtr& raw,
                         const PointCloud::Ptr filtered) const {
    CHECK_NOTNULL(raw.get());
    CHECK_NOTNULL(filtered.get());
    filtered->points.clear();
    std::unordered_map<unsigned long, const pcl::PointXYZ> registry;

    // Set up bounding box.
    const float min_bb = floor(limit_min_ * inverse_leaf_size_);
    const float max_bb = floor(limit_max_ * inverse_leaf_size_);
    const unsigned long side_bb = static_cast<unsigned long>(max_bb - min_bb);

    // Check that we're not going to overflow unsigned long.
    if (max_bb - min_bb >= static_cast<float>(ULONG_MAX) ||
        pow(max_bb - min_bb, 3) >= static_cast<float>(ULONG_MAX)) {
      ROS_WARN("Could overflow unsigned long. Did not filter.");
      pcl::copyPointCloud(*raw, *filtered);
      return false;
    }

    // Map every point in the raw cloud to a voxel, assign unique indices
    // for every voxel, and insert into the registry.
    for (unsigned long ii = 0; ii < raw->points.size(); ii++) {
      const pcl::PointXYZ p = raw->points[ii];

      // Check not out of bounds.
      if (isnan(p.x) || isnan(p.y) || isnan(p.z) ||
          isinf(p.x) || isinf(p.y) || isinf(p.z)) {
        ROS_WARN("Point was NAN or INF.");
        continue;
      }

      if (p.x >= limit_max_ || p.y >= limit_max_ || p.z >= limit_max_ ||
          p.x <= limit_min_ || p.y <= limit_min_ || p.z <= limit_min_) {
        ROS_WARN("Point out of bounds.");
        continue;
      }

      // Get 3D voxel coordinates by dividing by leaf size and flooring.
      const unsigned long voxel_x =
        static_cast<unsigned long>(floor(p.x * inverse_leaf_size_) - min_bb);
      const unsigned long voxel_y =
        static_cast<unsigned long>(floor(p.y * inverse_leaf_size_) - min_bb);
      const unsigned long voxel_z =
        static_cast<unsigned long>(floor(p.z * inverse_leaf_size_) - min_bb);

      // Get a single 3D index for these coordinates.
      const unsigned long voxel_idx =
        voxel_x + voxel_y * side_bb + voxel_z * side_bb * side_bb;

      // Add to registry.
      registry.insert({voxel_idx, p});
    }

    // Loop over every key in the registry and populate the filtered cloud.
    for (const auto& pair : registry) {
      filtered->points.push_back(pair.second);
    }

    return true;
  }

} // namespace atom
