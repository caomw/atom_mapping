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

#ifndef ATOM_MAPPING_ATOM_INDEX_H
#define ATOM_MAPPING_ATOM_INDEX_H

#include <atom_map/Atom.h>

#include <ros/ros.h>
#include <pcl/point_types.h>
#include <glog/logging.h>
#include <math.h>
#include <limits.h>
#include <functional>
#include <Eigen/Core>
#include <boost/functional/hash.hpp>

using Eigen::Vector3f;

namespace atom {
  struct AtomIndex {
    // Default constructor, destructor.
    ~AtomIndex() {}
    AtomIndex(int ii, int jj, int kk)
      : ii_(ii), jj_(jj), kk_(kk) {}

    // Construct from a 3D tuple of floats and a voxel side length.
    AtomIndex(float x, float y, float z, float voxel_size) {
      // Divide by voxel size.
      x /= voxel_size;
      y /= voxel_size;
      z /= voxel_size;

#ifdef ENABLE_DEBUG_MESSAGES
      // Check that we won't overflow long.
      float int_max = static_cast<float>(INT_MAX);
      float int_min = static_cast<float>(INT_MIN);

      if (x >= int_max || x <= int_min ||
          y >= int_max || y <= int_min ||
          z >= int_max || z <= int_min) {
        ROS_WARN("Overflowing int. Index may be inaccurate.");
      }
#endif

      // Quantize.
      const float ii = (x < 0) ? -std::ceil(std::abs(x)) : std::ceil(std::abs(x));
      const float jj = (y < 0) ? -std::ceil(std::abs(y)) : std::ceil(std::abs(y));
      const float kk = (z < 0) ? -std::ceil(std::abs(z)) : std::ceil(std::abs(z));
      ii_ = static_cast<int>(ii);
      jj_ = static_cast<int>(jj);
      kk_ = static_cast<int>(kk);
    }

    // Get the center of this bin.
    Vector3f GetBinCenter(float voxel_size) const {
      const float x = (ii_ < 0) ? static_cast<float>(ii_) + 0.5 : static_cast<float>(ii_) - 0.5;
      const float y = (jj_ < 0) ? static_cast<float>(jj_) + 0.5 : static_cast<float>(jj_) - 0.5;
      const float z = (kk_ < 0) ? static_cast<float>(kk_) + 0.5 : static_cast<float>(kk_) - 0.5;
      return voxel_size * Vector3f(x, y, z);
    }

    // Get the min corner coordinates of this bin.
    Vector3f GetMinCorner(float voxel_size) const {
      const float x = (ii_ < 0) ? static_cast<float>(ii_) : static_cast<float>(ii_) - 1.0;
      const float y = (jj_ < 0) ? static_cast<float>(jj_) : static_cast<float>(jj_) - 1.0;
      const float z = (kk_ < 0) ? static_cast<float>(kk_) : static_cast<float>(kk_) - 1.0;
      return voxel_size * Vector3f(x, y, z);
    }

    // Get the max corner coordinates of this bin.
    Vector3f GetMaxCorner(float voxel_size) const {
      const float x = (ii_ < 0) ? static_cast<float>(ii_) + 1.0 : static_cast<float>(ii_);
      const float y = (jj_ < 0) ? static_cast<float>(jj_) + 1.0 : static_cast<float>(jj_);
      const float z = (kk_ < 0) ? static_cast<float>(kk_) + 1.0 : static_cast<float>(kk_);
      return voxel_size * Vector3f(x, y, z);
    }

    // Need to overload the equality operator.
    bool operator==(const AtomIndex& other) const {
      return (ii_ == other.ii_ && jj_ == other.jj_ && kk_ == other.kk_);
    }

    // 3D coordinates.
    int ii_, jj_, kk_;
  }; // struct AtomIndex

  struct AtomIndexHasher {
    // Hash function for AtomIndex.
    size_t operator()(const AtomIndex& index) const {
      size_t combined_hash = 0;
      boost::hash_combine(combined_hash, index.ii_);
      boost::hash_combine(combined_hash, index.jj_);
      boost::hash_combine(combined_hash, index.kk_);
      return combined_hash;
    }
  }; // struct AtomIndexHasher
}

#endif
