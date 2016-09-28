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

#ifndef ATOM_MAPPING_ATOM_HASH_GRID_INDEX_H
#define ATOM_MAPPING_ATOM_HASH_GRID_INDEX_H

#include <atom_map/Atom.h>

#include <pcl/point_types.h>
#include <glog/logging.h>
#include <math.h>
#include <limits.h>
#include <functional>

namespace atom {
  struct AtomIndex {
    // Construct from a 3D tuple of floats and a voxel side length.
    ~AtomIndex() {}
    AtomIndex(float x, float y, float z, float voxel_size) {
      // Divide by voxel size.
      x /= voxel_size;
      y /= voxel_size;
      z /= voxel_size;

#ifdef ENABLE_DEBUG_MESSAGES
      // Check that we won't overflow long.
      float long_max = static_cast<float>(LONG_MAX);
      float long_min = static_cast<float>(LONG_MIN);

      if (x >= long_max || x <= long_min ||
          y >= long_max || y <= long_min ||
          z >= long_max || z <= long_min) {
        ROS_WARN("Overflowing long. Index may be inaccurate.")
      }
#endif

      // Quantize.
      ii_ = (x < 0) ? -ceil(abs(x)) : ceil(abs(x));
      jj_ = (y < 0) ? -ceil(abs(y)) : ceil(abs(y));
      kk_ = (z < 0) ? -ceil(abs(z)) : ceil(abs(z));
    }

    // Need to overload the equality operator.
    bool operator==(const AtomIndex& other) const {
      return (ii_ == other.ii_ && jj_ == other.jj_ && kk_ == other.kk_);
    }

    // A 3D tuple of longs.
    long ii_, jj_, kk_;
  }; // struct AtomIndex

  struct AtomIndexHasher {
    // Hash function for AtomIndex.
    size_t operator()(const AtomIndex& index) const {
      return (std::hash<long>(index.ii_) ^
              std::hash<long>(index.jj_) ^
              std::hash<long>(index.kk_));
    }
  }; // struct AtomIndexHasher
}

#endif
