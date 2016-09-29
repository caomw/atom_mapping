
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

#include <atom_map/AtomHashGrid.h>

#include <Eigen/Core>
#include <limits.h>
#include <math.h>
#include <glog/logging.h>

using Eigen::Vector3f;

namespace atom {
  AtomHashGrid::~AtomHashGrid() {}
  AtomHashGrid::AtomHashGrid()
    : initialized_(false),
      voxel_size_(1.0),
      max_distance_(-std::numeric_limits<float>::infinity()),
      min_distance_(std::numeric_limits<float>::infinity()) {}

  // Set atomic radius slightly smaller than tightest packing.
  // This will throw an "initialized" flag, without which
  // nothing else will work.
  void AtomHashGrid::SetAtomicRadius(float atomic_radius) {
    voxel_size_ = 1.9 * atomic_radius/sqrt(3.0);
    initialized_ = true;
  }

  // Nearest neighbor queries.
  bool AtomHashGrid::RadiusSearch(float x, float y, float z, float r,
                                  std::vector<Atom::Ptr>* neighbors) {
    if (!initialized_) {
#ifdef ENABLE_DEBUG_MESSAGES
      ROS_WARN("AtomHashGrid was not initialized. Radius search will fail.");
#endif
      return false;
    }

    CHECK_NOTNULL(neighbors);
    neighbors->clear();

    // Get indices for corners of bounding box.
    const AtomIndex min_corner(x - r, y - r, z - r, voxel_size_);
    const AtomIndex max_corner(x + r, y + r, z + r, voxel_size_);

    // Iterate over all voxel indices in the box, and collect all
    // atoms within range of the specified coordinates.
    for (long ii = min_corner.ii_; ii <= max_corner.ii_; ii++) {
      for (long jj = min_corner.jj_; jj <= max_corner.jj_; jj++) {
        for (long kk = min_corner.kk_; kk <= max_corner.kk_; kk++) {
          const AtomIndex query(ii, jj, kk);

          if (map_.count(query) > 0) {
            // If there is an Atom here, check proximity.
            const Atom::Ptr atom = registry_[map_.at(query)];
            const Vector3f p = atom->GetPosition();

            const float dx = p(0) - x;
            const float dy = p(1) - y;
            const float dz = p(2) - z;

            if (dx*dx + dy*dy + dz*dz <= r*r) {
              // Atom is in range, so add to neighbors.
              neighbors->push_back(atom);
            }
          }
        }
      }
    }

    return true;
  }

  bool AtomHashGrid::RadiusSearch(const pcl::PointXYZ& p, float r,
                                  std::vector<Atom::Ptr>* neighbors) {
    return RadiusSearch(p.x, p.y, p.z, r, neighbors);
  }

  // Insert a new Atom.
  bool AtomHashGrid::Insert(const Atom::Ptr& atom) {
    if (!initialized_) {
#ifdef ENABLE_DEBUG_MESSAGES
      ROS_WARN("AtomHashGrid was not initialized. Atom insertion will fail.");
#endif
        return false;
    }

    const Vector3f p = atom->GetPosition();
    const AtomIndex index(p(0), p(1), p(2), voxel_size_);

    // Check that this Atom does not land in the same bin as an existing Atom.
    if (map_.count(index) > 0) {
#ifdef ENABLE_DEBUG_MESSAGES
      ROS_WARN("AtomMap already contains an atom in this bin. Did not insert.");
#endif
      return false;
    } else {
      map_.insert({index, registry_.size()});
      registry_.push_back(atom);

      // Update max and min distances.
      const float sdf = atom->GetSignedDistance();
      if (sdf > max_distance_) {
        max_distance_ = sdf;
      } else if (sdf < min_distance_) {
        min_distance_ = sdf;
      }
    }

    return true;
  }

  // Return a list of all Atoms in the map.
  const std::vector<Atom::Ptr>& AtomHashGrid::GetAtoms() const {
    return registry_;
  }

  // Return the number of Atoms in the map.
  size_t AtomHashGrid::Size() const {
    return registry_.size();
  }

  // Return the maximum and minimum distances of any Atom to the surface.
  float AtomHashGrid::GetMaxDistance() const {
    return min_distance_;
  }

  float AtomHashGrid::GetMinDistance() const {
    return max_distance_;
  }
} // namespace atom
