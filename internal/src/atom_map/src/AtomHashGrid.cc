
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
#include <unordered_set>

using Eigen::Vector3f;

namespace atom {
  AtomHashGrid::~AtomHashGrid() {}
  AtomHashGrid::AtomHashGrid()
    : initialized_(false),
      voxel_size_(1.0),
      max_distance_(-std::numeric_limits<float>::infinity()),
      min_distance_(std::numeric_limits<float>::infinity()) {}

  // Set voxel side length to be larger than an atomic radius. This may be
  // tuned later, but preliminary timing suggests that it should be on the
  // order of 1-2 atomic radii.
  void AtomHashGrid::SetAtomicRadius(float atomic_radius) {
    voxel_size_ = 3.0 * atomic_radius; //1.9 * atomic_radius/sqrt(3.0);
    initialized_ = true;
  }

  // Populate a set of bins overlapping a particular ball.
  // Currently, this only works for r < voxel_size_.
  bool AtomHashGrid::GetOverlappingBins(float x, float y, float z, float r,
                                        std::vector<AtomIndex>* bins) const {
    if (!initialized_) {
#ifdef ENABLE_DEBUG_MESSAGES
      ROS_WARN("AtomHashGrid was not initialized. Overlapping bin computation will fail.");
#endif
      return false;
    }

    CHECK_NOTNULL(bins);
    bins->clear();

    // For a bin with center c and side length l, and ball of center a and
    // radius r, we know that the ball will extend beyond the boundary of the
    // bin along dimension i in the negative direction if c_i - a_i + r > l/2,
    // and similarly in the positive direction if a_i - c_i + r > l/2.
    const AtomIndex quant(x, y, z, voxel_size_);
    const Vector3f center = quant.GetBinCenter(voxel_size_);
    bins->push_back(quant);

    // Check each dimension, forward and backward.
    if (center(0) - x + r > 0.5 * voxel_size_)
      bins->push_back(AtomIndex(quant.ii_ - 1, quant.jj_, quant.kk_));
    if (x - center(0) + r > 0.5 * voxel_size_)
      bins->push_back(AtomIndex(quant.ii_ + 1, quant.jj_, quant.kk_));
    if (center(1) - y + r > 0.5 * voxel_size_)
      bins->push_back(AtomIndex(quant.ii_, quant.jj_ - 1, quant.kk_));
    if (y - center(1) + r > 0.5 * voxel_size_)
      bins->push_back(AtomIndex(quant.ii_, quant.jj_ + 1, quant.kk_));
    if (center(2) - z + r > 0.5 * voxel_size_)
      bins->push_back(AtomIndex(quant.ii_, quant.jj_, quant.kk_ - 1));
    if (z - center(2) + r > 0.5 * voxel_size_)
      bins->push_back(AtomIndex(quant.ii_, quant.jj_, quant.kk_ + 1));

    return true;
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

    // Find the set of all bins overlapping this search domain.
    std::vector<AtomIndex> search_bins;
    if (!GetOverlappingBins(x, y, z, r, &search_bins)) {
#ifdef ENABLE_DEBUG_MESSAGES
      ROS_WARN("Unable to get overlapping bins. Radius search will fail.");
#endif
      return false;
    }

    // Iterate over all search bins, and check all unique Atoms.
    std::unordered_set<unsigned int> unique_atom_ids;
    for (const auto& bin : search_bins) {
      if (map_.count(bin) == 0)
        continue;

      for (const auto& atom_id : map_.at(bin)) {
        // Check if we have seen this Atom before.
        if (unique_atom_ids.count(atom_id) > 0)
          continue;

        unique_atom_ids.insert(atom_id);

        // Check proximity.
        const Atom::Ptr atom = registry_[atom_id];
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

    // Get the bin IDs which overlap this Atom.
    const Vector3f p = atom->GetPosition();
    std::vector<AtomIndex> overlapping_bins;
    if (!GetOverlappingBins(p(0), p(1), p(2), atom->GetRadius(),
                            &overlapping_bins)) {
#ifdef ENABLE_DEBUG_MESSAGES
      ROS_WARN("Unable to get overlapping bins. Atom insertion will fail.");
#endif
      return false;
    }

    // Insert this Atom into all overlapping bins.
    for (const auto& bin : overlapping_bins) {
      if (map_.count(bin) == 0) {
        std::list<unsigned int> atom_index =
          {static_cast<unsigned int>(registry_.size())};
        map_.insert({bin, atom_index});
      } else {
        map_.at(bin).push_back(static_cast<unsigned int>(registry_.size()));
      }
    }

    // Add to registry.
    registry_.push_back(atom);

    // Update max and min distances.
    const float sdf = atom->GetSignedDistance();
    if (sdf > max_distance_) {
      max_distance_ = sdf;
    } else if (sdf < min_distance_) {
      min_distance_ = sdf;
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
