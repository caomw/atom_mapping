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

#include <atom_map/AtomKdtree.h>

#include <iostream>
#include <Eigen/Core>

namespace atom {
  AtomKdtree::AtomKdtree()
    : max_distance_(-std::numeric_limits<float>::infinity()),
      min_distance_(std::numeric_limits<float>::infinity()),
      occupancy_mode_(true) {}
  AtomKdtree::~AtomKdtree() {}

  // Set occupancy mode.
  void AtomKdtree::SetOccupancyMode(bool mode) { occupancy_mode_ = mode; }

  // Nearest neighbor queries.
  bool AtomKdtree::GetKNearestNeighbors(float x, float y, float z, size_t k,
                                        std::vector<Atom::Ptr>* neighbors) {
    CHECK_NOTNULL(neighbors);
    neighbors->clear();

    if (index_ == nullptr) {
      VLOG(1) << "Index has not been built. Points must be added before "
              <<  "querying the kd tree";
      return false;
    }

    // Convert the input point to the FLANN format.
    const int kNumColumns = 3;
    flann::Matrix<float> flann_query(new float[kNumColumns], 1, kNumColumns);
    flann_query[0][0] = x;
    flann_query[0][1] = y;
    flann_query[0][2] = z;

    // Search the kd tree for the nearest neighbor to the query.
    std::vector< std::vector<int> > query_match_indices;
    std::vector< std::vector<float> > query_distances;

    const int num_neighbors_found =
      index_->knnSearch(flann_query, query_match_indices,
                        query_distances, static_cast<int>(k),
                        flann::SearchParams(-1, 0.0, false));

    // Assign output.
    for (size_t ii = 0; ii < num_neighbors_found; ii++)
      neighbors->push_back(registry_[ query_match_indices[0][ii] ]);

    // Remember to delete query point.
    delete[] flann_query.ptr();

    return true;
  }

  // Nearest neighbor queries.
  bool AtomKdtree::GetKNearestNeighbors(const pcl::PointXYZ& p, size_t k,
                                        std::vector<Atom::Ptr>* neighbors) {
    return GetKNearestNeighbors(p.x, p.y, p.z, k, neighbors);
  }

  // Radius searching.
  bool AtomKdtree::RadiusSearch(float x, float y, float z, float r,
                                std::vector<Atom::Ptr>* neighbors) {
    CHECK_NOTNULL(neighbors);
    neighbors->clear();

    if (index_ == nullptr) {
      VLOG(1) << "Index has not been built. Points must be added before "
              <<  "querying the kd tree";
      return false;
    }

    // Convert the input point to the FLANN format.
    const int kNumColumns = 3;
    flann::Matrix<float> flann_query(new float[kNumColumns], 1, kNumColumns);
    flann_query[0][0] = x;
    flann_query[0][1] = y;
    flann_query[0][2] = z;

    // Search the kd tree for the nearest neighbor to the query.
    std::vector< std::vector<int> > query_match_indices;
    std::vector< std::vector<float> > query_distances;

    // FLANN checks Euclidean distance squared, so we pass in r * r.
    int num_neighbors_found =
      index_->radiusSearch(flann_query, query_match_indices,
                           query_distances, r * r,
                           flann::SearchParams(-1, 0.0, false));
    // Assign output.
    for (size_t ii = 0; ii < num_neighbors_found; ii++)
      neighbors->push_back(registry_[ query_match_indices[0][ii] ]);

    // Remember to delete query point.
    delete[] flann_query.ptr();

    return true;
  }

  // Radius searching.
  bool AtomKdtree::RadiusSearch(const pcl::PointXYZ& p, float r,
                                std::vector<Atom::Ptr>* neighbors) {
    return RadiusSearch(p.x, p.y, p.z, r, neighbors);
  }

  // Insert a new Atom.
  bool AtomKdtree::Insert(const Atom::Ptr& atom) {
    CHECK_NOTNULL(atom.get());

    // Copy the input point into FLANN's Matrix type.
    const int kNumColumns = 3;
    flann::Matrix<float> flann_point(atom->GetPosition().data(), 1, kNumColumns);

    // If this is the first point in the index, create the index and exit.
    if (index_ == nullptr) {
      // Single kd-tree.
      const int kMaxPointsPerLeaf = 20;
      index_.reset(new flann::KDTreeSingleIndex<flann::L2<float>>(
          flann_point, flann::KDTreeSingleIndexParams(kMaxPointsPerLeaf)));

      index_->buildIndex();
    } else {
      // If the index is already created, add the data point to the index.
      // Rebuild every time the index floats in size to occasionally rebalance
      // the kdtree.
      const int kRebuildThreshold = 2;
      index_->addPoints(flann_point, kRebuildThreshold);
    }

    // Add point to registry.
    registry_.push_back(atom);

    // Update max and min distances.
    if (!occupancy_mode_) {
      const float sdf =
        std::static_pointer_cast<SdfAtom>(atom)->GetSignedDistance();
      if (sdf > max_distance_)
        max_distance_ = sdf;
      else if (sdf < min_distance_)
        min_distance_ = sdf;
    }

    return true;
  }

  // Return a list of all Atoms in the map.
  const std::vector<Atom::Ptr>& AtomKdtree::GetAtoms() const {
    return registry_;
  }

  // Return the size of this tree.
  size_t AtomKdtree::Size() const {
    return registry_.size();
  }

  // Return the maximum and minimum distances of any Atom to the surface.
  float AtomKdtree::GetMaxDistance() const {
    return max_distance_;
  }

  float AtomKdtree::GetMinDistance() const {
    return min_distance_;
  }

} // namespace atom
