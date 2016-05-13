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

namespace atom {
  AtomKdtree::AtomKdtree() :
    max_distance_(-std::numeric_limits<double>::infinity()),
    min_distance_(std::numeric_limits<double>::infinity()) {}
  AtomKdtree::~AtomKdtree() {
    // Free memory from points in the kd tree.
    if (index_ != nullptr) {
      for (size_t ii = 0; ii < index_->size(); ++ii) {
        double* point = index_->getPoint(ii);
        delete[] point;
      }
    }
  }

  // Nearest neighbor queries.
  bool AtomKdtree::GetKNearestNeighbors(double x, double y, double z, size_t k,
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
    flann::Matrix<double> flann_query(new double[kNumColumns], 1, kNumColumns);
    flann_query[0][0] = x;
    flann_query[0][1] = y;
    flann_query[0][2] = z;

    // Search the kd tree for the nearest neighbor to the query.
    std::vector< std::vector<int> > query_match_indices;
    std::vector< std::vector<double> > query_distances;

    int num_neighbors_found =
      index_->knnSearch(flann_query, query_match_indices,
                        query_distances, static_cast<int>(k),
                        flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));

    // Assign output.
    for (size_t ii = 0; ii < num_neighbors_found; ii++)
      neighbors->push_back(registry_[ query_match_indices[0][ii] ]);

    return true;
  }

  // Nearest neighbor queries.
  bool AtomKdtree::GetKNearestNeighbors(const pcl::PointXYZ& p, size_t k,
                                        std::vector<Atom::Ptr>* neighbors) {
    return GetKNearestNeighbors(p.x, p.y, p.z, k, neighbors);
  }

  // Radius searching.
  bool AtomKdtree::RadiusSearch(double x, double y, double z, double r,
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
    flann::Matrix<double> flann_query(new double[kNumColumns], 1, kNumColumns);
    flann_query[0][0] = x;
    flann_query[0][1] = y;
    flann_query[0][2] = z;

    // Search the kd tree for the nearest neighbor to the query.
    std::vector< std::vector<int> > query_match_indices;
    std::vector< std::vector<double> > query_distances;

    int num_neighbors_found =
      index_->radiusSearch(flann_query, query_match_indices,
                           query_distances, static_cast<float>(r),
                           flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));

    // Assign output.
    for (size_t ii = 0; ii < num_neighbors_found; ii++)
      neighbors->push_back(registry_[ query_match_indices[0][ii] ]);

    return true;
  }

  // Radius searching.
  bool AtomKdtree::RadiusSearch(const pcl::PointXYZ& p, double r,
                                std::vector<Atom::Ptr>* neighbors) {
    return RadiusSearch(p.x, p.y, p.z, r, neighbors);
  }

  // Insert a new Atom.
  bool AtomKdtree::Insert(Atom::Ptr atom) {
    CHECK_NOTNULL(atom.get());

    // Copy the input point into FLANN's Matrix type.
    const int kNumColumns = 3;
    flann::Matrix<double> flann_point(new double[kNumColumns], 1, kNumColumns);
    gu::Vec3 pos = atom->GetPosition();
    flann_point[0][0] = pos(0);
    flann_point[0][1] = pos(1);
    flann_point[0][2] = pos(2);

    // If this is the first point in the index, create the index and exit.
    if (index_ == nullptr) {
      // Single kd-tree. No approximation.
      const int kNumRandomizedKDTrees = 1;
      index_.reset(new flann::Index< flann::L2<double> >(
                   flann_point, flann::KDTreeIndexParams(kNumRandomizedKDTrees)));
      index_->buildIndex();
    } else {
      // Find all neighbors and add them to this Atom's neighbors list.
      if (!SetNeighbors(atom))
        return false;

      // Update those neighbors' lists to include this Atom.
      UpdateNeighbors(atom);

      // If the index is already created, add the data point to the index. Rebuild
      // every time the index doubles in size to occasionally rebalance the kdtree.
      const int kRebuildThreshold = 2;
      index_->addPoints(flann_point, kRebuildThreshold);
    }

    // Add point to registry.
    registry_.push_back(atom);

    // Update max and min distances.
    const double sdf = atom->GetSignedDistance();
    if (sdf > max_distance_) {
      max_distance_ = sdf;
    } else if (sdf < min_distance_) {
      min_distance_ = sdf;
    }

    return true;
  }

  // Find all neighbors for an Atom and set that Atom's neighbors_ field.
  bool AtomKdtree::SetNeighbors(Atom::Ptr atom) {
    CHECK_NOTNULL(atom.get());

    // Find all neighbors within a radius of one atomic diameter.
    std::vector<Atom::Ptr> neighbors;
    gu::Vec3 pos = atom->GetPosition();
    if (!RadiusSearch(pos(0), pos(1), pos(2), atom->GetRadius(), &neighbors))
      return false;

    // Add all neighbors. Check if any are within one atomic radius. If so,
    // clear neighbors and return false.
    atom->ClearNeighbors();
    for (size_t ii = 0; ii < neighbors.size(); ii++) {
      if (neighbors[ii]->Contains(atom)) {
        atom->ClearNeighbors();
        return false;
      }

      atom->AddNeighbor(neighbors[ii]);
    }

    return true;
  }

  // Update neighbors' lists of neighboring Atoms to include a new Atom.
  void AtomKdtree::UpdateNeighbors(Atom::Ptr atom) {
    CHECK_NOTNULL(atom.get());

    // Add this atom to each neighbor.
    std::vector<Atom::Ptr> neighbors = atom->GetNeighbors();
    for (size_t ii = 0; ii < neighbors.size(); ii++)
      neighbors[ii]->AddNeighbor(atom);
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
  double AtomKdtree::GetMaxDistance() const {
    return max_distance_;
  }

  double AtomKdtree::GetMinDistance() const {
    return min_distance_;
  }

} // namespace atom
