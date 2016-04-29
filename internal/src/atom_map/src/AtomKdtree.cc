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

namespace atom {
  AtomKdtree::~AtomKdtree() { kd_free(tree); }
  AtomKdtree::AtomKdtree() {
    tree = kd_create(3); // 3D data.
    CHECK_NOTNULL(tree);
  };

  // Nearest neighbor queries.
  bool AtomKdtree::GetKNearestNeighbors(double x, double y, double z, double k,
                                        std::vector<Atom::Ptr>* neighbors) {
    // TODO!
    return true;
  }

  // Nearest neighbor queries.
  bool AtomKdtree::GetKNearestNeighbors(const pcl::PointXYZ& p, double k,
                                        std::vector<Atom::Ptr>* neighbors) {
    return GetKNearestNeighbors(p.x, p.y, p.z, k, neighbors);
  }

  // Radius searching.
  bool AtomKdtree::RadiusSearch(double x, double y, double z, double r,
                                std::vector<Atom::Ptr>* neighbors) {
    CHECK_NOTNULL(neighbors);
    neighbors->clear();

    // Run radius search.
    Kdsearch* result = kd_nearest_range3(tree, x, y, z, r);
    if (!result)
      return false;

    // Extract Atoms and add them to neighbors.
    while (!kd_res_end(result)) {
      Atom::Ptr atom = *(Atom::Ptr*) kd_res_item_data(result);
      neighbors->push_back(atom);

      kd_res_next(result);
    }

    kd_res_free(result);
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

    // Find all neighbors and add them to this Atom's neighbors list.
    if (!SetNeighbors(atom))
      return false;

    // Update those neighbors' lists to include this Atom.
    UpdateNeighbors(atom);

    // Insert.
    gu::Vec3 pos = atom->GetPosition();
    kd_insert3(tree, pos(0), pos(1), pos(2), (void *) &atom);
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
}
