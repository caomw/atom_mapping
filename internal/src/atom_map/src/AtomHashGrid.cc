
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

using Eigen::Vector3f;

namespace atom {
  AtomHashGrid::~AtomHashGrid() {}
  AtomHashGrid::AtomHashGrid(double atomic_radius)
    : voxel_size_(1.9 * atomic_radius/sqrt(3.0)) {}

  // Nearest neighbor queries.
  bool AtomHashGrid::RadiusSearch(float x, float y, float z, float r,
                                  std::vector<Atom::Ptr>* neighbors);
  bool AtomHashGrid::RadiusSearch(const pcl::PointXYZ& p, float r,
                                  std::vector<Atom::Ptr>* neighbors);

  // Insert a new Atom.
  bool AtomHashGrid::Insert(Atom::Ptr atom) {
    const Vector3f p = atom->GetPosition();
    AtomIndex index(p(0), p(1), p(2), voxel_size_);

  }

  // Return a list of all Atoms in the map.
  const std::vector<Atom::Ptr>& AtomHashGrid::GetAtoms() const;

  // Return the number of Atoms in the map.
  size_t AtomHashGrid::Size() const;

  // Return the maximum and minimum distances of any Atom to the surface.
  float AtomHashGrid::GetMaxDistance() const;
  float AtomHashGrid::GetMinDistance() const;
}

#endif
