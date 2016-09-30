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

#ifndef ATOM_MAPPING_ATOM_HASH_GRID_H
#define ATOM_MAPPING_ATOM_HASH_GRID_H

#include <atom_map/Atom.h>
#include <atom_map/AtomIndex.h>

#include <pcl/point_types.h>
#include <unordered_map>

namespace atom {
  class AtomHashGrid {
  public:
    AtomHashGrid();
    ~AtomHashGrid();

    // Set atomic radius. This will throw an "initialized" flag, without which
    // nothing else will work.
    void SetAtomicRadius(float atomic_radius);

    // Nearest neighbor queries.
    bool RadiusSearch(float x, float y, float z, float r,
                      std::vector<Atom::Ptr>* neighbors);
    bool RadiusSearch(const pcl::PointXYZ& p, float r,
                      std::vector<Atom::Ptr>* neighbors);

    // Insert a new Atom.
    bool Insert(const Atom::Ptr& atom);

    // Return a list of all Atoms in the map.
    const std::vector<Atom::Ptr>& GetAtoms() const;

    // Return the number of Atoms in the map.
    size_t Size() const;

    // Return the maximum and minimum distances of any Atom to the surface.
    float GetMaxDistance() const;
    float GetMinDistance() const;

  private:
    // Populate a set of bins overlapping a particular ball.
    bool GetOverlappingBins(float x, float y, float z, float r,
                            std::vector<AtomIndex>* bins) const;

    // A hash map to store all Atoms, where the keys are grid indices on
    // an implicit grid of side length sufficiently small that no two Atoms
    // occupy the same grid cell, and the values are lists of indices into a
    // vector storing all Atoms.
    std::unordered_map<AtomIndex, std::vector<unsigned int>, AtomIndexHasher> map_;
    std::vector<Atom::Ptr> registry_;

    // Voxel side length -- slightly smaller than side length which separates
    // two Atoms on opposing corners of a cube.
    float voxel_size_;
    bool initialized_;

    // Keep track of maximum and minimum signed distances to the surface.
    float max_distance_;
    float min_distance_;
  };
}

#endif
