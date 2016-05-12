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

#ifndef ATOM_MAPPING_ATOM_MAP_H
#define ATOM_MAPPING_ATOM_MAP_H

#include <atom_map/Atom.h>

#include <flann/flann.h>
#include <pcl/point_types.h>
#include <glog/logging.h>
#include <math.h>

namespace gu = geometry_utils;

namespace atom {
  class AtomKdtree {
  public:
    AtomKdtree();
    ~AtomKdtree();

    // Nearest neighbor queries.
    bool GetKNearestNeighbors(double x, double y, double z, size_t k,
                              std::vector<Atom::Ptr>* neighbors);
    bool GetKNearestNeighbors(const pcl::PointXYZ& p, size_t k,
                              std::vector<Atom::Ptr>* neighbors);
    bool RadiusSearch(double x, double y, double z, double r,
                      std::vector<Atom::Ptr>* neighbors);
    bool RadiusSearch(const pcl::PointXYZ& p, double r,
                      std::vector<Atom::Ptr>* neighbors);

    // Insert a new Atom.
    bool Insert(Atom::Ptr atom);

    // Return a list of all Atoms in the map.
    const std::vector<Atom::Ptr>& GetAtoms() const;

    // Return the size of this tree.
    size_t Size() const;

  private:
    // A Flann kdtree to hold all the Atoms. Searches in this tree return
    // indices, which are then mapped to Atom::Ptr types in an array.
    std::shared_ptr< flann::Index< flann::L2<double> > > index_;
    std::vector<Atom::Ptr> registry_;

    // Find all neighbors for an Atom and set that Atom's neighbors_ field.
    bool SetNeighbors(Atom::Ptr atom);

    // Update neighbors' lists of neighboring Atoms to include a new Atom.
    void UpdateNeighbors(Atom::Ptr atom);
  };
}

#endif
