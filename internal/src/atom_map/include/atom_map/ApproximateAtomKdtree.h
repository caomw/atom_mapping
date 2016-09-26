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

#ifndef ATOM_MAPPING_APPROX_ATOM_KDTREE_H
#define ATOM_MAPPING_APPROX_ATOM_KDTREE_H

#include <atom_map/Atom.h>

#include <flann/flann.h>
#include <pcl/point_types.h>
#include <glog/logging.h>
#include <math.h>

namespace atom {
  class ApproximateAtomKdtree {
  public:
    ApproximateAtomKdtree();
    ~ApproximateAtomKdtree();

    // Nearest neighbor queries.
    bool GetKNearestNeighbors(float x, float y, float z, size_t k,
                              std::vector<Atom::Ptr>* neighbors);
    bool GetKNearestNeighbors(const pcl::PointXYZ& p, size_t k,
                              std::vector<Atom::Ptr>* neighbors);
    bool RadiusSearch(float x, float y, float z, float r,
                      std::vector<Atom::Ptr>* neighbors);
    bool RadiusSearch(const pcl::PointXYZ& p, float r,
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
    std::shared_ptr< flann::KDTreeIndex< flann::L2<float> > > index_;
    std::vector<Atom::Ptr> registry_;
  }; // class ApproximateAtomKdtree
} // namespace atom

#endif
