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

///////////////////////////////////////////////////////////////////////////////
//
// The OccupancyAtom class derives from the Atom class, and includes only
// signed distance estimation (and NOT occupancy).
//
///////////////////////////////////////////////////////////////////////////////

#ifndef ATOM_MAPPING_SDF_ATOM_H
#define ATOM_MAPPING_SDF_ATOM_H

#include <atom_map/Atom.h>

namespace atom {
  class SdfAtom : public Atom {
  public:
    ~SdfAtom();

    // Factory method.
    static Atom::Ptr Create(const Vector3f& p);

    // Getters.
    float GetSignedDistance() const;
    float GetSignedDistanceVariance() const;

    // Setters.
    void SetSignedDistance(float d);

    // Update the signed distance function for this atom. Variance is implicitly
    // the reciprocal of absolute distance, i.e. we trust measurements that say
    // the atom is closer to a surface more than those that say it is far away.
    // As above, weight is the overlap fraction between this atom and the one
    // which caused this update. Currently, we scale variance ~1/(0.5 + weight).
    // This scaling, however, is completely arbitrary.
    void UpdateSignedDistance(float sdf_update, float weight = 1.0);

  private:
    // Private constructor.
    SdfAtom(const Vector3f& p);

    // Signed distance estimate. By convention, this will be a positive number for
    // atoms that are in free space, and negative for those that are within obstacles.
    float sdf_mean_;

    // Uncertainty of the signed distance estimate. Do a simple maximum likelihood
    // update with each new measurement.
    float sdf_variance_;
  }; //\class SdfAtom
} //\namespace atom

#endif
