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

#include <atom_map/SdfAtom.h>

#include <glog/logging.h>
#include <iostream>
#include <math.h>

namespace atom {

  // Constructor/destructor.
  SdfAtom::~SdfAtom() {}
  SdfAtom::SdfAtom(const Vector3f& p)
    : Atom(p) {}

  // Factory method.
  Atom::Ptr SdfAtom::Create(const Vector3f& p) {
    Atom::Ptr ptr(new SdfAtom(p));
    return ptr;
  }

  // Getters.
  float SdfAtom::GetSignedDistance() const { return sdf_mean_; }
  float SdfAtom::GetSignedDistanceVariance() const { return sdf_variance_; }

  // Setters.
  void SdfAtom::SetSignedDistance(float d) {
    sdf_mean_ = d;
    sdf_variance_ = ToVariance(d);
  }

  // Update the signed distance function for this atom. Variance is implicitly
  // the reciprocal of absolute distance, i.e. we trust measurements that say
  // the atom is closer to a surface more than those that say it is far away.
  // As above, weight is the overlap fraction between this atom and the one
  // which caused this update. Currently, we scale variance ~1/(0.5 + weight).
  // This scaling, however, is completely arbitrary.
  void SdfAtom::UpdateSignedDistance(float sdf_update, float weight) {
#ifdef ENABLE_DEBUG_MESSAGES
    if (weight < 0.0 || weight > 1.0) {
      VLOG(1) << "Weight is not a number between [0, 1]: "
              << weight << ".";
    }
#endif
    const float sdf_variance_update = ToVariance(sdf_update) / (0.5 + weight);
    const float k = sdf_variance_ / (sdf_variance_ + sdf_variance_update);
    sdf_mean_ += k * (sdf_update - sdf_mean_);
    sdf_variance_ *= 1.0 - k;
  }
} //\namespace atom
