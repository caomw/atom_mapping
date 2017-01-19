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
 * Authors: Erik Nelson            ( eanelson@eecs.berkeley.edu )
 *          David Fridovich-Keil   ( dfk@eecs.berkeley.edu )
 */

#include <atom_map/Atom.h>

#include <glog/logging.h>
#include <iostream>
#include <math.h>

namespace atom {

  // Set static radius to something sensible.
  float Atom::radius_ = 0.5;

  // Protected constructory/destructor.
  Atom::~Atom() {}
  Atom::Atom(const Vector3f& p)
    : position_(p) {}

  // Getters.
  float Atom::GetRadius() const { return radius_; }
  Vector3f& Atom::GetPosition() { return position_; }

  // Static setter.
  void Atom::SetRadius(float r) { radius_ = r; }

  // Check if this Atom contains a point.
  bool Atom::Contains(float x, float y, float z) const {
    return GetDistanceTo(x, y, z) <= radius_;
  }

  bool Atom::Contains(const Atom::Ptr& atom) const {
    return GetDistanceTo(atom) <= radius_;
  }

  bool Atom::Contains(const pcl::PointXYZ& p) const {
    return GetDistanceTo(p) <= radius_;
  }

  // Compute the distance to a point.
  float Atom::GetDistanceTo(float x, float y, float z) const {
    const float dx = x - position_(0);
    const float dy = y - position_(1);
    const float dz = z - position_(2);
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }

  float Atom::GetDistanceTo(const Atom::Ptr& atom) const {
    CHECK_NOTNULL(atom.get());
    Vector3f p = atom->GetPosition();

    const float dx = p(0) - position_(0);
    const float dy = p(1) - position_(1);
    const float dz = p(2) - position_(2);
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }

  float Atom::GetDistanceTo(const pcl::PointXYZ& p) const {
    const float dx = p.x - position_(0);
    const float dy = p.y - position_(1);
    const float dz = p.z - position_(2);
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }

  float Atom::ComputeOverlapFraction(const Atom::Ptr& atom) const {
    const float distance = GetDistanceTo(atom);

    return 1 - (0.25 / pow(radius_, 3)) * distance *
      (3.0 * pow(radius_, 2) - 0.25 * pow(distance, 2));
  }

  float Atom::ComputeOverlapFraction(const pcl::PointXYZ& p) const {
    const float distance = GetDistanceTo(p);

    return 1 - (0.25 / pow(radius_, 3)) * distance *
      (3.0 * pow(radius_, 2) - 0.25 * pow(distance, 2));
  }

  float Atom::ComputeOverlapFraction(float x, float y, float z) const {
    const float distance = GetDistanceTo(x, y, z);

    return 1 - (0.25 / pow(radius_, 3)) * distance *
      (3.0 * pow(radius_, 2) - 0.25 * pow(distance, 2));
  }

  // Helpers to convert between log odds and probability.
  float ToProbability(float log_odds) {
#ifdef ENABLE_DEBUG_MESSAGES
    if (log_odds < 0.0) {
      VLOG(1) << "Input log-odds value is less than zero: " << log_odds << ".";
    }
#endif
    return 1.0 - (1.0 / (1.0 + exp(log_odds)));
  }

  float ToLogOdds(float probability) {
#ifdef ENABLE_DEBUG_MESSAGES
    if (probability < 0.0 || probability > 1.0) {
      VLOG(1) << "Input Probability value is not in [0, 1]: " << probability
              << ".";
    }
#endif
    return log(probability / (1.0 - probability));
  }

  // Helper to compute 'variance' of an SDF update.
  float ToVariance(float sdf_update) {
    return sdf_update * sdf_update; // This is arbitrary.
  }

} //\namespace atom
