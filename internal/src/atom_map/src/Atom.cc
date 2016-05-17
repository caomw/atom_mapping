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

namespace gu = geometry_utils;

namespace atom {

  Atom::~Atom() {}
  Atom::Atom(double r)
    : log_odds_(ToLogOdds(0.5)), radius_(r), sdf_mean_(0.0),
      sdf_variance_(std::numeric_limits<double>::infinity()),
      position_(0.0, 0.0, 0.0) {}

  Atom::Ptr Atom::Create(double r) {
    Atom::Ptr atom(new Atom(r));
    return atom;
  }

  double Atom::GetProbability() const {
    return ToProbability(log_odds_);
  }

  double Atom::GetLogOdds() const {
    return log_odds_;
  }

  double Atom::GetSignedDistance() const {
    return sdf_mean_;
  }

  double Atom::GetSignedDistanceVariance() const {
    return sdf_variance_;
  }

  double Atom::GetRadius() const {
    return radius_;
  }

  gu::Vec3 Atom::GetPosition() const {
    return position_;
  }

  void Atom::SetProbability(double p) {
    log_odds_ = ToLogOdds(p);
  }

  void Atom::SetLogOdds(double l) {
    log_odds_ = l;
  }

  void Atom::SetSignedDistance(double d) {
    sdf_mean_ = d;
    sdf_variance_ = ToVariance(d);
  }

  void Atom::SetPosition(const gu::Vec3& p) {
    position_ = p;
  }

  void Atom::UpdateProbability(double probability_update, double weight) {
#ifdef ENABLE_DEBUG_MESSAGES
    if (probability_update < 0.0 || probability_update > 1.0) {
      VLOG(1) << "Probability update is not a probability in [0, 1]: "
              << probability_update << ".";
    }
    if (weight < 0.0 || weight > 1.0) {
      VLOG(1) << "Weight is not a number between [0, 1]: "
              << weight << ".";
    }
#endif
    log_odds_ += weight * ToLogOdds(probability_update);
  }

  void Atom::UpdateLogOdds(double log_odds_update, double weight) {
#ifdef ENABLE_DEBUG_MESSAGES
    if (log_odds_update < 0.0) {
      VLOG(1) << "Log-odds update is less than zero: " << log_odds_update << ".";
    }
    if (weight < 0.0 || weight > 1.0) {
      VLOG(1) << "Weight is not a number between [0, 1]: "
              << weight << ".";
    }
#endif
    log_odds_ += weight * log_odds_update;
  }

  void Atom::UpdateSignedDistance(double sdf_update, double weight) {
#ifdef ENABLE_DEBUG_MESSAGES
    if (weight < 0.0 || weight > 1.0) {
      VLOG(1) << "Weight is not a number between [0, 1]: "
              << weight << ".";
    }
#endif
    const double sdf_variance_update = ToVariance(sdf_update) / (0.5 + weight);
    const double k = sdf_variance_ / (sdf_variance_ + sdf_variance_update);
    sdf_mean_ += k * (sdf_update - sdf_mean_);
    sdf_variance_ *= 1.0 - k;
  }

  bool Atom::Contains(double x, double y, double z) const {
    return GetDistanceTo(x, y, z) <= radius_;
  }

  bool Atom::Contains(const Atom::Ptr& atom) const {
    return GetDistanceTo(atom) <= radius_;
  }

  bool Atom::Contains(const pcl::PointXYZ& p) const {
    return GetDistanceTo(p) <= radius_;
  }

  double Atom::GetDistanceTo(double x, double y, double z) const {
    const double dx = x - position_(0);
    const double dy = y - position_(1);
    const double dz = z - position_(2);
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }

  double Atom::GetDistanceTo(const Atom::Ptr& atom) const {
    CHECK_NOTNULL(atom.get());
    gu::Vec3 p = atom->GetPosition();

    const double dx = p(0) - position_(0);
    const double dy = p(1) - position_(1);
    const double dz = p(2) - position_(2);
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }

  double Atom::GetDistanceTo(const pcl::PointXYZ& p) const {
    const double dx = p.x - position_(0);
    const double dy = p.y - position_(1);
    const double dz = p.z - position_(2);
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }

  double Atom::ComputeOverlapFraction(const Atom::Ptr& atom) const {
    const double distance = GetDistanceTo(atom);

    return 1 - (0.25 / pow(radius_, 3)) * distance *
      (3.0 * pow(radius_, 2) - 0.25 * pow(distance, 2));
  }

  void Atom::AddNeighbor(Atom::Ptr neighbor) {
    CHECK_NOTNULL(neighbor.get());
    neighbors_.push_back(neighbor);
  }

  void Atom::ClearNeighbors() {
    neighbors_.clear();
  }

  const std::vector<Atom::Ptr>& Atom::GetNeighbors() const {
    return neighbors_;
  }

  double ToProbability(double log_odds) {
#ifdef ENABLE_DEBUG_MESSAGES
    if (log_odds < 0.0) {
      VLOG(1) << "Input log-odds value is less than zero: " << log_odds << ".";
    }
#endif
    return 1.0 - (1.0 / (1.0 + exp(log_odds)));
  }

  double ToLogOdds(double probability) {
#ifdef ENABLE_DEBUG_MESSAGES
    if (probability < 0.0 || probability > 1.0) {
      VLOG(1) << "Input Probability value is not in [0, 1]: " << probability
              << ".";
    }
#endif
    return log(probability / (1.0 - probability));
  }

  double ToVariance(double sdf_update) {
    return sdf_update * sdf_update; // This is arbitrary.
  }

} //\namespace atom
