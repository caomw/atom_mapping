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

#include <atom_map/OccupancyAtom.h>

#include <glog/logging.h>
#include <iostream>
#include <math.h>

namespace atom {

  // Set the static log odds clamps to reasonable values.
  float OccupancyAtom::log_odds_clamp_low_ = ToLogOdds(0.03);
  float OccupancyAtom::log_odds_clamp_high_ = ToLogOdds(0.97);

  // Constructor/destructor.
  OccupancyAtom::~OccupancyAtom() {}
  OccupancyAtom::OccupancyAtom(const Vector3f& p)
    : Atom(p) {}

  // Factory method.
  Atom::Ptr OccupancyAtom::Create(const Vector3f& p) {
    Atom::Ptr ptr(new OccupancyAtom(p));
    return ptr;
  }

  // Getters.
  float OccupancyAtom::GetLogOdds() const { return log_odds_; }
  float OccupancyAtom::GetProbability() const {
    return ToProbability(log_odds_);
  }

  // Setters.
  void OccupancyAtom::SetProbabilityClamps(float low, float high) {
    SetLogOddsClamps(ToLogOdds(low), ToLogOdds(high));
  }

  void OccupancyAtom::SetLogOddsClamps(float low, float high) {
    log_odds_clamp_low_ = low;
    log_odds_clamp_high_ = high;
  }

  void OccupancyAtom::SetProbability(float p) {
    log_odds_ = ToLogOdds(p);
  }

  void OccupancyAtom::SetLogOdds(float l) {
    log_odds_ = l;
  }

  // Update the probability value stored in this atom. The weight parameter
  // is the fraction of overlap between two atoms. Currently, we use this
  // as a linear scaling between 0.5 (0) and the update probability (log-odds).
  void OccupancyAtom::UpdateProbability(float probability_update, float weight) {
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

    // Enforce clamping.
    if (log_odds_ > log_odds_clamp_high_) log_odds_ = log_odds_clamp_high_;
    else if (log_odds_ < log_odds_clamp_low_) log_odds_ = log_odds_clamp_low_;
  }

  void OccupancyAtom::UpdateLogOdds(float log_odds_update, float weight) {
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

    // Enforce clamping.
    if (log_odds_ > log_odds_clamp_high_) log_odds_ = log_odds_clamp_high_;
    else if (log_odds_ < log_odds_clamp_low_) log_odds_ = log_odds_clamp_low_;
  }
} //\namespace atom
