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
// occupancy probability estimation (and NOT signed distance function).
//
///////////////////////////////////////////////////////////////////////////////

#ifndef ATOM_MAPPING_OCCUPANCY_ATOM_H
#define ATOM_MAPPING_OCCUPANCY_ATOM_H

#include <atom_map/Atom.h>

namespace atom {
  class OccupancyAtom : public Atom {
  public:
    ~OccupancyAtom();

    // Factory method.
    static Atom::Ptr Create(const Vector3f& p);

    // Getters.
    float GetProbability() const;
    float GetLogOdds() const;

    // Setters.
    static void SetProbabilityClamps(float low, float high);
    static void SetLogOddsClamps(float low, float high);
    void SetProbability(float p);
    void SetLogOdds(float l);

    // Update the probability value stored in this atom. The weight parameter
    // is the fraction of overlap between two atoms. Currently, we use this
    // as a linear scaling between 0.5 (0) and the update probability (log-odds).
    void UpdateProbability(float probability_update, float weight = 1.0);
    void UpdateLogOdds(float log_odds_update, float weight = 1.0);

  private:
    // Private constructor.
    OccupancyAtom(const Vector3f& p);

    // Log-odds probability that this chunk of space is occupied. We use
    // log-odds to avoid numerical instability when updating.
    float log_odds_;

    // Clamping thresholds. Above or below these, do not update probability.
    static float log_odds_clamp_low_;
    static float log_odds_clamp_high_;
  }; //\class OccupancyAtom

} //\namespace atom

#endif
