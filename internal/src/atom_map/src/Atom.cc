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
#include <math.h>

namespace gu = geometry_utils;

namespace atom {

Atom::Atom() : log_odds_(ToLogOdds(0.5)), radius_(0.0) {}
Atom::~Atom() {}

double Atom::GetProbability() const {
  return ToProbability(log_odds_);
}

double Atom::GetLogOdds() const {
  return log_odds_;
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

void Atom::SetRadius(double r) {
  radius_ = r;
}

void Atom::SetPosition(const gu::Vec3& p) {
  position_ = p;
}

void Atom::UpdateProbability(double probability_update) {
#ifdef ENABLE_DEBUG_MESSAGES
  if (probability_update < 0.0 || probability_update > 1.0) {
    VLOG(1) << "Probability update is not a probability in [0, 1]: "
            << probability_update << ".";
  }
#endif
  log_odds_ += ToLogOdds(probability_update);
}

void Atom::UpdateLogOdds(double log_odds_update) {
#ifdef ENABLE_DEBUG_MESSAGES
  if (log_odds_update < 0.0) {
    VLOG(1) << "Log-odds update is less than zero: " << log_odds_update << ".";
  }
#endif
  log_odds_ += log_odds_update;
}

void Atom::AddNeighbor(Atom* neighbor) {
  CHECK_NOTNULL(neighbor);
  neighbors_.push_back(std::shared_ptr<Atom>(neighbor));
}

const std::vector< std::shared_ptr<Atom> >& Atom::GetNeighbors() const {
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

} //\namespace atom
