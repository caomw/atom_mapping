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

///////////////////////////////////////////////////////////////////////////////
//
// The Atom class defines a single atom in the map. These can be though of as
// small bubbles of a specified radius centered at a 2D or 3D point, with a
// probability of being occupied. Additional information can be added, such as
// signed distance to nearest surface. Atoms are the basic building blocks of
// atom maps.
//
///////////////////////////////////////////////////////////////////////////////

#include <geometry_utils/Vector3.h>

#include <pcl/point_types.h>
#include <memory>
#include <vector>

#ifndef ATOM_H
#define ATOM_H

namespace atom {
class Atom {
 public:
  ~Atom();

  typedef std::shared_ptr<Atom> Ptr;

  // Factory method.
  static Ptr Create();

  // Getters.
  double GetProbability() const;
  double GetLogOdds() const;
  double GetSignedDistance() const;
  double GetSignedDistanceVariance() const;
  double GetRadius() const;
  geometry_utils::Vec3 GetPosition() const;

  // Setters.
  void SetProbability(double p);
  void SetLogOdds(double l);
  void SetSignedDistance(double d);
  void SetRadius(double r);
  void SetPosition(const geometry_utils::Vec3& p);

  // Update the probability value stored in this atom.
  void UpdateProbability(double probability_update);
  void UpdateLogOdds(double log_odds_update);

  // Update the signed distance function for this atom. Variance is implicitly
  // the reciprocal of absolute distance, i.e. we trust measurements that say
  // the atom is closer to a surface more than those that say it is far away.
  void UpdateSignedDistance(double sdf_update);

  // Check if this Atom contains a point.
  bool Contains(double x, double y, double z);
  bool Contains(pcl::PointXYZ& p);

  // Check distance to a point.
  double GetDistanceTo(double x, double y, double z);
  double GetDistanceTo(pcl::PointXYZ& p);

  // Add a neighboring atom to this one.
  void AddNeighbor(Atom::Ptr neighbor);
  const std::vector<Ptr>& GetNeighbors() const;

 private:
  // Log-odds probability that this chunk of space is occupied. We use log-odds
  // to avoid numerical instability when multiplying, e.g., 1e-4 to itself
  // multiple times.
  double log_odds_;

  // Signed distance estimate. By convention, this will be a positive number for
  // atoms that are in free space, and negative for those that are within obstacles.
  double sdf_mean_;

  // Uncertainty of the signed distance estimate. Do a simple maximum likelihood
  // update with each new measurement.
  double sdf_variance_;

  // Position in 3D space.
  geometry_utils::Vec3 position_;

  // Atomic radius. Other atoms cannot be inserted into the map within this
  // radius.
  double radius_;

  // Pointers to neighboring atoms. This list is incrementally updated when new
  // atoms are added to the map, and begins empty.
  std::vector<Ptr> neighbors_;

  // Private constructor.
  Atom();
}; //\class Atom


// Conversion from a probability in [0, 1] to a log-odds probability in [0,
// infty).
double ToProbability(double log_odds);

// Conversion from a log-odds probability in [0, infty) to a probability in [0,
// 1].
double ToLogOdds(double probability);

// An arbitrary map from signed distance values to variances. The general idea is
// that smaller values (by magnitude) are more reliable because there is a smaller
// chance that there exists a nearer surface.
double ToVariance(double sdf_update);

} //\namespace atom

#endif
