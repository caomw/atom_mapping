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
// Test script for the Atom and AtomKdtree classes.
//
///////////////////////////////////////////////////////////////////////////////

#include <gtest/gtest.h>

#include <atom_map/Atom.h>
#include <atom_map/OccupancyAtom.h>
#include <atom_map/SdfAtom.h>
#include <atom_map/VoxelGrid.h>
#include <atom_map/AtomKdtree.h>

#include <Eigen/Core>
#include <math.h>
#include <random>
#include <iostream>

using namespace atom;
using Eigen::Vector3f;

// Test the VoxelGrid class.
TEST(VoxelGrid, TestVoxelGrid) {
  PointCloud::Ptr raw(new PointCloud);
  PointCloud::Ptr filtered(new PointCloud);

  // Put a bunch of random points in the raw point cloud.
  const size_t kNumPoints = 1000;
  const float kLowerBound = -1000.0;
  const float kUpperBound = 1000.0;
  std::uniform_real_distribution<float> unif(kLowerBound, kUpperBound);
  std::default_random_engine rng;

  for (size_t ii = 0; ii < kNumPoints; ii++) {
    raw->points.push_back(pcl::PointXYZ(unif(rng), unif(rng), unif(rng)));
  }

  // Run the VoxelGrid filter.
  const float kLeafSize = 0.05;
  const float kMaxSize = 10000.0;
  VoxelGrid vg(kLeafSize, -kMaxSize, kMaxSize);
  EXPECT_TRUE(vg.Filter(raw, filtered));
}

// Test that we can update probability/log odds.
TEST(OccupancyAtom, TestUpdateProbability) {
  // Set params.
  const float kProbability = 0.5;
  const float kClampLow = 0.01;
  const float kClampHigh = 0.99;
  const float kSignedDistance = 1.0;
  const float kRadius = 1.0;
  const Vector3f kPosition(1.0, 1.0, 1.0);

  Atom::SetRadius(kRadius);
  OccupancyAtom::SetProbabilityClamps(kClampLow, kClampHigh);
  Atom::Ptr atom = OccupancyAtom::Create(kPosition);
  CHECK_NOTNULL(atom.get());

  std::static_pointer_cast<OccupancyAtom>(atom)->SetProbability(kProbability);

  // Update probability a bunch of times with the same small value and
  // make sure that the resulting probability is approximately zero.
  const size_t kNumUpdates = 30;
  const float kProbabilityUpdate = 0.1;
  for (size_t ii = 0; ii < kNumUpdates; ii++)
    std::static_pointer_cast<OccupancyAtom>(atom)->
      UpdateProbability(kProbabilityUpdate);

  EXPECT_NEAR(std::static_pointer_cast<OccupancyAtom>(atom)->GetProbability(),
              kClampLow, 1e-4);
}

// Test that we can update signed distance.
TEST(SdfAtom, TestUpdateSignedDistance) {
  // Set params.
  const float kRadius = 0.1;
  const float kProbability = 0.5;
  const float kSignedDistance = 1.0;
  const Vector3f kPosition(1.0, 1.0, 1.0);

  Atom::SetRadius(kRadius);
  Atom::Ptr atom = SdfAtom::Create(kPosition);
  CHECK_NOTNULL(atom.get());

  std::static_pointer_cast<SdfAtom>(atom)->SetSignedDistance(kSignedDistance);

  // Update signed distance a bunch of times with the same value and make sure
  // that the resulting signed distance is approximately the new value.
  const size_t kNumUpdates = 100;
  const float kSignedDistanceUpdate = 0.3;
  for (size_t ii = 0; ii < kNumUpdates; ii++)
    std::static_pointer_cast<SdfAtom>(atom)->
      UpdateSignedDistance(kSignedDistanceUpdate);

  EXPECT_NEAR(std::static_pointer_cast<SdfAtom>(atom)->GetSignedDistance(),
              kSignedDistanceUpdate, 1e-3);
}

// Test that we can create an AtomKdtree and insert lots of points.
TEST(AtomKdtree, TestAtomKdtreeInsertion) {
  AtomKdtree tree;

  // Set params for random point generation.
  const float kRadius = 0.0001;
  Atom::SetRadius(kRadius);
  const size_t kNumPoints = 500;
  const float kLowerBound = 0.0;
  const float kUpperBound = 1.0;
  std::random_device rd;
  std::default_random_engine rng(rd());
  std::uniform_real_distribution<float> unif(kLowerBound, kUpperBound);

  // Generate a bunch of random points and add to the tree.
  std::vector<Atom::Ptr> atoms;
  for (size_t ii = 0; ii < kNumPoints; ii++) {
    const Vector3f pos(unif(rng), unif(rng), unif(rng));
    Atom::Ptr atom = OccupancyAtom::Create(pos);

    // Insert.
    ASSERT_TRUE(tree.Insert(atom));
    atoms.push_back(atom);
  }

  // Check nearest neighbors.
  for (size_t ii = 0; ii < atoms.size(); ii++) {
    Atom::Ptr atom = atoms[ii];

    // Nearest neighbor search.
    std::vector<Atom::Ptr> neighbors;
    const size_t kNearestNeighbors = 1;
    ASSERT_TRUE(tree.GetKNearestNeighbors(atom->GetPosition()(0),
                                          atom->GetPosition()(1),
                                          atom->GetPosition()(2),
                                          kNearestNeighbors, &neighbors));

    ASSERT_EQ(neighbors.size(), 1);

    // Check that the nearest neighbor matches.
    EXPECT_EQ(neighbors[0]->GetPosition()(0), atom->GetPosition()(0));
    EXPECT_EQ(neighbors[0]->GetPosition()(1), atom->GetPosition()(1));
    EXPECT_EQ(neighbors[0]->GetPosition()(2), atom->GetPosition()(2));
  }
}

// Test radius search by making sure all points returned are actually within
// the specified radius.
TEST(AtomKdtree, TestAtomKdtreeRadiusSearch) {
  AtomKdtree tree;

  // Set params for random point generation.
  const float kRadius = 0.01;
  Atom::SetRadius(kRadius);
  const size_t kNumPoints = 1000;
  const float kLowerBound = -10.0;
  const float kUpperBound = 10.0;
  std::random_device rd;
  std::default_random_engine rng(rd());
  std::uniform_real_distribution<float> unif(kLowerBound, kUpperBound);

  // Generate a bunch of random points and add to the tree.
  for (size_t ii = 0; ii < kNumPoints; ii++) {
    const Vector3f pos(unif(rng), unif(rng), unif(rng));
    Atom::Ptr atom = OccupancyAtom::Create(pos);

    // Insert.
    ASSERT_TRUE(tree.Insert(atom));
  }

  // Check nearest neighbors with a radius search.
  const size_t kNumChecks = 1000;
  for (size_t ii = 0; ii < kNumChecks; ii++) {
    std::vector<Atom::Ptr> neighbors;
    const float kSearchRadius = 2.0 * kRadius;
    const pcl::PointXYZ p(unif(rng), unif(rng), unif(rng));
    ASSERT_TRUE(tree.RadiusSearch(p.x, p.y, p.z,
                                  kSearchRadius, &neighbors));

    for (size_t jj = 0; jj < neighbors.size(); jj++)
      EXPECT_TRUE(neighbors[jj]->GetDistanceTo(p) < kSearchRadius);
  }
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
