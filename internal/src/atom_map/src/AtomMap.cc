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

#include <atom_map/AtomMap.h>

#include <Eigen/Dense>

namespace gu = geometry_utils;
namespace pu = parameter_utils;

namespace atom {

  AtomMap::AtomMap() : initialized_(false) {}
  AtomMap::~AtomMap() {}

  // Initialize.
  bool AtomMap::Initialize(const ros::NodeHandle& n) {
    name_ = ros::names::append(n.getNamespace(), "atom_map");

    if (!LoadParameters(n)) {
      ROS_ERROR("%s: Failed to load parameters.", name_.c_str());
      return false;
    }

    if (!RegisterCallbacks(n)) {
      ROS_ERROR("%s: Failed to register callbacks.", name_.c_str());
      return false;
    }

    initialized_ = true;
    return true;
  }

  // Getters.
  double AtomMap::GetSignedDistance(double x, double y, double z) const {
    std::vector<Atom::Ptr> neighbors;
    if (!map_.GetKNearestNeighbors(x, y, z, num_neighbors_, neighbors)) {
      ROS_WARN("%s: Error in extracting nearest neighbors.", name_.c_str());
      return std::numeric_limits<double>::infinity();
    }

    // Generate covariance matrix for the training data.
    const size_t kNumNeighbors = neighbors.size();
    Eigen::Matrix
  }

  double AtomMap::GetProbability(double x, double y, double z) const;

  // Updates.
  void AtomMap::Update(const pcl::PointXYZ& point,
                       const pcl::PointXYZ& robot);
  void AtomMap::Update(const PointCloud::ConstPtr& cloud,
                       const pcl::PointXYZ& robot);

  // Load parameters and register callbacks.
  bool AtomMap::RegisterCallbacks(const ros::NodeHandle& n) { return true; }
  bool AtomMap::LoadParameters(const ros::NodeHandle& n) {
    if (!pu::Get("atom/gamma", gamma_)) return false;
    if (!pu::Get("atom/radius", radius_)) return false;
    if (!pu::Get("atom/thickness", max_surface_thickness_)) return false;

    return true;
  }

  // Sample a ray. Given a robot position and a measured point, discretize the
  // ray from sensor to observation and insert/update atoms along the way.
  void AtomMap::SampleRay(const pcl::PointXYZ& point,
                          const pcl::PointXYZ& robot) {}

  // Apply the covariance kernel function.
  double AtomMap::CovarianceKernel(const pcl::PointXYZ& p1,
                                   const pcl::PointXYZ& p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;
    return std::exp(-gamma_ * (dx*dx + dy*dy + dz*dz));
  }
}
