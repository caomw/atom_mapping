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

#include <atom_map_example/atom_map_example.h>
#include <parameter_utils/ParameterUtils.h>

namespace pu = parameter_utils;

namespace atom {
  AtomMapExample::AtomMapExample() {}
  AtomMapExample::~AtomMapExample() {}

  // Initialize.
  bool AtomMapExample::Initialize(const ros::NodeHandle& n) {
    name_ = ros::names::append(n.getNamespace(), "atom_map_example");

    if (!map_.Initialize(n)) {
      ROS_ERROR("%s: Failed to initialize AtomMap.", name_.c_str());
      return false;
    }

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

  // Load parameters and register callbacks.
  bool AtomMapExample::LoadParameters(const ros::NodeHandle& n) {
    if (!pu::Get("atom_example/data_topic", data_topic_)) return false;

    return true;
  }

  bool AtomMapExample::RegisterCallbacks(const ros::NodeHandle& n) {
    ros::NodeHandle node(n);

    // Register listener.
    point_cloud_subscriber_ =
      node.subscribe<PointCloud>(data_topic_.c_str(), 10,
                                 &AtomMapExample::AddPointCloudCallback, this);

    return true;
  }

  // Callback to process point clouds.
  void AtomMapExample::AddPointCloudCallback(const PointCloud::ConstPtr& cloud) {
    pcl::PointXYZ p(0.0, 0.0, 0.0);
    map_.Update(cloud, p);

    // Publish.
    map_.Publish();
  }

} // namespace atom
