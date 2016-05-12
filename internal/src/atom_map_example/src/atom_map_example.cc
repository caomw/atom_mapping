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

namespace atom {
  AtomMapExample::AtomMapExample() : tf_listener_(tf_buffer_), initialized_(false) {}
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
    if (!pu::Get("atom_example/robot_frame", tf_robot_frame_)) return false;
    if (!pu::Get("atom_example/world_frame", tf_world_frame_)) return false;
    if (!pu::Get("atom_example/filter_leaf_size", filter_leaf_size_)) return false;

    return true;
  }

  bool AtomMapExample::RegisterCallbacks(const ros::NodeHandle& n) {
    ros::NodeHandle node(n);

    // Register point cloud subscriber callback.
    point_cloud_subscriber_ =
      node.subscribe<PointCloud>(data_topic_.c_str(), 10,
                                 &AtomMapExample::AddPointCloudCallback, this);

    return true;
  }

  // Callback to process point clouds.
  void AtomMapExample::AddPointCloudCallback(const PointCloud::ConstPtr& cloud) {
    // Get transform.
    geometry_msgs::TransformStamped tf;
    try {
      tf = tf_buffer_.lookupTransform(tf_world_frame_.c_str(),
                                      tf_robot_frame_.c_str(), ros::Time(0));
    } catch(tf2::TransformException &ex) {
      ROS_WARN("%s: %s", name_.c_str(), ex.what());
      ROS_WARN("%s: Did not insert this scan.", name_.c_str());
      ros::Duration(0.1).sleep();
      return;
    }

    // Voxel grid filter.
    PointCloud::Ptr filtered_cloud(new PointCloud);
    pcl::VoxelGrid<pcl::PointXYZ> grid_filter;
    grid_filter.setInputCloud(cloud);
    grid_filter.setLeafSize(filter_leaf_size_, filter_leaf_size_, filter_leaf_size_);
    grid_filter.filter(*filtered_cloud);

    // Transform point cloud into world frame.
    const gu::Transform3 pose = gr::FromROS(tf.transform);
    const Eigen::Matrix3d rotation = pose.rotation.Eigen();
    const Eigen::Vector3d translation = pose.translation.Eigen();

    Eigen::Matrix4d Rt = Eigen::Matrix4d::Identity();
    Rt.block(0, 0, 3, 3) = rotation;
    Rt.block(0, 3, 3, 1) = translation;

    PointCloud::Ptr transformed_cloud(new PointCloud);
    pcl::transformPointCloud(*filtered_cloud, *transformed_cloud, Rt);

    // Run map update.
    pcl::PointXYZ p(translation(0), translation(1), translation(2));
    map_.Update(transformed_cloud, p);

    // Publish.
    map_.PublishFull();
  }

} // namespace atom
