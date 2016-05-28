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
#include <atom_map/VoxelGrid.h>
#include <pcl/conversions.h>

namespace atom {
  AtomMapExample::AtomMapExample()
    : tf_listener_(tf_buffer_), initialized_(false) {}
  AtomMapExample::~AtomMapExample() {}

  // Optionally save on close.
  void AtomMapExample::Save() const {
    if (save_on_close_)
      map_.Save(file_to_save_);
  }

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
    if (!pu::Get("atom_example/filtered_cloud_topic", filtered_cloud_topic_))
      return false;
    if (!pu::Get("atom_example/filter_leaf_size", filter_leaf_size_)) return false;
    if (!pu::Get("atom_example/buffer_all", buffer_all_)) return false;
    if (!pu::Get("atom_example/pose_topic", pose_topic_)) return false;
    if (!pu::Get("atom_example/fixed_frame", fixed_frame_)) return false;
    if (!pu::Get("atom_example/save_on_close", save_on_close_)) return false;
    if (!pu::Get("atom_example/file_to_save", file_to_save_)) return false;

    return true;
  }

  bool AtomMapExample::RegisterCallbacks(const ros::NodeHandle& n) {
    ros::NodeHandle node(n);

    // Set buffer size. (Octomap uses 5.)
    const int kBufferSize = (buffer_all_) ? 100 : 5;

    // Register point cloud subscriber callback.
    point_cloud_subscriber_ =
      node.subscribe<PointCloud>(data_topic_.c_str(), kBufferSize,
                                 &AtomMapExample::AddPointCloudCallback, this);

    // Register robot pose subscriber callback.
    if (buffer_all_) {
      pose_subscriber_ =
        node.subscribe<geometry_msgs::PoseStamped>(pose_topic_.c_str(), kBufferSize,
                                                   &AtomMapExample::AddPoseCallback,
                                                   this);
    }

    // Set up filtered cloud publisher.
    filtered_cloud_publisher_ =
      node.advertise<PointCloud>(filtered_cloud_topic_.c_str(), 0);

    return true;
  }

  // Process a point cloud, pose pair.
  void AtomMapExample::ProcessPointCloud(const PointCloud::ConstPtr& cloud,
                                         const Eigen::Matrix4d& pose) {
    // Voxel grid filter.
    PointCloud::Ptr filtered_cloud(new PointCloud);
    VoxelGrid grid_filter(filter_leaf_size_);
    grid_filter.Filter(cloud, filtered_cloud);

    // Set frame name and publish.
    filtered_cloud->header = cloud->header;
    PublishFilteredCloud(filtered_cloud);

    // Transform point cloud into world frame.
    PointCloud::Ptr transformed_cloud(new PointCloud);
    pcl::transformPointCloud(*filtered_cloud, *transformed_cloud, pose);

    // Run map update.
    pcl::PointXYZ p(pose(0, 3), pose(1, 3), pose(2, 3));
    map_.Update(transformed_cloud, p);

    // Publish.
    map_.PublishFullOccupancy();
    map_.PublishFullSignedDistance();
  }

  // Publish the filtered cloud.
  void AtomMapExample::PublishFilteredCloud(const PointCloud::ConstPtr& filtered) {
    CHECK_NOTNULL(filtered.get());

    if (filtered_cloud_publisher_.getNumSubscribers() <= 0) return;
    filtered_cloud_publisher_.publish(filtered);
  }


  // Callback to process incoming pose messages.
  void AtomMapExample::AddPoseCallback(const geometry_msgs::PoseStamped::ConstPtr& pose) {
    // Extract pose as a 4x4 matrix in SE(3).
    const gu::Transform3 tf = gr::FromROS(pose->pose);
    const Eigen::Matrix3d rotation = tf.rotation.Eigen();
    const Eigen::Vector3d translation = tf.translation.Eigen();

    Eigen::Matrix4d Rt = Eigen::Matrix4d::Identity();
    Rt.block(0, 0, 3, 3) = rotation;
    Rt.block(0, 3, 3, 1) = translation;

    // Add to queue if there is no available point cloud.
    if (point_cloud_queue_.empty())
      pose_queue_.push(Rt);

    // Otherwise, go ahead and process this pair.
    else {
      const PointCloud::ConstPtr cloud = point_cloud_queue_.front();
      point_cloud_queue_.pop();
      ProcessPointCloud(cloud, Rt);
    }
  }

  // Callback to process incoming point clouds.
  void AtomMapExample::AddPointCloudCallback(const PointCloud::ConstPtr& cloud) {
    Eigen::Matrix4d pose = Eigen::Matrix4d::Identity();

    if (buffer_all_) {
      // Add to queue if there is no available pose.
      if (pose_queue_.empty())
        point_cloud_queue_.push(cloud);

      // Otherwise, go ahead and process this pair.
      else {
        pose = pose_queue_.front();
        pose_queue_.pop();
      }
    } else {
      // Get transform.
      geometry_msgs::TransformStamped tf;
      try {
        tf = tf_buffer_.lookupTransform(fixed_frame_.c_str(),
                                        cloud->header.frame_id,
                                        pcl_conversions::fromPCL(cloud->header.stamp));
      } catch(tf2::TransformException &ex) {
        ROS_WARN("%s: %s", name_.c_str(), ex.what());
        ROS_WARN("%s: Did not insert this scan.", name_.c_str());
        ros::Duration(0.1).sleep();
        return;
      }

      // Transform point cloud into world frame.
      const gu::Transform3 transform = gr::FromROS(tf.transform);
      const Eigen::Matrix3d rotation = transform.rotation.Eigen();
      const Eigen::Vector3d translation = transform.translation.Eigen();

      pose.block(0, 0, 3, 3) = rotation;
      pose.block(0, 3, 3, 1) = translation;
    }

    // Process the point cloud.
    ProcessPointCloud(cloud, pose);
  }

} // namespace atom
