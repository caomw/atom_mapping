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

#ifndef ATOM_MAPPING_ATOM_MAP_EXAMPLE_H
#define ATOM_MAPPING_ATOM_MAP_EXAMPLE_H

#include <atom_map/AtomMap.h>
#include <atom_map/AStarPlanner.h>

#include <pcl/point_types.h>
#include <pcl/common/transforms.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl_ros/point_cloud.h>
#include <tf2_ros/transform_listener.h>
#include <geometry_msgs/TransformStamped.h>
#include <geometry_msgs/PoseStamped.h>
#include <Eigen/Core>
#include <queue>

typedef pcl::PointCloud<pcl::PointXYZ> PointCloud;

using Eigen::Vector3f;

namespace atom {
  class AtomMapExample {
  public:
    AtomMapExample();
    ~AtomMapExample();

    // Initialize.
    bool Initialize(const ros::NodeHandle& n);

    // Optionally save.
    void MaybeSave() const;

    // Optionally compute an A* path from start point to finish point and publish.
    void MaybePublishPath();

  private:
    AtomMap map_;

    // Keep track of first and current positions.
    bool first_pose_;
    Vector3f initial_position_;
    Vector3f current_position_;

    // Subscribers and queues.
    ros::Subscriber point_cloud_subscriber_;
    ros::Subscriber pose_subscriber_;
    std::queue<const PointCloud::ConstPtr> point_cloud_queue_;
    std::queue<const Eigen::Matrix4d> pose_queue_;

    // Voxel grid filter leaf size. This should be on the order of the atomic
    // radius, in order to avoid oversampling.
    float filter_leaf_size_;

    // Republish filtered point clouds if there are any listeners.
    ros::Publisher filtered_cloud_publisher_;
    std::string filtered_cloud_topic_;

    // Publish A* path from initial position to final position. Planner can optionally
    // prefer paths near the specified manifold if path_sdf_manifold_ is positive.
    ros::Publisher path_publisher_;
    std::string path_topic_;
    float path_sdf_manifold_;

    // Topics to listen to. Pose topic is only used if buffer flag is set. If
    // buffer flag is not set, then we need to know the fixed frame id and use
    // ROS's transform tree.
    bool buffer_all_;
    std::string data_topic_;
    std::string pose_topic_;
    std::string fixed_frame_;
    tf2_ros::Buffer tf_buffer_;
    tf2_ros::TransformListener tf_listener_;

    // Optionally save the map on close.
    bool save_on_close_;
    std::string file_to_save_;

    // Name and initialization.
    bool initialized_;
    std::string name_;

    // Load parameters and register callbacks.
    bool LoadParameters(const ros::NodeHandle& n);
    bool RegisterCallbacks(const ros::NodeHandle& n);

    // Process a point cloud, pose pair.
    void ProcessPointCloud(const PointCloud::ConstPtr& cloud,
                           const Eigen::Matrix4d& pose);

    // Publishing.
    void PublishFilteredCloud(const PointCloud::ConstPtr& filtered);

    // Callbacks.
    void AddPointCloudCallback(const PointCloud::ConstPtr& cloud);
    void AddPoseCallback(const geometry_msgs::PoseStamped::ConstPtr& pose);
  }; // class AtomMapExample
} // namespace atom

#endif
