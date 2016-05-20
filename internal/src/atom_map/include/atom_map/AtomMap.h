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

#ifndef ATOM_MAPPING_ATOM_MAP_H
#define ATOM_MAPPING_ATOM_MAP_H

#include <atom_map/Atom.h>
#include <atom_map/AtomKdtree.h>
#include <atom_map/RaySamples.h>
#include <parameter_utils/ParameterUtils.h>
#include <geometry_utils/GeometryUtilsROS.h>

#include <ros/ros.h>
#include <std_msgs/ColorRGBA.h>
#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>
#include <pcl_ros/point_cloud.h>
#include <glog/logging.h>
#include <math.h>

typedef pcl::PointCloud<pcl::PointXYZ> PointCloud;

namespace atom {
  class AtomMap {
  public:
    AtomMap();
    ~AtomMap();

    // Initialize.
    bool Initialize(const ros::NodeHandle& n);

    // Getters.
    void GetSignedDistance(float x, float y, float z,
                           float* distance, float* variance);
    float GetProbability(float x, float y, float z);

    // Update.
    void Update(const PointCloud::ConstPtr& cloud, const pcl::PointXYZ& robot);

    // Publishing.
    void PublishFullOccupancy() const;
    void PublishFullSignedDistance() const;

  private:
    // A kdtree to hold all the Atoms.
    AtomKdtree map_;

    // Atomic radius.
    float radius_;

    // Min and max scan ranges.
    float min_scan_range_;
    float max_scan_range_;

    // Optionally update signed distance and/or occupancy.
    bool update_occupancy_;
    bool update_signed_distance_;

    // Radius for nearest neighbor search for surface normal extraction.
    float surface_normal_radius_;

    // Maximum distance to trace a ray inside of an object. This should be small,
    // in order to ensure we don't accidentally trace all the way through. However,
    // it is definitely necessary in order to ensure proper surface detection.
    float max_occupied_backoff_;

    // Maximum distance to trace the surface normal inside (supposedly) free space.
    float max_normal_backoff_;

    // Maximum number of Atoms to lay down along the normal vector into free space.
    int max_samples_normal_;

    // Sensor angular resolution. This is the smallest anglue between two range
    // measurements. We also provide a flag to turn on angular interleaving.
    float angular_resolution_;
    bool angular_interleaving_;
    float lambda_;

    // Probability of hits and misses for occupancy updates.
    float probability_hit_;
    float probability_miss_;

    // Probability clamping values.
    float probability_clamp_high_;
    float probability_clamp_low_;

    // Number of nearest neighbors to examine for GP surface distance regression.
    int num_neighbors_;

    // Characteristic parameter for the radial basis function used as a
    // covariance kernel for signed distance function estimation. Larger gamma
    // means that covariance decreases faster as points get farther apart.
    float gamma_;

    // Noise variance to add to covariance matrices. Assume isotropic noise.
    float noise_variance_;

    // Visualization parameters and publishers.
    std::string fixed_frame_id_;
    std::string full_occupancy_topic_;
    std::string full_sdf_topic_;
    bool only_show_occupied_;
    float occupied_threshold_; // Above this is considered occupied.
    float sdf_threshold_;      // Smaller than this is considered occupied.
    ros::Publisher full_occupancy_publisher_;
    ros::Publisher full_sdf_publisher_;

    // Initialization and naming.
    bool initialized_;
    std::string name_;

    // Load parameters and register callbacks.
    bool LoadParameters(const ros::NodeHandle& n);
    bool RegisterCallbacks(const ros::NodeHandle& n);

    // Try to add a single atom.
    void MaybeInsertAtom(const pcl::PointXYZ& position, float sdf);
    void MaybeInsertAtom(const Atom::Ptr& atom);

    // Copy parameters from a different AtomMap.
    void CopyParametersFrom(const AtomMap& reference);

    // Return a list of all Atoms in the map.
    const std::vector<Atom::Ptr>& GetAtoms() const;

    // Sample a ray and do a probabilistic and signed distance update.
    // Given a robot position and a measured point, discretize the ray from sensor
    // to observation and return vectors of points and distances. Append to the
    // output RaySamples argument.
    void SampleRay(const pcl::PointXYZ& point, const pcl::Normal& normal,
                   const pcl::PointXYZ& robot, RaySamples* samples);
#if 0
    void Update(const pcl::PointXYZ& point, const pcl::Normal& normal,
                const pcl::PointXYZ& robot);
#endif

    // Apply the covariance kernel function.
    float CovarianceKernel(const pcl::PointXYZ& p1, const pcl::PointXYZ& p2);

    // Convert a probability of occupancy to a ROS color.
    std_msgs::ColorRGBA ProbabilityToRosColor(float probability) const;

    // Convert a signed distnace value to a ROS color.
    std_msgs::ColorRGBA SignedDistanceToRosColor(float sdf) const;
  };
}

#endif
