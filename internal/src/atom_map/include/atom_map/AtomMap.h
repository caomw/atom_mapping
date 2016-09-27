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
#include <atom_map/CsvWriter.h>

#include <ros/ros.h>
#include <std_msgs/ColorRGBA.h>
#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>
#include <pcl_ros/point_cloud.h>
#include <glog/logging.h>
#include <Eigen/Core>
#include <math.h>

using Eigen::Vector3f;

typedef pcl::PointCloud<pcl::PointXYZ> PointCloud;

namespace atom {
  class AtomMap {
  public:
    AtomMap();
    ~AtomMap();

    // Initialize.
    bool Initialize(const ros::NodeHandle& n);

    // Getters.
#if 0
    void GetSignedDistance(float x, float y, float z,
                           float* distance, float* variance);
    float GetProbability(float x, float y, float z);
#endif

    Atom::Ptr GetNearestAtom(float x, float y, float z);
    const std::vector<Atom::Ptr>& GetAtoms() const;

    // Get the neighbors of an Atom in the implicit graph. Returns false
    // if the Atom is not itself in the map.
    bool GetConnectedNeighbors(Atom::Ptr& atom,
                               std::vector<Atom::Ptr>* connected);

    // Update.
    void Update(const PointCloud::ConstPtr& cloud, const pcl::PointXYZ& robot);

    // Publishing.
    void PublishOccupancy() const;
    void PublishSignedDistance() const;
    void PublishPointCloud() const;

    // Save to '.csv' file. First line contains just the number of Atoms in the
    // map. Subsequent lines contain x, y, z coordinates followed by sdf.
    void Save(const std::string& filename = "atoms.csv") const;

  private:
    // A kdtree to hold all the Atoms.
    AtomKdtree map_;

    // Atomic radius.
    float radius_;

    // Connectedness radius. Atoms within this distance of any given Atom are considered
    // to be 'connected' in the implicit graph.
    float connectedness_radius_;

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

    // Flag for applying voxel grid filter to sampled Atoms for each new scan.
    bool voxel_grid_;

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
    std::string occupancy_topic_;
    std::string sdf_topic_;
    std::string pcld_topic_;
    bool only_show_occupied_;
    float occupied_threshold_; // Above this is considered occupied.
    float free_threshold_;     // Below this is considered free space.
    float sdf_threshold_;      // Smaller than this is considered near a surface.
    ros::Publisher occupancy_pub_;
    ros::Publisher sdf_pub_;
    ros::Publisher pcld_pub_;

    // Publish on a timer.
    ros::Timer timer;

    // Initialization and naming.
    bool initialized_;
    std::string name_;

    // Load parameters and register callbacks.
    bool LoadParameters(const ros::NodeHandle& n);
    bool RegisterCallbacks(const ros::NodeHandle& n);

    // Timer callback. Call all publishers on a timer.
    void TimerCallback(const ros::TimerEvent& event) const;

    // Try to add a single atom.
    void MaybeInsertAtom(const Atom::Ptr& atom);

    // Apply the covariance kernel function.
    float CovarianceKernel(const pcl::PointXYZ& p1, const pcl::PointXYZ& p2);

    // Convert a probability of occupancy to a ROS color.
    std_msgs::ColorRGBA ProbabilityToRosColor(float probability) const;

    // Convert a signed distnace value to a ROS color.
    std_msgs::ColorRGBA SignedDistanceToRosColor(float sdf) const;
  };
}

#endif
