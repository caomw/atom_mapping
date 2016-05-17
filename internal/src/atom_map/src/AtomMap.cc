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

#include <visualization_msgs/Marker.h>
#include <Eigen/Dense>
#include <iostream>
#include <math.h>

namespace gu = geometry_utils;
namespace gr = gu::ros;
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
  void AtomMap::GetSignedDistance(double x, double y, double z,
                                  double* distance, double* variance) {
    CHECK_NOTNULL(distance);
    CHECK_NOTNULL(variance);

    // Find nearest neighbors.
    std::vector<Atom::Ptr> neighbors;
    if (!map_.GetKNearestNeighbors(x, y, z, num_neighbors_, &neighbors)) {
      ROS_WARN("%s: Error in extracting nearest neighbors.", name_.c_str());
      *distance = std::numeric_limits<double>::infinity();
      *variance = std::numeric_limits<double>::infinity();
    }

    // Generate covariance matrix for the training data.
    const size_t kNumNeighbors = neighbors.size();
    Eigen::MatrixXd K11(kNumNeighbors, kNumNeighbors);
    for (size_t ii = 0; ii < kNumNeighbors; ii++) {
      pcl::PointXYZ p1;
      p1.x = neighbors[ii]->GetPosition()(0);
      p1.y = neighbors[ii]->GetPosition()(1);
      p1.z = neighbors[ii]->GetPosition()(2);

      for (size_t jj = 0; jj < ii; jj++) {
        pcl::PointXYZ p2;
        p2.x = neighbors[jj]->GetPosition()(0);
        p2.y = neighbors[jj]->GetPosition()(1);
        p2.z = neighbors[jj]->GetPosition()(2);

        double cov = CovarianceKernel(p1, p2);
        K11(ii, jj) = cov;
        K11(jj, ii) = cov;
      }

      K11(ii, ii) = 1.0 + noise_variance_;
    }

    // Generate query cross covariance and training distance column vectors.
    Eigen::VectorXd K12(kNumNeighbors);
    Eigen::VectorXd training_dists(kNumNeighbors);
    pcl::PointXYZ q(x, y, z);
    for (size_t ii = 0; ii < kNumNeighbors; ii++) {
      pcl::PointXYZ p;
      p.x = neighbors[ii]->GetPosition()(0);
      p.y = neighbors[ii]->GetPosition()(1);
      p.z = neighbors[ii]->GetPosition()(2);

      K12(ii) = CovarianceKernel(p, q);
      training_dists(ii) = neighbors[ii]->GetSignedDistance();
    }

    // Gaussian conditioning. Use LDLT (stable Cholesky) for speed and because
    // we cannot guarantee that the covariance matrix of the training data
    // K11 is strictly positive definite (it may be semidefinite or nearly
    // semidefinite, in which case regular Cholesky is unstable).
    Eigen::LDLT<Eigen::MatrixXd> cholesky(K11);
    *distance = K12.transpose() * cholesky.solve(training_dists);
    *variance = 1.0 - K12.transpose() * cholesky.solve(K12);
  }

  // Find probability of occupancy. Return -1 if error or if this point does
  // not lie within an Atom.
  double AtomMap::GetProbability(double x, double y, double z) {
    std::vector<Atom::Ptr> neighbors;
    if (!map_.RadiusSearch(x, y, z, radius_, &neighbors)) {
      ROS_WARN("%s: Error in radius search during GetProbability().", name_.c_str());
      return -1.0;
    }

    if (neighbors.size() == 0) {
      ROS_WARN("%s: Nearest neighbor is too far away.", name_.c_str());
      return -1.0;
    }

    // There can only be one neighbor.
    return neighbors[0]->GetProbability();
  }

  // Update the map for this observation.
  void AtomMap::Update(const pcl::PointXYZ& point, const pcl::Normal& normal,
                       const pcl::PointXYZ& robot) {
    // Sample the ray.
    SampledRay samples;
    SampleRay(point, normal, robot, &samples);

    // Try to insert occupied Atoms.
    for (size_t ii = 0; ii < samples.occupied_points_.size(); ii++) {
      MaybeInsertAtom(samples.occupied_points_[ii],
                      samples.occupied_distances_[ii]);
    }

    // Try to insert Atoms along the ray.
    if (update_occupancy_) {
      for (size_t ii = 0; ii < samples.ray_points_.size(); ii++) {
        MaybeInsertAtom(samples.ray_points_[ii],
                        samples.ray_distances_[ii]);
      }
    }

    // Try to insert Atoms along the normal.
    if (update_signed_distance_) {
      for (size_t ii = 0; ii < samples.normal_points_.size(); ii++) {
        MaybeInsertAtom(samples.normal_points_[ii],
                        samples.normal_distances_[ii]);
      }
    }
  }

  // Insert Atom at this position into the tree. Handle options regarding
  // updating occupancy and signed distance.
  void AtomMap::MaybeInsertAtom(const pcl::PointXYZ& position, double sdf) {
    Atom::Ptr atom = Atom::Create(radius_);
    atom->SetPosition(gu::Vec3(position.x, position.y, position.z));

    // If the AtomKdtree is empty, just insert this atom.
    if (map_.Size() == 0) {
      // Set probability of occupancy.
      if (update_occupancy_) {
        if (sdf > 0.0)
          atom->SetProbability(probability_miss_);
        else
          atom->SetProbability(probability_hit_);
      }

      // Set signed distance.
      if (update_signed_distance_)
        atom->SetSignedDistance(sdf);

      // Insert.
      if (!map_.Insert(atom))
        ROS_WARN("%s: Error inserting a new Atom.", name_.c_str());

      return;
    }

    std::vector<Atom::Ptr> neighbors;
    if (!map_.RadiusSearch(position, 2.0 * radius_ - 1e-4, &neighbors)) {
      ROS_WARN("%s: Error in radius search during Update().", name_.c_str());
      return;
    }

    // Handle case where sample lies more than twice the atomic radius
    // from its nearest neighbor.
    if (neighbors.size() == 0) {
      // Set probability of occupancy.
      if (update_occupancy_) {
        if (sdf > 0.0)
          atom->SetProbability(probability_miss_);
        else
          atom->SetProbability(probability_hit_);
      }

      // Set signed distance.
      if (update_signed_distance_)
        atom->SetSignedDistance(sdf);

      // Insert into kdtree. Insertion here automatically updates neighbors
      // in the implicit graph structure of the kdtree.
      if (!map_.Insert(atom))
        ROS_WARN("%s: Error inserting a new Atom.", name_.c_str());
    }

    // Handle case where sample lies inside an existing Atom.
    else {
      for (size_t jj = 0; jj < neighbors.size(); jj++) {
        Atom::Ptr neighbor = neighbors[jj];

        // Compute overlap fraction.
        const double weight = atom->ComputeOverlapFraction(neighbor);

        // Update occupancy.
        if (update_occupancy_) {
          if (sdf > 0.0)
            neighbor->UpdateProbability(probability_miss_, weight);
          else
            neighbor->UpdateProbability(probability_hit_, weight);
        }

        // Update signed distance.
        if (update_signed_distance_) {
          if (weight < 0.0 || weight > 1.0)
            ROS_WARN("Weight was out of bounds (%5.4f). Distance between atoms was %5.4f.",
                     weight, neighbor->GetDistanceTo(atom));
          else
            neighbor->UpdateSignedDistance(sdf, weight);
        }
      }
    }
  }

  // Update the map given a set of observations.
  void AtomMap::Update(const PointCloud::ConstPtr& cloud,
                       const pcl::PointXYZ& robot) {
    // Get surface normals.
    pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
    pcl::search::KdTree<pcl::PointXYZ>::Ptr
      tree(new pcl::search::KdTree<pcl::PointXYZ>());
    ne.setInputCloud(cloud);
    ne.setSearchMethod(tree);
    ne.setRadiusSearch(surface_normal_radius_);
    ne.setViewPoint(robot.x, robot.y, robot.z);

    pcl::PointCloud<pcl::Normal>::Ptr
      normals(new pcl::PointCloud<pcl::Normal>);
    ne.compute(*normals);

    // Ensure that we have the same number of normals as we do input points.
    if (normals->points.size() != cloud->points.size()) {
      ROS_WARN("%s: Error calculating surface normals. Incorrect number of points.",
               name_.c_str());
      return;
    }

    // Update each point individually.
    for (size_t ii = 0; ii < cloud->points.size(); ii++) {
      pcl::Normal normal = normals->points[ii];
      Update(cloud->points[ii], normal, robot);
    }
  }

  // Load parameters and register callbacks.
  bool AtomMap::RegisterCallbacks(const ros::NodeHandle& n) {
    ros::NodeHandle node(n);

    // Publishers.
    full_occupancy_publisher_ =
      node.advertise<visualization_msgs::Marker>(full_occupancy_topic_.c_str(), 0);
    full_sdf_publisher_ =
      node.advertise<visualization_msgs::Marker>(full_sdf_topic_.c_str(), 0);

    return true;
  }

  bool AtomMap::LoadParameters(const ros::NodeHandle& n) {
    if (!pu::Get("atom/gamma", gamma_)) return false;
    if (!pu::Get("atom/radius", radius_)) return false;
    if (!pu::Get("atom/noise", noise_variance_)) return false;
    if (!pu::Get("atom/probability_hit", probability_hit_)) return false;
    if (!pu::Get("atom/probability_miss", probability_miss_)) return false;
    if (!pu::Get("atom/full_occupancy_topic", full_occupancy_topic_)) return false;
    if (!pu::Get("atom/full_sdf_topic", full_sdf_topic_)) return false;
    if (!pu::Get("atom/fixed_frame_id", fixed_frame_id_)) return false;
    if (!pu::Get("atom/only_show_occupied", only_show_occupied_)) return false;
    if (!pu::Get("atom/occupied_threshold", occupied_threshold_)) return false;
    if (!pu::Get("atom/sdf_threshold", sdf_threshold_)) return false;
    if (!pu::Get("atom/surface_normal_radius", surface_normal_radius_)) return false;
    if (!pu::Get("atom/max_occupied_backoff", max_occupied_backoff_)) return false;
    if (!pu::Get("atom/max_normal_backoff", max_normal_backoff_)) return false;
    if (!pu::Get("atom/num_neighbors", num_neighbors_)) return false;
    if (!pu::Get("atom/update_occupancy", update_occupancy_)) return false;
    if (!pu::Get("atom/update_signed_distance", update_signed_distance_)) return false;

    return true;
  }

  // Sample a ray. Given a robot position and a measured point, discretize the
  // ray from sensor to observation and insert/update atoms along the way.
  // In particular, discretize such that the atoms closest to the surface are
  // tangent to it. Behind the surface, walk along the surface normal (which by
  // default points toward the robot's side of the surface). Also optionally
  // walk along the surface normal but on the unoccupied side of the surface.
  void AtomMap::SampleRay(const pcl::PointXYZ& point, const pcl::Normal& normal,
                          const pcl::PointXYZ& robot, SampledRay* samples) {
    CHECK_NOTNULL(samples);
    samples->ClearAll();

    // Compute the range to the observed point and the unit direction.
    double dx = robot.x - point.x;
    double dy = robot.y - point.y;
    double dz = robot.z - point.z;
    const double range = std::sqrt(dx*dx + dy*dy + dz*dz);

    // Handle non-returns, i.e. range == 0.
    if (range < 1e-6) {
      return;
    }

    dx /= range; dy /= range; dz /= range;

    // Unpack normal vector. Not const because it will be set to dx/dy/dz if NAN.
    double nx = normal.normal_x;
    double ny = normal.normal_y;
    double nz = normal.normal_z;

    if (isnan(nx) || isnan(ny) || isnan(nz)) {
      nx = dx; ny = dy; nz = dz;
    }

    // Start at the surface and walk away from the robot, ideally along the normal
    // vector, if it is not NAN. Otherwise just trace along the ray.
    const size_t num_samples_back =
      static_cast<size_t>(0.5 * max_occupied_backoff_ / radius_);
    const double step_size_back =
      0.5 * max_occupied_backoff_ / static_cast<double>(num_samples_back);
    for (size_t ii = 0; ii < num_samples_back; ii++) {
      const double backoff = static_cast<double>(2 * ii + 1) * step_size_back;

      pcl::PointXYZ p;
      p.x = point.x - backoff * nx;
      p.y = point.y - backoff * ny;
      p.z = point.z - backoff * nz;

      samples->occupied_points_.push_back(p);
      samples->occupied_distances_.push_back(-backoff);
    }

    if (update_signed_distance_) {
      // Start at the surface and walk along the surface normal (toward free space).
      const size_t num_samples_normal =
        static_cast<size_t>(0.5 * max_normal_backoff_ / radius_);
      const double step_size_normal =
        0.5 * max_normal_backoff_ / static_cast<double>(num_samples_normal);
      for (size_t ii = 0; ii < num_samples_normal; ii++) {
        const double backoff = static_cast<double>(2 * ii + 1) * step_size_normal;

        pcl::PointXYZ p;
        p.x = point.x + backoff * nx;
        p.y = point.y + backoff * ny;
        p.z = point.z + backoff * nz;

        samples->normal_points_.push_back(p);
        samples->normal_distances_.push_back(backoff);
      }
    }

    if (update_occupancy_) {
      // Start at the surface and walk toward the robot. Initially, backoff just
      // enough to be tangent to the surface.
      const double front_initial_backoff = radius_ / (dx * nx + dy * ny + dz * nz);

      // If this backoff distance is greater than the range to the robot, just
      // ignore this point.
      if (front_initial_backoff < range) {
        const size_t num_samples_front =
          static_cast<size_t>(0.5 * (range - front_initial_backoff) / radius_);
        const double step_size_front =
          0.5 * (range - front_initial_backoff) / static_cast<double>(num_samples_front);
        for (size_t ii = 0; ii < num_samples_front; ii++) {
          const double backoff =
            static_cast<double>(2 * ii + 1) * step_size_front + front_initial_backoff;

          pcl::PointXYZ p;
          p.x = point.x + backoff * dx;
          p.y = point.y + backoff * dy;
          p.z = point.z + backoff * dz;

          samples->ray_points_.push_back(p);
          samples->ray_distances_.push_back(backoff);
        }
      }
    }
  }

  // Apply the covariance kernel function. This is just the simplest option, but
  // there are many other valid kernels to consider. For example one easy change
  // would be to scale covariance by the 'variance' estimate of signed distance
  // at each Atom.
  double AtomMap::CovarianceKernel(const pcl::PointXYZ& p1,
                                   const pcl::PointXYZ& p2) {
    const double dx = p1.x - p2.x;
    const double dy = p1.y - p2.y;
    const double dz = p1.z - p2.z;
    return std::exp(-gamma_ * (dx*dx + dy*dy + dz*dz));
  }

  // Publish the full AtomMap colored by occupancy probability. Optionally,
  // only show the occupied atoms.
  void AtomMap::PublishFullOccupancy() const {
    if (full_occupancy_publisher_.getNumSubscribers() <= 0)
      return;

    // Initialize marker.
    visualization_msgs::Marker m;
    m.header.frame_id = fixed_frame_id_;
    // m.header.stamp = ros::Time();
    m.ns = fixed_frame_id_;
    m.id = 0;
    m.action = visualization_msgs::Marker::ADD;
    m.type = visualization_msgs::Marker::SPHERE_LIST;
    m.color.r = 0.0;
    m.color.g = 0.4;
    m.color.b = 0.8;
    m.color.a = 0.2;
    m.scale.x = 2.0 * radius_;
    m.scale.y = 2.0 * radius_;
    m.scale.z = 2.0 * radius_;
    m.pose = gr::ToRosPose(gu::Transform3::Identity());

    // Loop over all atoms and add to marker.
    const std::vector<Atom::Ptr> atoms = map_.GetAtoms();
    ROS_INFO("%s: Publishing %lu atoms.", name_.c_str(), atoms.size());

    for (size_t ii = 0; ii < atoms.size(); ii++) {
      const double probability_occupied = atoms[ii]->GetProbability();

      // Maybe only show if probably occupied.
      if (!only_show_occupied_ || probability_occupied > occupied_threshold_) {
        gu::Vec3 p = atoms[ii]->GetPosition();
        m.points.push_back(gr::ToRosPoint(p));
        m.colors.push_back(ProbabilityToRosColor(probability_occupied));
      }
    }

    full_occupancy_publisher_.publish(m);
  }

  // Publish the full AtomMap colored by signed distance.
  void AtomMap::PublishFullSignedDistance() const {
    if (full_sdf_publisher_.getNumSubscribers() <= 0)
      return;

    // Initialize marker.
    visualization_msgs::Marker m;
    m.header.frame_id = fixed_frame_id_;
    // m.header.stamp = ros::Time();
    m.ns = fixed_frame_id_;
    m.id = 0;
    m.action = visualization_msgs::Marker::ADD;
    m.type = visualization_msgs::Marker::SPHERE_LIST;
    m.color.r = 0.0;
    m.color.g = 0.4;
    m.color.b = 0.8;
    m.color.a = 0.2;
    m.scale.x = 2.0 * radius_;
    m.scale.y = 2.0 * radius_;
    m.scale.z = 2.0 * radius_;
    m.pose = gr::ToRosPose(gu::Transform3::Identity());

    // Loop over all atoms and add to marker.
    const std::vector<Atom::Ptr> atoms = map_.GetAtoms();
    ROS_INFO("%s: Publishing %lu atoms.", name_.c_str(), atoms.size());

    for (size_t ii = 0; ii < atoms.size(); ii++) {
      const double sdf = atoms[ii]->GetSignedDistance();

      // Maybe only show if probably occupied.
      if (!only_show_occupied_ || std::abs(sdf) < sdf_threshold_) {
        gu::Vec3 p = atoms[ii]->GetPosition();
        m.points.push_back(gr::ToRosPoint(p));
        m.colors.push_back(SignedDistanceToRosColor(sdf));
      }
    }

    full_sdf_publisher_.publish(m);
  }

  // Convert a probability of occupancy to a ROS color. Red is more likely to be
  // occupied, blue is more likely to be free.
  std_msgs::ColorRGBA AtomMap::ProbabilityToRosColor(double probability) const {
#ifdef ENABLE_DEBUG_MESSAGES
    if (probability < 0.0) {
      ROS_ERROR("%s: Probability is out of bounds.", name_.c_str());
      probability = 0.0;
    }
    if (probability > 1.0) {
      ROS_ERROR("%s: Probability is out of bounds.", name_.c_str());
      probability = 1.0;
    }
#endif
    std_msgs::ColorRGBA c;
    c.r = probability;
    c.g = 0.0;
    c.b = 1.0 - probability;
    c.a = 0.2;

    return c;
  }

  // Convert a signed distnace to a ROS color. Red is probably close to a surface,
  // and blue is probably far from a surface.
  std_msgs::ColorRGBA AtomMap::SignedDistanceToRosColor(double sdf) const {
    // const double min_dist = map_.GetMinDistance();
    // const double max_dist = map_.GetMaxDistance();
    const double min_dist = -sdf_threshold_;
    const double max_dist = sdf_threshold_;

    std_msgs::ColorRGBA c;
    c.r = 1.0 - (sdf - min_dist) / (max_dist - min_dist);
    c.g = 0.0;
    c.b = (sdf - min_dist) / (max_dist - min_dist);
    c.a = 1.0;

    return c;
  }

} // namespace atom
