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
  void AtomMap::Update(const pcl::PointXYZ& point,
                       const pcl::PointXYZ& robot) {
    std::vector<pcl::PointXYZ> samples;
    std::vector<double> signed_distances;

    // Sample the ray.
    SampleRay(point, robot, &samples, &signed_distances);

    // For each sample, update occupancy and signed distance function
    // in existing Atoms, or add new Atoms to the map.
    for (size_t ii = 0; ii < samples.size(); ii++) {
      const double sdf = signed_distances[ii];

      // If the AtomKdtree is empty, just insert this atom.
      if (map_.Size() == 0) {
        Atom::Ptr atom = Atom::Create(radius_);
        atom->SetPosition(gu::Vec3(samples[ii].x, samples[ii].y, samples[ii].z));

        // Set probability of occupancy.
        if (sdf > 0.0)
          atom->SetProbability(probability_miss_);
        else
          atom->SetProbability(probability_hit_);

        // Set signed distance.
        atom->SetSignedDistance(sdf);

        // Insert.
        if (!map_.Insert(atom))
          ROS_WARN("%s: Error inserting a new Atom.", name_.c_str());

        continue;
      }

      std::vector<Atom::Ptr> neighbors;
      if (!map_.RadiusSearch(samples[ii], 2.0 * radius_ - 1e-6, &neighbors)) {
        ROS_WARN("%s: Error in radius search during Update().", name_.c_str());
        continue;
      }

      // Handle case where sample lies more than twice the atomic radius
      // from its nearest neighbor.
      else if (neighbors.size() == 0) {
        Atom::Ptr atom = Atom::Create(radius_);
        atom->SetPosition(gu::Vec3(samples[ii].x, samples[ii].y, samples[ii].z));

        // Set probability of occupancy.
        if (sdf > 0.0)
          atom->SetProbability(probability_miss_);
        else
          atom->SetProbability(probability_hit_);

        // Set signed distance.
        atom->SetSignedDistance(sdf);

        // Insert into kdtree. Insertion here automatically updates neighbors
        // in the implicit graph structure of the kdtree.
        if (!map_.Insert(atom))
          ROS_WARN("%s: Error inserting a new Atom.", name_.c_str());
      }

      // Handle case where sample lies inside an existing Atom.
      // TODO: Handle partial overlaps.
      else {
        for (size_t jj = 0; jj < neighbors.size(); jj++) {
          if (neighbors[jj]->Contains(samples[ii])) {
            // Update probability of occupancy.
            if (sdf > 0.0)
              neighbors[jj]->UpdateProbability(probability_miss_);
            else
              neighbors[jj]->UpdateProbability(probability_hit_);

            // Update signed distance.
            neighbors[jj]->UpdateSignedDistance(sdf);
            break;
          }
        }
      }
    }
  }

  // Update the map given a set of observations.
  void AtomMap::Update(const PointCloud::ConstPtr& cloud,
                       const pcl::PointXYZ& robot) {
    for (size_t ii = 0; ii < cloud->points.size(); ii++) {
      Update(cloud->points[ii], robot);
    }
  }

  // Load parameters and register callbacks.
  bool AtomMap::RegisterCallbacks(const ros::NodeHandle& n) {
    ros::NodeHandle node(n);
    atom_publisher_ =
      node.advertise<visualization_msgs::Marker>(visualization_topic_.c_str(), 0);
    return true;
  }

  bool AtomMap::LoadParameters(const ros::NodeHandle& n) {
    if (!pu::Get("atom/gamma", gamma_)) return false;
    if (!pu::Get("atom/radius", radius_)) return false;
    if (!pu::Get("atom/thickness", max_surface_thickness_)) return false;
    if (!pu::Get("atom/noise", noise_variance_)) return false;
    if (!pu::Get("atom/probability_hit", probability_hit_)) return false;
    if (!pu::Get("atom/probability_miss", probability_miss_)) return false;
    if (!pu::Get("atom/visualization_topic", visualization_topic_)) return false;
    if (!pu::Get("atom/fixed_frame_id", fixed_frame_id_)) return false;

    return true;
  }

  // Sample a ray. Given a robot position and a measured point, discretize the
  // ray from sensor to observation and insert/update atoms along the way.
  // In particular, discretize such that the atoms closest to the surface are
  // tangent to it.
  void AtomMap::SampleRay(const pcl::PointXYZ& point,
                          const pcl::PointXYZ& robot,
                          std::vector<pcl::PointXYZ>* samples,
                          std::vector<double>* signed_distances) {
    CHECK_NOTNULL(samples);
    CHECK_NOTNULL(signed_distances);

    samples->clear();
    signed_distances->clear();

    // Compute the range to the observed point and the unit direction.
    double dx = robot.x - point.x;
    double dy = robot.y - point.y;
    double dz = robot.z - point.z;
    double range = std::sqrt(dx*dx + dy*dy + dz*dz);

    // Handle non-returns, i.e. range == 0.
    if (range < 1e-6) {
      return;
    }

    dx /= range; dy /= range; dz /= range;

    // Start at the surface and walk toward the robot.
    const size_t num_samples_front = static_cast<size_t>(range / (2.0 * radius_));
    for (size_t ii = 0; ii < num_samples_front; ii++) {
      const double backoff = static_cast<double>(2 * ii + 1) * radius_;

      pcl::PointXYZ p;
      p.x = point.x + backoff * dx;
      p.y = point.y + backoff * dy;
      p.z = point.z + backoff * dz;

      samples->push_back(p);
      signed_distances->push_back(backoff);
    }

    // Start at the surface and walk away from the robot.
    const size_t num_samples_back =
      static_cast<size_t>(max_surface_thickness_ / (2.0 * radius_));
    for (size_t ii = 0; ii < num_samples_back; ii++) {
      const double backoff = static_cast<double>(2 * ii + 1) * radius_;

      pcl::PointXYZ p;
      p.x = point.x - backoff * dx;
      p.y = point.y - backoff * dy;
      p.z = point.z - backoff * dz;

      samples->push_back(p);
      signed_distances->push_back(-backoff);
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

  // Publish all atoms.
  void AtomMap::Publish() const {
    if (atom_publisher_.getNumSubscribers() <= 0)
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
      gu::Vec3 p = atoms[ii]->GetPosition();

      m.points.push_back(gr::ToRosPoint(p));
      m.colors.push_back(ProbabilityToRosColor(atoms[ii]->GetProbability()));
    }

    atom_publisher_.publish(m);
  }

  // Convert a probability of occupancy to a ROS color.
  std_msgs::ColorRGBA AtomMap::ProbabilityToRosColor(double probability) const {
    if (probability < 0.0) {
      ROS_ERROR("%s: Probability is out of bounds.", name_.c_str());
      probability = 0.0;
    }

    if (probability > 1.0) {
      ROS_ERROR("%s: Probability is out of bounds.", name_.c_str());
      probability = 1.0;
    }

    std_msgs::ColorRGBA c;
    c.r = probability;
    c.g = 0.0;
    c.b = 1.0 - probability;
    c.a = 0.2;

    return c;
  }

} // namespace atom
