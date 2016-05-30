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
#include <atom_map/ApproximateAtomMap.h>
#include <atom_map/AtomMapParameters.h>
#include <atom_map/CsvWriter.h>

#include <visualization_msgs/Marker.h>
#include <algorithm>
#include <Eigen/Dense>
#include <iostream>
#include <random>
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
Atom::Ptr AtomMap::GetNearestAtom(float x, float y, float z) {
  std::vector<Atom::Ptr> neighbors;
  if (!map_.GetKNearestNeighbors(x, y, z, 1, &neighbors))
    return nullptr;

  if (neighbors.size() != 1) return nullptr;
  return neighbors[0];
}

void AtomMap::GetSignedDistance(float x, float y, float z, float* distance,
                                float* variance) {
  CHECK_NOTNULL(distance);
  CHECK_NOTNULL(variance);

  // Find nearest neighbors.
  std::vector<Atom::Ptr> neighbors;
  if (!map_.GetKNearestNeighbors(x, y, z, num_neighbors_, &neighbors)) {
    ROS_WARN("%s: Error in extracting nearest neighbors.", name_.c_str());
    *distance = std::numeric_limits<float>::infinity();
    *variance = std::numeric_limits<float>::infinity();
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

      float cov = CovarianceKernel(p1, p2);
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
float AtomMap::GetProbability(float x, float y, float z) {
  std::vector<Atom::Ptr> neighbors;
  if (!map_.RadiusSearch(x, y, z, radius_, &neighbors)) {
    ROS_WARN("%s: Error in radius search during GetProbability().",
             name_.c_str());
    return -1.0;
  }

  if (neighbors.size() == 0) {
    ROS_WARN("%s: Nearest neighbor is too far away.", name_.c_str());
    return -1.0;
  }

  // There can only be one neighbor.
  return neighbors[0]->GetProbability();
}

void AtomMap::MaybeInsertAtom(const Atom::Ptr& atom) {
  // If the AtomKdtree is empty, just insert this atom.
  if (map_.Size() == 0) {
    if (!map_.Insert(atom))
      ROS_WARN("%s: Error inserting a new Atom.", name_.c_str());
    return;
  }

  std::vector<Atom::Ptr> neighbors;
  pcl::PointXYZ position;
  position.x = atom->GetPosition().X();
  position.y = atom->GetPosition().Y();
  position.z = atom->GetPosition().Z();

  if (!map_.RadiusSearch(position, 2.0 * radius_ - 1e-4, &neighbors)) {
    ROS_WARN("%s: Error in radius search during Update().", name_.c_str());
    return;
  }

  if (neighbors.size() == 0) {
    if (!map_.Insert(atom))
      ROS_WARN("%s: Error inserting a new Atom.", name_.c_str());
    return;
  }

  // Handle case where sample lies inside an existing Atom.
  else {
    for (size_t jj = 0; jj < neighbors.size(); jj++) {
      Atom::Ptr neighbor = neighbors[jj];

      if (atom->GetDistanceTo(neighbor) > 2.0 * radius_) {
        ROS_WARN("%s: Neighbor is too far away: %lf > %lf.", name_.c_str(),
                 atom->GetDistanceTo(neighbor), 2.0 * radius_);
        continue;
      }

      // Compute overlap fraction.
      const float weight = atom->ComputeOverlapFraction(neighbor);
      if (weight >= 0.0 && weight <= 1.0) {
        // Update occupancy.
        if (update_occupancy_)
          neighbor->UpdateProbability(atom->GetProbability(), weight);

        // Update signed distance.
        if (update_signed_distance_) {
          neighbor->UpdateSignedDistance(atom->GetSignedDistance(), weight);
        }
      } else
        ROS_WARN("%s: Weight was out of bounds(%lf). Distance between atoms was %lf.",
                 name_.c_str(), weight, neighbor->GetDistanceTo(atom));
    }
  }
}

// Return a list of all Atoms in the map.
const std::vector<Atom::Ptr>& AtomMap::GetAtoms() const {
  return map_.GetAtoms();
}

// Get the neighbors of an Atom in the implicit graph. Returns false
// if the Atom is not itself in the map.
bool AtomMap::GetConnectedNeighbors(Atom::Ptr& atom,
                                    std::vector<Atom::Ptr>* connected) {
  CHECK_NOTNULL(connected);
  connected->clear();

  // Unpack.
  const gu::Vec3f p = atom->GetPosition();

  // Do a radius search.
  std::vector<Atom::Ptr> neighbors;
  if (!map_.RadiusSearch(p(0), p(1), p(2), connectedness_radius_, &neighbors))
    return false;

  // Filter out the given Atom if it exists. Otherwise return false.
  // Obviously only add if the atom is considered free.
  bool contains_query = false;
  for (size_t ii = 0; ii < neighbors.size(); ii++) {
    Atom::Ptr neighbor = neighbors[ii];

    if (!contains_query && neighbor->GetDistanceTo(atom) < 1e-4)
      contains_query = true;
    else if (neighbor->GetProbability() < free_threshold_)
      connected->push_back(neighbor);
  }

  if (!contains_query) return false;
  return true;
}

// Update the map given a set of observations.
void AtomMap::Update(const PointCloud::ConstPtr& cloud,
                     const pcl::PointXYZ& robot) {
  // Create an ApproximateAtomMap from this cloud.
  AtomMapParameters params;
  params.radius_ = radius_;
  params.min_scan_range_ = min_scan_range_;
  params.max_scan_range_ = max_scan_range_;
  params.update_occupancy_ = update_occupancy_;
  params.update_signed_distance_ = update_signed_distance_;
  params.surface_normal_radius_ = surface_normal_radius_;
  params.max_occupied_backoff_ = max_occupied_backoff_;
  params.max_normal_backoff_ = max_normal_backoff_;
  params.max_samples_normal_ = max_samples_normal_;
  params.angular_resolution_ = angular_resolution_;
  params.angular_interleaving_ = angular_interleaving_;
  params.lambda_ = lambda_;
  params.probability_hit_ = probability_hit_;
  params.probability_miss_ = probability_miss_;
  params.name_ = name_ + "/single_scan";

  ApproximateAtomMap small_map(params, cloud, robot);

  // Merge miniature map with the larger atom map.
  for (const auto& atom : small_map.GetAtoms())
    MaybeInsertAtom(atom);
}

// Load parameters and register callbacks.
bool AtomMap::RegisterCallbacks(const ros::NodeHandle& n) {
  ros::NodeHandle node(n);

  // Publishers.
  full_occupancy_publisher_ = node.advertise<visualization_msgs::Marker>(
      full_occupancy_topic_.c_str(), 0);
  full_sdf_publisher_ =
      node.advertise<visualization_msgs::Marker>(full_sdf_topic_.c_str(), 0);

  return true;
}

bool AtomMap::LoadParameters(const ros::NodeHandle& n) {
  // General parameters.
  if (!pu::Get("atom/radius", radius_)) return false;
  if (!pu::Get("atom/connectedness_radius", connectedness_radius_)) return false;
  if (!pu::Get("atom/update_occupancy", update_occupancy_)) return false;
  if (!pu::Get("atom/update_signed_distance", update_signed_distance_))
    return false;
  if (!pu::Get("atom/only_show_occupied", only_show_occupied_)) return false;
  if (!pu::Get("atom/occupied_threshold", occupied_threshold_)) return false;
  if (!pu::Get("atom/free_threshold", free_threshold_)) return false;
  if (!pu::Get("atom/sdf_threshold", sdf_threshold_)) return false;
  if (!pu::Get("atom/probability_hit", probability_hit_)) return false;
  if (!pu::Get("atom/probability_miss", probability_miss_)) return false;
  if (!pu::Get("atom/probability_clamp_high", probability_clamp_high_))
    return false;
  if (!pu::Get("atom/probability_clamp_low", probability_clamp_low_))
    return false;
  if (!pu::Get("atom/surface_normal_radius", surface_normal_radius_))
    return false;
  if (!pu::Get("atom/full_occupancy_topic", full_occupancy_topic_))
    return false;
  if (!pu::Get("atom/full_sdf_topic", full_sdf_topic_)) return false;
  if (!pu::Get("atom/fixed_frame_id", fixed_frame_id_)) return false;

  // Raytracing parameters.
  if (!pu::Get("atom/raytracing/min_scan_range", min_scan_range_)) return false;
  if (!pu::Get("atom/raytracing/max_scan_range", max_scan_range_)) return false;
  if (!pu::Get("atom/raytracing/max_occupied_backoff", max_occupied_backoff_))
    return false;
  if (!pu::Get("atom/raytracing/max_normal_backoff", max_normal_backoff_)) return false;
  if (!pu::Get("atom/raytracing/max_samples_normal", max_samples_normal_)) return false;
  if (!pu::Get("atom/raytracing/angular_resolution", angular_resolution_)) return false;
  if (!pu::Get("atom/raytracing/angular_interleaving", angular_interleaving_))
    return false;
  if (!pu::Get("atom/raytracing/lambda", lambda_)) return false;

  // Gaussian process regression parameters.
  if (!pu::Get("atom/gp/gamma", gamma_)) return false;
  if (!pu::Get("atom/gp/noise", noise_variance_)) return false;
  if (!pu::Get("atom/gp/num_neighbors", num_neighbors_)) return false;

  // Radius and clamping are constant across all atoms. Call static setters.
  Atom::SetRadius(radius_);
  Atom::SetProbabilityClamps(probability_clamp_low_,
                             probability_clamp_high_);

  return true;
}

// Apply the covariance kernel function. This is just the simplest option, but
// there are many other valid kernels to consider. For example one easy change
// would be to scale covariance by the 'variance' estimate of signed distance
// at each Atom.
float AtomMap::CovarianceKernel(const pcl::PointXYZ& p1,
                                 const pcl::PointXYZ& p2) {
  const float dx = p1.x - p2.x;
  const float dy = p1.y - p2.y;
  const float dz = p1.z - p2.z;
  return std::exp(-gamma_ * (dx * dx + dy * dy + dz * dz));
}

// Publish the full AtomMap colored by occupancy probability. Optionally,
// only show the occupied atoms.
void AtomMap::PublishFullOccupancy() const {
  if (full_occupancy_publisher_.getNumSubscribers() <= 0) return;

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
    const float probability_occupied = atoms[ii]->GetProbability();

    // Maybe only show if probably occupied.
    if (!only_show_occupied_ || probability_occupied > occupied_threshold_) {
      const gu::Vec3f p = atoms[ii]->GetPosition();
      m.points.push_back(gr::ToRosPoint(p));
      m.colors.push_back(ProbabilityToRosColor(probability_occupied));
    }
  }

  full_occupancy_publisher_.publish(m);
}

// Publish the full AtomMap colored by signed distance.
void AtomMap::PublishFullSignedDistance() const {
  if (full_sdf_publisher_.getNumSubscribers() <= 0) return;

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
    const float sdf = atoms[ii]->GetSignedDistance();

    // Maybe only show if probably occupied.
    if (!only_show_occupied_ || std::abs(sdf) < sdf_threshold_) {
      const gu::Vec3f p = atoms[ii]->GetPosition();
      m.points.push_back(gr::ToRosPoint(p));
      m.colors.push_back(SignedDistanceToRosColor(sdf));
    }
  }

  full_sdf_publisher_.publish(m);
}

// Convert a probability of occupancy to a ROS color. Red is more likely to be
// occupied, blue is more likely to be free.
std_msgs::ColorRGBA AtomMap::ProbabilityToRosColor(float probability) const {
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
  c.a = 1.0;

  return c;
}

// Convert a signed distnace to a ROS color. Red is probably close to a surface,
// and blue is probably far from a surface.
std_msgs::ColorRGBA AtomMap::SignedDistanceToRosColor(float sdf) const {
  const float min_dist = (only_show_occupied_) ? -sdf_threshold_ : 0.0;
  const float max_dist = (only_show_occupied_) ? sdf_threshold_
    : std::max(map_.GetMaxDistance(), std::abs(map_.GetMinDistance()));

  std_msgs::ColorRGBA c;
  c.r = (sdf - min_dist) / (max_dist - min_dist);
  c.g = 0.0;
  c.b = 1.0 - (sdf - min_dist) / (max_dist - min_dist);
  c.a = 1.0;

  return c;
}

// Save to '.csv' file. First line contains just the number of Atoms in the
// map. Subsequent lines contain x, y, z coordinates followed by sdf.
void AtomMap::Save(const std::string& filename) const {
  ROS_INFO("%s: Saving to %s.", name_.c_str(), filename.c_str());
  const std::vector<Atom::Ptr> atoms = map_.GetAtoms();

  // Open a file.
  file::CsvWriter writer(filename);

  // Write the number of Atoms.
  std::vector<int> first_line;
  first_line.push_back(static_cast<int>(atoms.size()));
  writer.WriteLine(first_line);

  // Write each atom.
  for (size_t ii = 0; ii < atoms.size(); ii++) {
    Atom::Ptr atom = atoms[ii];
    const gu::Vec3f position = atom->GetPosition();
    const double sdf = atom->GetSignedDistance();

    // Pack.
    std::vector<double> data;
    data.push_back(position(0));
    data.push_back(position(1));
    data.push_back(position(2));
    data.push_back(sdf);

    // Write.
    writer.WriteLine(data);
  }

  return;
}


}  // namespace atom
