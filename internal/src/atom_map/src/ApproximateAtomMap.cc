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

#include <atom_map/ApproximateAtomMap.h>
#include <atom_map/VoxelGrid.h>

#include <Eigen/Dense>
#include <iostream>
#include <math.h>

namespace gu = geometry_utils;
namespace gr = gu::ros;

namespace atom {

  ApproximateAtomMap::~ApproximateAtomMap() {}
  ApproximateAtomMap::ApproximateAtomMap(const AtomMapParameters& params,
                                         const PointCloud::ConstPtr& cloud,
                                         const pcl::PointXYZ& robot)
    : radius_(params.radius_),
      min_scan_range_(params.min_scan_range_),
      max_scan_range_(params.max_scan_range_),
      update_occupancy_(params.update_occupancy_),
      update_signed_distance_(params.update_signed_distance_),
      surface_normal_radius_(params.surface_normal_radius_),
      max_occupied_backoff_(params.max_occupied_backoff_),
      max_normal_backoff_(params.max_normal_backoff_),
      max_samples_normal_(params.max_samples_normal_),
      angular_resolution_(params.angular_resolution_),
      angular_interleaving_(params.angular_interleaving_),
      lambda_(params.lambda_),
      voxel_grid_(params.voxel_grid_),
      probability_hit_(params.probability_hit_),
      probability_miss_(params.probability_miss_),
      name_(params.name_) {
    // Steps:
    // (1) Get surface normals for incoming points.
    // (2) Sample all rays along the point cloud.
    // (3) Shuffle points randomly.
    // (4) Insert points into the map.
    // (5) Voxel grid filter.

    // (1) Get surface normals for incoming points.
    pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
    pcl::search::KdTree<pcl::PointXYZ>::Ptr
      tree(new pcl::search::KdTree<pcl::PointXYZ>());
    ne.setInputCloud(cloud);
    ne.setSearchMethod(tree);
    ne.setRadiusSearch(surface_normal_radius_);
    ne.setViewPoint(robot.x, robot.y, robot.z);

    pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
    ne.compute(*normals);

    // Ensure that we have the same number of normals as we do input points.
    if (normals->points.size() != cloud->points.size()) {
      ROS_WARN("%s: Error calculating surface normals. Incorrect number of points.",
               name_.c_str());
      return;
    }

    // (2) Sample all rays along the point cloud.
    RaySamples samples;
    for (size_t ii = 0; ii < cloud->points.size(); ii++) {
      pcl::Normal normal = normals->points[ii];
      SampleRay(cloud->points[ii], normal, robot, &samples);
    }

    // (3) Shuffle points randomly if no angular interleaving.
    std::vector<size_t> occupied_indices(samples.occupied_points_.size());
    std::iota(occupied_indices.begin(), occupied_indices.end(), 0);

    std::vector<size_t> ray_indices(samples.ray_points_.size());
    std::iota(ray_indices.begin(), ray_indices.end(), 0);

    std::vector<size_t> normal_indices(samples.normal_points_.size());
    std::iota(normal_indices.begin(), normal_indices.end(), 0);

    if (!angular_interleaving_) {
      std::shuffle(occupied_indices.begin(), occupied_indices.end(),
                   std::mt19937{std::random_device{}()});
      std::shuffle(ray_indices.begin(), ray_indices.end(),
                   std::mt19937{std::random_device{}()});
      std::shuffle(normal_indices.begin(), normal_indices.end(),
                   std::mt19937{std::random_device{}()});
    }

    // (4) Insert points into this map.
    std::vector<Atom::Ptr> raw;
    std::vector<Atom::Ptr>* ptr = (voxel_grid_) ? &raw : &atoms_;
    for (const auto& idx : occupied_indices)
      InsertAtom(samples.occupied_points_[idx],
                 samples.occupied_distances_[idx], ptr);

    if (update_occupancy_) {
      for (const auto& idx : ray_indices)
        InsertAtom(samples.ray_points_[idx],
                   samples.ray_distances_[idx], ptr);
    }
    if (update_signed_distance_) {
      for (const auto& idx : normal_indices)
        InsertAtom(samples.normal_points_[idx],
                   samples.normal_distances_[idx], ptr);
    }

    // (5) Voxel grid filter.
    if (voxel_grid_) {
      VoxelGrid grid_filter(2.0 * radius_);
      grid_filter.Filter(raw, &atoms_);
    }
  }

  // Insert Atom into the provided list.
  void ApproximateAtomMap::InsertAtom(const pcl::PointXYZ& position, float sdf,
                                      std::vector<Atom::Ptr>* raw) {
    CHECK_NOTNULL(raw);

    Atom::Ptr atom = Atom::Create();
    atom->SetPosition(gu::Vec3f(position.x, position.y, position.z));

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
    raw->push_back(atom);

    return;
  }

  // Return a list of all Atoms in the map.
  const std::vector<Atom::Ptr>& ApproximateAtomMap::GetAtoms() const {
    return atoms_;
  }

  // Sample a ray. Given a robot position and a measured point, discretize the
  // ray from sensor to observation and insert/update atoms along the way.
  // In particular, discretize such that the atoms closest to the surface are
  // tangent to it. Behind the surface, walk along the surface normal (which by
  // default points toward the robot's side of the surface). Also optionally
  // walk along the surface normal but on the unoccupied side of the surface.
  // Moreover, along the ray to the robot, optionally interleave Atoms as they
  // approach the sensor.
  void ApproximateAtomMap::SampleRay(const pcl::PointXYZ& point, const pcl::Normal& normal,
                          const pcl::PointXYZ& robot, RaySamples* samples) {
    CHECK_NOTNULL(samples);

    // Unpack point for convenience.
    float px = point.x;
    float py = point.y;
    float pz = point.z;

    // Compute the range to the observed point and the unit direction.
    float dx = robot.x - px;
    float dy = robot.y - py;
    float dz = robot.z - pz;
    float range = std::sqrt(dx * dx + dy * dy + dz * dz);

    // Handle non-returns, i.e. range == 0.
    if (range < 1e-6) return;

    // Stop if range is below lower bound.
    if (range < min_scan_range_) return;

    dx /= range;
    dy /= range;
    dz /= range;

    // Unpack normal vector. Not const because it will be set to dx/dy/dz if NAN
    // or if out of range.
    float nx = normal.normal_x;
    float ny = normal.normal_y;
    float nz = normal.normal_z;

    if (isnan(nx) || isnan(ny) || isnan(nz) || range > max_scan_range_) {
      nx = dx;
      ny = dy;
      nz = dz;
    }

    // Only update occupied and normal directions if range is below upper bound.
    if (range <= max_scan_range_) {
      // Start at the surface and walk away from the robot, ideally along the normal
      // vector, if it is not NAN. Otherwise just trace along the ray.
      const size_t num_samples_back =
        static_cast<size_t>(0.5 * max_occupied_backoff_ / radius_);
      const float step_size_back =
        0.5 * max_occupied_backoff_ / static_cast<float>(num_samples_back);
      for (size_t ii = 0; ii < num_samples_back; ii++) {
        const float backoff = static_cast<float>(2 * ii + 1) * step_size_back;

        pcl::PointXYZ p;
        p.x = px - backoff * nx;
        p.y = py - backoff * ny;
        p.z = pz - backoff * nz;

        samples->occupied_points_.push_back(p);
        samples->occupied_distances_.push_back(-backoff);
      }

      if (update_signed_distance_) {
        // Start at the surface and walk along the surface normal (toward free
        // space).
        const size_t dense_samples_normal =
          static_cast<size_t>(0.5 * max_normal_backoff_ / radius_);
        const size_t num_samples_normal =
          (dense_samples_normal < max_samples_normal_)
          ? dense_samples_normal
          : static_cast<size_t>(max_samples_normal_);
        const float step_size_normal =
          0.5 * max_normal_backoff_ / static_cast<float>(num_samples_normal);
        for (size_t ii = 0; ii < num_samples_normal; ii++) {
          const float backoff = static_cast<float>(2 * ii + 1) * step_size_normal;

          pcl::PointXYZ p;
          p.x = px + backoff * nx;
          p.y = py + backoff * ny;
          p.z = pz + backoff * nz;

          samples->normal_points_.push_back(p);
          samples->normal_distances_.push_back(backoff);
        }
      }
    }

    if (update_occupancy_) {
      // Start at the surface and walk toward the robot. Initially, backoff just
      // enough to be tangent to the surface.
      float ray_initial_backoff = radius_ / (dx * nx + dy * ny + dz * nz);

      // If range is outside the max, then add the difference to the initial backoff.
      if (range > max_scan_range_)
        ray_initial_backoff += range - max_scan_range_;

      // If this backoff distance is greater than the range to the robot, just
      // ignore this point.
      if (ray_initial_backoff >= range) return;

      // Compute number of samples.
      const size_t num_samples_ray =
        static_cast<size_t>(0.5 * (range - ray_initial_backoff) / radius_);
      const float step_size_ray = 0.5 * (range - ray_initial_backoff) /
        static_cast<float>(num_samples_ray);

      // This is the counter to track which 'delta' we are at: i.e. how many rays
      // can possibly intersect here.
      float delta_number = 1.0;
      //float delta = radius_ / tan(0.5 * angular_resolution_);
      float delta = radius_ * sqrt(2.0 / (1.0 - cos(angular_resolution_)));
      float probability = (1.0 + lambda_) / (lambda_ + 2.0 * delta_number + 1.0);

      // Random number generation. This is really important to get right.
      std::default_random_engine rng(rd_());
      std::uniform_real_distribution<float> unif(0.0, 1.0);
      for (size_t ii = 0; ii < num_samples_ray; ii++) {
        const float backoff =
          static_cast<float>(2 * ii + 1) * step_size_ray + ray_initial_backoff;

        // Updating delta parameters.
        if (angular_interleaving_ && range - backoff < delta) {
          //delta_number =
          //  std::ceil(2.0 * atan2(radius_, range - backoff) / angular_resolution_);
          //delta = radius_ / tan(0.5 * delta_number * angular_resolution_);
          delta_number =
            std::ceil(acos(1 - 2.0 * radius_ * radius_ /
                           ((range - backoff) * (range - backoff))) / angular_resolution_);
          delta = radius_ * sqrt(2.0 / (1.0 - cos(delta_number * angular_resolution_)));
          probability = (1.0 + lambda_) / (lambda_ + 2.0 * delta_number - 1.0);
        }

        // Interleaving.
        if (!angular_interleaving_ || unif(rng) < probability) {
          pcl::PointXYZ p;
          p.x = px + backoff * dx;
          p.y = py + backoff * dy;
          p.z = pz + backoff * dz;

          samples->ray_points_.push_back(p);
          samples->ray_distances_.push_back(backoff);
        }
      }
    }
  }
}  // namespace atom
