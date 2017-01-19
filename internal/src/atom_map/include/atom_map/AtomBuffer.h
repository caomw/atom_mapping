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

#ifndef ATOM_MAPPING_ATOM_BUFFER_H
#define ATOM_MAPPING_ATOM_BUFFER_H

#include <atom_map/Atom.h>
#include <atom_map/AtomMapParameters.h>
#include <atom_map/RaySamples.h>

#include <ros/ros.h>
#include <std_msgs/ColorRGBA.h>
#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>
#include <pcl_ros/point_cloud.h>
#include <glog/logging.h>
#include <random>
#include <math.h>

typedef pcl::PointCloud<pcl::PointXYZ> PointCloud;

namespace atom {
  class AtomBuffer {
  public:
    ~AtomBuffer();
    AtomBuffer(const AtomMapParameters& params,
               const PointCloud::ConstPtr& cloud,
               const pcl::PointXYZ& robot);

    // Return a list of all Atoms in the map.
    const std::vector<Atom::Ptr>& GetAtoms() const;

  private:
    // An simple buffer to hold all the Atoms. Since we are only
    // inserting one scan (ideally with interleaving), and we are going
    // to then use an exact kdtree to merge this with the big map, we
    // only need to keep a list of Atoms to insert later.
    std::vector<Atom::Ptr> atoms_;

    // Random number generation.
    std::random_device rd_;

    // Atomic radius.
    const float radius_;

    // Occupancy mode flag.
    const bool occupancy_mode_;

    // Min and max scan ranges.
    const float min_scan_range_;
    const float max_scan_range_;

    // Radius for nearest neighbor search for surface normal extraction.
    const float surface_normal_radius_;

    // Maximum distance to trace a ray inside of an object. This should be small,
    // in order to ensure we don't accidentally trace all the way through. However,
    // it is definitely necessary in order to ensure proper surface detection.
    const float max_occupied_backoff_;

    // Maximum distance to trace the surface normal inside (supposedly) free space.
    const float max_normal_backoff_;

    // Maximum number of Atoms to lay down along the normal vector into free space.
    const int max_samples_normal_;

    // Sensor angular resolution. This is the smallest anglue between two range
    // measurements. We also provide a flag to turn on angular interleaving.
    const float angular_resolution_;
    const bool angular_interleaving_;
    const float lambda_;

    // Flag for applying voxel grid filter to sampled Atoms.
    const bool voxel_grid_;

    // Probability of hits and misses for occupancy updates.
    const float probability_hit_;
    const float probability_miss_;

    // Name, for ROS warnings and debug messages.
    const std::string name_;

    // Add a single atom.
    void InsertAtom(const pcl::PointXYZ& position, float sdf,
                    std::vector<Atom::Ptr>* raw);

    // Sample a ray and do a probabilistic and signed distance update.
    // Given a robot position and a measured point, discretize the ray from sensor
    // to observation and return vectors of points and distances. Append to the
    // output RaySamples argument.
    void SampleRay(const pcl::PointXYZ& point, const pcl::Normal& normal,
                   const pcl::PointXYZ& robot, RaySamples* samples);
  }; // class AtomBuffer
} // namespace atom

#endif
