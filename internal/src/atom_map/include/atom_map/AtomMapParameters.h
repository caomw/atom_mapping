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

#ifndef ATOM_MAPPING_ATOM_MAP_PARAMS_H
#define ATOM_MAPPING_ATOM_MAP_PARAMS_H

namespace atom {
  struct AtomMapParameters {
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

    // Flag for applying voxel grid filter to sampled Atoms.
    bool voxel_grid_;

    // Name, for ROS warnings and debug messages.
    std::string name_;
  }; // struct AtomMapParameters
} // namespace atom

#endif
