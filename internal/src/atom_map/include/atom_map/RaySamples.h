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

#ifndef ATOM_MAPPING_RAY_SAMPLES_H
#define ATOM_MAPPING_RAY_SAMPLES_H

#include <pcl/point_types.h>
#include <vector>

namespace atom {
  struct RaySamples {
    // List of points and signed distances along rays between the observer
    // and the surface.
    std::vector<pcl::PointXYZ> ray_points_;
    std::vector<double> ray_distances_;

    // List of points and signed distances normal to the surface, on the side of
    // and the observer.
    std::vector<pcl::PointXYZ> normal_points_;
    std::vector<double> normal_distances_;

    // List of points and signed distances normal to the surface and on the
    // 'occupied' side of the surface.
    std::vector<pcl::PointXYZ> occupied_points_;
    std::vector<double> occupied_distances_;

    // Clear all of these lists.
    void ClearAll() {
      ray_points_.clear();
      ray_distances_.clear();
      normal_points_.clear();
      normal_distances_.clear();
      occupied_points_.clear();
      occupied_distances_.clear();
    }
  }; // struct RaySamples
} // namespace atom

#endif
