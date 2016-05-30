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

#ifndef ATOM_MAPPING_ATOM_PATH_H
#define ATOM_MAPPING_ATOM_PATH_H

#include <atom_map/Atom.h>
#include <geometry_utils/GeometryUtilsROS.h>

#include <ros/ros.h>
#include <visualization_msgs/Marker.h>
#include <vector>

namespace gu = geometry_utils;
namespace gr = gu::ros;

namespace atom {
  struct AtomPath {
    // A list to store the Atoms.
    std::vector<Atom::Ptr> atoms_;

    // Compute length.
    double Length() const {
      double length = 0.0;

      for (size_t ii = 1; ii < atoms_.size(); ii++)
        length += atoms_[ii - 1]->GetDistanceTo(atoms_[ii]);

      return length;
    }

    // A function to publish a visualization of this path.
    void Visualize(ros::Publisher& pub, const std::string& fixed_frame_id) const {
      if (pub.getNumSubscribers() <= 0) return;

      // Don't do anything if there are no atoms in the path.
      if (atoms_.size() == 0) return;

      // Make empty marker messages.
      std_msgs::ColorRGBA line_color;
      line_color.r = 0.0;
      line_color.g = 0.4;
      line_color.b = 0.8;
      line_color.a = 1.0;

      visualization_msgs::Marker line_marker;
      line_marker.header.frame_id = fixed_frame_id;
      line_marker.ns = fixed_frame_id;
      line_marker.id = 1;
      line_marker.action = visualization_msgs::Marker::ADD;
      line_marker.type = visualization_msgs::Marker::LINE_STRIP;
      line_marker.color = line_color;
      line_marker.scale.x = 0.25 * atoms_[0]->GetRadius();
      line_marker.pose = gr::ToRosPose(gu::Transform3::Identity());

      std_msgs::ColorRGBA atom_color;
      atom_color.r = 0.0;
      atom_color.g = 0.8;
      atom_color.b = 0.4;
      line_color.a = 1.0;

      visualization_msgs::Marker atom_marker;
      atom_marker.header.frame_id = fixed_frame_id;
      atom_marker.ns = fixed_frame_id;
      atom_marker.id = 2;
      atom_marker.action = visualization_msgs::Marker::ADD;
      atom_marker.type = visualization_msgs::Marker::SPHERE_LIST;
      atom_marker.color = atom_color;
      atom_marker.scale.x = 2.0 * atoms_[0]->GetRadius();
      atom_marker.scale.y = 2.0 * atoms_[0]->GetRadius();
      atom_marker.scale.z = 2.0 * atoms_[0]->GetRadius();
      atom_marker.pose = gr::ToRosPose(gu::Transform3::Identity());

      // If only one Atom, only populate the Atom marker.
      if (atoms_.size() == 1) {
        const gu::Vec3f p = atoms_[0]->GetPosition();
        atom_marker.points.push_back(gr::ToRosPoint(p));
        atom_marker.colors.push_back(atom_color);

        pub.publish(atom_marker);
      } else {
        // Loop over all atoms and add to marker.
        for (size_t ii = 0; ii < atoms_.size(); ii++) {
          const gu::Vec3f p = atoms_[ii]->GetPosition();
          line_marker.points.push_back(gr::ToRosPoint(p));
          line_marker.colors.push_back(line_color);

          atom_marker.points.push_back(gr::ToRosPoint(p));
          atom_marker.colors.push_back(atom_color);
        }

        // Publish.
        pub.publish(line_marker);
        pub.publish(atom_marker);
      }
    }
  }; // struct AtomPath
} // namespace atom

#endif
