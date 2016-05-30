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

#ifndef ATOM_MAPPING_ASTAR_PLANNER_H
#define ATOM_MAPPING_ASTAR_PLANNER_H

#include <atom_map/Planner.h>
#include <atom_map/Atom.h>
#include <atom_map/AtomMap.h>
#include <atom_map/AtomPath.h>
#include <geometry_utils/Vector3.h>

#include <glog/logging.h>
#include <math.h>
#include <queue>
#include <vector>

namespace gu = geometry_utils;

namespace atom {
  // A functor class for use in the priority queue.
  struct AtomPathComparitor {
    bool operator()(AtomPath& path1, AtomPath& path2) {
      return path1.ExpectedTotalLength() <= path2.ExpectedTotalLength();
    }
  };


  class AStarPlanner : public Planner {
  public:
    AStarPlanner(AtomMap* map, size_t max_iters = 500)
      : Planner(map), max_iters_(max_iters) {}
    ~AStarPlanner() {}

    // Plan a path.
    bool Plan(const gu::Vec3f& start_position, const gu::Vec3f& goal_position,
              AtomPath* path) const;

  private:
    // Maximum number of iterations during A* search.
    const size_t max_iters_;
  };

  // ------------------------------- IMPLEMENTATION --------------------------- //
  bool AStarPlanner::Plan(const gu::Vec3f& start_position, const gu::Vec3f& goal_position,
                     AtomPath* path) const {
    CHECK_NOTNULL(path);

    // Find nearest Atoms to start and goal positions.
    Atom::Ptr start =
      map_->GetNearestAtom(start_position(0), start_position(1), start_position(2));
    Atom::Ptr goal =
      map_->GetNearestAtom(goal_position(0), goal_position(1), goal_position(2));

    if (start == nullptr || goal == nullptr) return false;

    // Create a comparitor.
    std::priority_queue<AtomPath, std::vector<AtomPath>, AtomPathComparitor> pq;

    // Steps:
    // (1) Start from the 'start' Atom. Push it into the queue.
    // (2) Pop from the queue. If goal, end.
    // (3) Otherwise, add all neighbors to the queue.
    // (4) Goto (2).

    // (1) Put the 'start' Atom on the queue.
    AtomPath initial_path;
    initial_path.goal_ = goal;
    initial_path.atoms_.push_back(start);
    pq.push(initial_path);

    // (2) Loop till top of the queue is the goal.
    size_t iters = 0;
    while (iters++ < max_iters_ && pq.size() > 0 &&
           pq.top().atoms_.back()->GetDistanceTo(goal) > 1e-4) {
      ROS_INFO("Iteration : %lu", iters);
      AtomPath shortest = pq.top();
      const size_t num_atoms = shortest.atoms_.size();
      if (num_atoms == 0) return false;

      // Extract most recent Atom.
      ROS_INFO("Getting most recent Atom.");
      Atom::Ptr current_atom = shortest.atoms_[num_atoms - 1];

      // Maybe extract previous Atom.
      ROS_INFO("Getting previous Atom.");
      Atom::Ptr previous_atom =
        (num_atoms > 1) ? shortest.atoms_[num_atoms - 2] : nullptr;

      // (3) Add all neighboring paths that do not include the previous Atom.
      std::vector<Atom::Ptr> neighbors;
      if (!map_->GetConnectedNeighbors(current_atom, &neighbors))
        return false;
      ROS_INFO("Adding %lu neighboring paths.", neighbors.size());

      for (size_t ii = 0; ii < neighbors.size(); ii++) {
        Atom::Ptr neighbor = neighbors[ii];

        if (previous_atom == nullptr ||
            neighbor->GetDistanceTo(previous_atom) > 1e-4) {
          AtomPath next_path;
          next_path.goal_ = goal;
          next_path.atoms_ = shortest.atoms_;
          next_path.atoms_.push_back(neighbor);
          pq.push(next_path);
        }
      }

      // Pop.
      pq.pop();
    }

    if (iters >= max_iters_) {
      ROS_WARN("A* planner ran out of iterations.");
      return false;
    }

    if (pq.size() == 0) {
      ROS_WARN("A* planner attempted to pop an empty queue.");
      return false;
    }

    // Return the path at the top of the queue.
    path->atoms_ = pq.top().atoms_;
    return true;
  }

}
#endif
