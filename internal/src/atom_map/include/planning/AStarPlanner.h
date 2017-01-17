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

#include <atom_map/Atom.h>
#include <atom_map/AtomMap.h>
#include <planning/Planner.h>
#include <planning/AtomPath.h>
#include <planning/ShortestPathsTree.h>

#include <Eigen/Core>
#include <glog/logging.h>
#include <math.h>
#include <queue>
#include <vector>

namespace sp = atom::shortest_paths;

using Eigen::Vector3f;

namespace atom {
  class AStarPlanner : public Planner {
  public:
    AStarPlanner(AtomMap* map,
                 float sdf_manifold = -1.0,
                 size_t max_iters = std::numeric_limits<size_t>::max())
      : Planner(map), max_iters_(max_iters), sdf_manifold_(sdf_manifold) {}

    ~AStarPlanner() {}

    // Plan a path.
    bool Plan(const Vector3f& start_position, const Vector3f& goal_position,
              AtomPath* path) const;

  private:
    // Maximum number of iterations during A* search.
    const size_t max_iters_;

    // Signed distance value on which to plan. A* will try to plan a path
    // that stays close to this manifold.
    const float sdf_manifold_;
  }; // class AStarPlanner

  // ------------------------------- IMPLEMENTATION --------------------------- //
  bool AStarPlanner::Plan(const Vector3f& start_position,
                          const Vector3f& goal_position,
                          AtomPath* path) const {
    CHECK_NOTNULL(path);

    // Find nearest Atoms to start and goal positions.
    Atom::Ptr start = map_->GetAtomContaining(start_position(0),
                                              start_position(1), start_position(2));
    Atom::Ptr goal = map_->GetAtomContaining(start_position(0),
                                             start_position(1), start_position(2));

    if (start == nullptr || goal == nullptr) return false;

    // Create a comparitor.
    std::priority_queue<sp::Node::Ptr, std::vector<sp::Node::Ptr>,
                        sp::Tree::NodeComparitor> pq;

    // Create an empty ShortestPathsTree.
    sp::Tree tree;
    sp::Tree::SetGoal(goal);
    sp::Tree::SetSDFManifold(sdf_manifold_);

    // Steps:
    // (1) Start from the 'start' Atom. Push it into the queue.
    // (2) Pop from the queue. If goal, end.
    // (3) Otherwise, add all neighbors to the queue.
    // (4) Goto (2).

    // (1) Put the 'start' Atom on the queue.
    sp::Node::Ptr start_node = sp::Node::Create(start);
    tree.SetRoot(start_node);
    pq.push(start_node);

    // (2) Loop till top of the queue is the goal.
    size_t iters = 0;
    while (iters++ < max_iters_ && pq.size() > 0 &&
           pq.top()->GetAtom()->GetDistanceTo(goal) > 1e-4) {
      sp::Node::Ptr shortest = pq.top();
      Atom::Ptr current_atom = shortest->GetAtom();
#if ENABLE_DEBUG_MESSAGES
      if (iters % 100 == 0)
        ROS_INFO("A* iteration %lu: length = %f",
                 iters, shortest->GetPathLength());
#endif

      // (3) Add all neighboring Atoms that are not already in the tree.
      std::vector<Atom::Ptr> neighbors;
      if (!map_->GetConnectedNeighbors(current_atom, &neighbors))
        return false;

      for (size_t ii = 0; ii < neighbors.size(); ii++) {
        Atom::Ptr neighbor = neighbors[ii];

        if (!tree.Contains(neighbor)) {
          sp::Node::Ptr next_node = sp::Node::Create(neighbor);
          tree.Attach(shortest, next_node);
          pq.push(next_node);
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
    pq.top()->GetPath(path);
    return true;
  }

}
#endif
