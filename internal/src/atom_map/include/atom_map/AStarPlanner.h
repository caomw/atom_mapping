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

#ifndef ATOM_MAPPING_PLANNER_H
#define ATOM_MAPPING_PLANNER_H

#include <atom_map/Planner.h>
#include <atom_map/Atom.h>
#include <atom_map/AtomMap.h>
#include <atom_map/AtomPath.h>

#include <glog/logging.h>
#include <math.h>
#include <queue>
#include <vector>

namespace atom {
  // A functor class for use in the priority queue. Comparitor should return
  // true if path1 is shorter than path2 (i.e. it should prefer shorter paths).
  class AtomPathComparitor {
    bool operator()(AtomPath& path1, AtomPath& path2) {
      return path1.Length() < path2.Length();
    }
  }


  class AStarPlanner : public Planner {
  public:
    AStarPlanner(AtomMap* map, size_t max_iters)
      : Planner(map), max_iters_(max_iters) {}
    ~AStarPlanner() {}

    // Plan a path.
    bool Plan(Atom::Ptr& start, Atom::Ptr& goal, AtomPath* path) const {
      CHECK_NOTNULL(path);
      std::priority_queue<AtomPath, std::vector<AtomPath>, AtomPathComparitor> pq;

      // Steps:
      // (1) Start from the 'start' Atom. Push it into the queue.
      // (2) Pop from the queue. If goal, end.
      // (3) Otherwise, add all neighbors to the queue.
      // (4) Goto (2).

      // (1) Put the 'start' Atom on the queue.
      AtomPath initial_path;
      initial_path.atoms_.push_back(start);
      pq.push(initial_path);

      // (2) Loop till top of the queue is the goal.
      size_t iters = 0;
      while (iters++ < max_iters_ &&
             pq.top().atoms_.back()->GetDistanceTo(goal) > 1e-4) {
        const AtomPath shortest = pq.top();
        const size_t num_atoms = shortest.atoms_.size();
        if (num_atoms == 0) return false;
        pq.pop();

        // Extract most recent Atom.
        const Atom::Ptr current_atom = shortest.atoms_[num_atoms - 1];

        // Maybe extract previous Atom.
        const Atom::Ptr previous_atom =
          (num_atoms > 1) ? shortest.atoms_[num_atoms - 2] : nullptr;

        // (3) Add all neighboring paths that do not include the previous Atom.
        const std::vector<Atom::Ptr> neighbors =
          map->GetConnectedNeighbors(current_atom);

        for (size_t ii = 0; ii < neighbors.size(); ii++) {
          const Atom::Ptr neighbor = neighbors[ii];

          if (previous_atom != nullptr && neighbor->GetDistanceTo(previous_atom) > 1e-4) {
            AtomPath next_path;
            next_path.atoms_ = shortest.atoms_;
            next_path.atoms_.push_back(neighbor);
            pq.push(next_path);
          }
        }
      }

      // Return the path at the top of the queue.
      path->atoms_ = pq.top().atoms_;
      return true;
    };

  private:
    // Maximum number of iterations during A* search.
    const size_t max_iters_;
  };
}

#endif
