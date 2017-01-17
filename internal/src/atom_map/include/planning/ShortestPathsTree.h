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

#ifndef ATOM_MAPPING_SHORTEST_PATHS_TREE_H
#define ATOM_MAPPING_SHORTEST_PATHS_TREE_H

#include <atom_map/Atom.h>
#include <planning/AtomPath.h>

#include <glog/logging.h>
#include <math.h>
#include <unordered_map>
#include <memory>
#include <vector>

namespace atom {
namespace shortest_paths {
  class Node {
  public:
    typedef std::shared_ptr<Node> Ptr;
    ~Node() {}

    // Factory method.
    static Ptr Create(Atom::Ptr& atom) {
      Ptr node(new Node(atom));
      return node;
    }

    // Setter.
    void SetParent(Node::Ptr& parent) {
      parent_ = parent;
      path_length_ =
        parent->path_length_ + atom_->GetDistanceTo(parent->atom_);
    }

    // Getters.
    float GetPathLength() const { return path_length_; }
    Node::Ptr GetParent() { return parent_; }
    Atom::Ptr GetAtom() { return atom_; }

    // Get the path back to the root from a given node.
    void GetPath(AtomPath* path) {
      CHECK_NOTNULL(path);
      path->atoms_.clear();

      Node::Ptr parent = parent_;

      // Handle case where we're already at the root.
      if (parent == nullptr) {
        path->atoms_.push_back(atom_);
        return;
      }

      // Set up a vector to hold the reversed path as we trace up the tree.
      std::vector<Atom::Ptr> reversed_path;
      reversed_path.push_back(atom_);

      while (parent != nullptr) {
        reversed_path.push_back(parent->atom_);
        parent = parent->parent_;
      }

      // Reverse the list.
      for (size_t ii = 1; ii < reversed_path.size() - 1; ii ++) {
        path->atoms_.push_back(reversed_path[reversed_path.size() - ii]);
      }

      return;
    }

  private:
    Atom::Ptr atom_;
    Node::Ptr parent_;
    float path_length_;

    Node(Atom::Ptr& atom) : atom_(atom), path_length_(0.0) {}
  }; // class Node

  class Tree {
  public:
    Tree() {}
    ~Tree() {}

    // Static setters for the goal atom and sdf manifold.
    // Compiler insists on static because NodeComparitor needs to see it. Why?
    static void SetGoal(Atom::Ptr& goal) { goal_ = goal; }
    static void SetSDFManifold(float sdf_manifold) { sdf_manifold_ = sdf_manifold; }

    // Set root.
    void SetRoot(Node::Ptr& root) { registry_.insert({root->GetAtom(), root}); }

    // Attach two nodes.
    void Attach(Node::Ptr& parent, Node::Ptr& child) {
      child->SetParent(parent);
      registry_.insert({child->GetAtom(), child});
    }

    // Check if a node is contained in the tree.
    bool Contains(Atom::Ptr& node) { return registry_.count(node) > 0; }

    // Comparitor. If sdf_manifold_ is specified, then add the perpendicular distance
    // from the node to the manifold into the score.
    struct NodeComparitor {
      bool operator()(Node::Ptr& node1, Node::Ptr& node2) {
        const float manifold_dist1 =
          (sdf_manifold_ < 0.0) ? 0.0
          : std::abs(node1->GetAtom()->GetSignedDistance() - sdf_manifold_);
        const float manifold_dist2 =
          (sdf_manifold_ < 0.0) ? 0.0
          : std::abs(node2->GetAtom()->GetSignedDistance() - sdf_manifold_);

        const float score1 = manifold_dist1 +
          node1->GetPathLength() + node1->GetAtom()->GetDistanceTo(goal_);
        const float score2 = manifold_dist2 +
          node2->GetPathLength() + node2->GetAtom()->GetDistanceTo(goal_);

        return score1 > score2;
      }
    }; // class NodeComparitor

  protected:
    static Atom::Ptr goal_;
    static float sdf_manifold_;

  private:
    std::unordered_map<Atom::Ptr, Node::Ptr> registry_;

  }; // class Tree

  // ------------------------------ IMPLEMENTATION ---------------------------------- //
  Atom::Ptr Tree::goal_ = nullptr;
  float Tree::sdf_manifold_ = -1.0;

} // namespace shortest_paths
} // namespace atom

#endif
