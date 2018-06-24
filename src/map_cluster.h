//
// Created by Jason Shen on 5/21/18.
//

#ifndef HIERARCHICAL_MAP_MAP_CLUSTER_H
#define HIERARCHICAL_MAP_MAP_CLUSTER_H

#include <nlohmann/json.hpp>
#include <unordered_map>
#include <unordered_set>
#include "graph.h"
#include "psdd_manager.h"
#include "types.h"
using json = nlohmann::json;
using std::string;
using std::unordered_map;
using std::unordered_set;
class LeafConstraintHandler;
class MapCluster {
 public:
  MapCluster(ClusterSize cluster_index, std::string cluster_name,
             std::unordered_set<NodeSize> nodes, MapCluster *left_child,
             MapCluster *right_child);
  void SetInternalExternalEdgesFromEdgeList(const std::vector<Edge *> &edges);
  static MapCluster *MapClusterFromNetworkJson(
      const json &json_spec, ClusterSize cluster_index,
      const std::string &cluster_name,
      const std::unordered_map<std::string, MapCluster *>
          &constructed_clusters);
  const std::unordered_set<NodeSize> &nodes() const;
  MapCluster *left_child() { return left_child_; }
  MapCluster *right_child() { return right_child_; }
  const std::unordered_map<Edge *, NodeSize> &external_edges() const {
    return external_edges_;
  }
  const std::unordered_set<Edge *> &internal_edges() const {
    return internal_edges_;
  }
  const std::string &cluster_name() { return cluster_name_; }
  void InitConstraint(
      const std::unordered_map<MapCluster *, SddLiteral> *cluster_variable_map,
      const std::unordered_map<Edge *, SddLiteral> *edge_variable_map,
      const std::unordered_map<Edge *, MapCluster *> *edge_cluster_map,
      PsddManager *manager, Vtree *local_vtree,
      LeafConstraintHandler *leaf_constraint_handler);
  Vtree *GenerateLocalVtree(
      const std::unordered_map<MapCluster *, SddLiteral> *cluster_variable_map,
      const std::unordered_map<Edge *, SddLiteral> *edge_variable_map) const;
  static Vtree *GenerateVtreeUsingMinFillOfPrimalGraph(
      const unordered_set<Edge *> &edges,
      const unordered_map<Edge *, SddLiteral> *edge_variable_map);
  PsddNode *empty_path_constraint() const;
  PsddNode *internal_path_constraint() const;
  PsddNode *terminal_path_constraint(Edge *entering_edge) const;
  PsddNode *non_terminal_path_constraint(Edge *entering_edge_1,
                                         Edge *entering_edge_2) const;

 protected:
  PsddNode *ConstructEmptyPathConstraintForInternalCluster() const;
  PsddNode *ConstructInternalPathConstraintForInternalCluster() const;
  PsddNode *ConstructTerminalPathConstraintForInternalCluster(Edge *edge) const;
  PsddNode *ConstructNonTerminalPathConstraintForInternalCluster(
      Edge *edge_1, Edge *edge_2) const;
  ClusterSize cluster_index_;
  string cluster_name_;
  unordered_set<NodeSize> nodes_;
  MapCluster *left_child_;
  MapCluster *right_child_;
  unordered_map<Edge *, NodeSize> external_edges_;
  unordered_set<Edge *> internal_edges_;
  // compilation members
  const unordered_map<MapCluster *, SddLiteral> *cluster_variable_map_;
  const unordered_map<Edge *, SddLiteral> *edge_variable_map_;
  const unordered_map<Edge *, MapCluster *> *edge_cluster_map_;
  Vtree *local_vtree_;
  PsddManager *pm_;
  // compilation results for internal clusters
  PsddNode *internal_path_constraint_;
  PsddNode *empty_path_constraint_;
  std::map<Edge *, PsddNode *> terminal_path_constraint_;
  std::map<std::pair<Edge *, Edge *>, PsddNode *> non_terminal_path_constraint_;
};
#endif  // HIERARCHICAL_MAP_MAP_CLUSTER_H
