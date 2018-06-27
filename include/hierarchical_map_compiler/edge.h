//
// Created by Jason Shen on 5/21/18.
//

#ifndef HIERARCHICAL_MAP_EDGE_H
#define HIERARCHICAL_MAP_EDGE_H

#include <hierarchical_map_compiler/types.h>
#include <cstdint>
#include <nlohmann/json.hpp>
#include <string>

using nlohmann::json;
class Edge {
 public:
  Edge(std::string edge_name, NodeSize x_node_index, NodeSize y_node_index);
  NodeSize x_node_index() const;
  NodeSize y_node_index() const;
  NodeSize OtherEnd(NodeSize one_end) const;
  const std::string& edge_name() const;
  std::pair<NodeSize, NodeSize> OrderedNodePair() const;
  json to_json() const;

 private:
  std::string edge_name_;
  NodeSize x_node_index_;
  NodeSize y_node_index_;
};

#endif  // HIERARCHICAL_MAP_EDGE_H
