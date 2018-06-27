//
// Created by Jason Shen on 5/21/18.
//

#ifndef HIERARCHICAL_MAP_GRAPH_H
#define HIERARCHICAL_MAP_GRAPH_H

#include <hierarchical_map_compiler/edge.h>
#include <nlohmann/json.hpp>
#include <set>
#include <unordered_map>
#include <vector>
extern "C" {
#include <sdd/sddapi.h>
}

using json = nlohmann::json;
class Graph {
 public:
  ~Graph();
  const std::vector<Edge*>& edges() const;
  static Graph* GraphFromJsonEdgeList(const json& edge_list);
  static Graph* GraphFromEdgeList(std::vector<Edge*> edge_list);
  static Graph* GraphFromStolenEdgeList(std::vector<Edge*> edge_list);
  std::vector<Edge*> GreedyEdgeOrder() const;
  std::unordered_map<NodeSize, std::set<NodeSize>> AdjacencyMap() const;
  std::map<std::pair<NodeSize, NodeSize>, std::vector<Edge*>> NodePairToEdges()
      const;
  std::unordered_map<NodeSize, EdgeSize> NodeDegrees() const;
  std::set<NodeSize> Vertices() const;
  NodeSize node_size() const { return node_size_; }
  json to_json() const;

 private:
  explicit Graph(std::vector<Edge*> edges, bool stolen_edges);
  std::vector<Edge*> edges_;
  NodeSize node_size_;
  bool stolen_edges_;
};
#endif  // HIERARCHICAL_MAP_GRAPH_H
