//
// Created by Jason Shen on 6/2/18.
//

#ifndef HIERARCHICAL_MAP_SDD_SIMPLE_PATH_COMPILER_H
#define HIERARCHICAL_MAP_SDD_SIMPLE_PATH_COMPILER_H
extern "C" {
#include <sdd/sddapi.h>
}
#include <hierarchical_map_compiler/graph.h>
#include <hierarchical_map_compiler/types.h>
#include <unordered_map>

class SddSimplePathCompiler {
 public:
  static SddSimplePathCompiler *GetSddSimplePathCompiler(
      const std::vector<Edge *> &edges, NodeSize src_node, NodeSize dst_node,
      SddManager *sdd_manager);
  virtual ~SddSimplePathCompiler() = default;
  virtual SddNode *Compile() = 0;
};

#endif  // HIERARCHICAL_MAP_SDD_SIMPLE_PATH_COMPILER_H
