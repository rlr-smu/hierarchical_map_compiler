//
// Created by Jason Shen on 6/2/18.
//

#ifndef HIERARCHICAL_MAP_SDD_SIMPLE_PATH_COMPILER_H
#define HIERARCHICAL_MAP_SDD_SIMPLE_PATH_COMPILER_H
extern "C" {
#include <sddapi.h>
}
#include <unordered_map>
#include "graph.h"
#include "types.h"

class SddSimplePathCompiler {
 public:
  static SddSimplePathCompiler *GetSddSimplePathCompiler(
      const std::vector<Edge *> &edges, NodeSize src_node, NodeSize dst_node,
      SddManager *sdd_manager);
  virtual ~SddSimplePathCompiler() = default;
  virtual SddNode *Compile() = 0;
};

#endif  // HIERARCHICAL_MAP_SDD_SIMPLE_PATH_COMPILER_H
