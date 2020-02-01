# Top down compiler for binary hierarchical map.

This package implements the paper [Structured Bayesian Networks: From Inference to Learning with Routes](https://www.aaai.org/ojs/index.php/AAAI/article/view/4796). It compiles an SDD representing the underlying constraint of a binary hierarchical map. Given training data and the compiled SDD, one can learn a PSDD that induces a joint distribution over binary hierarchical routes, using packages [PSDD](https://github.com/hahaXD/psdd), [PyPSDD](https://github.com/art-ai/pypsdd).

### Dependency
* [GMP](https://gmplib.org/manual/C_002b_002b-Interface-General.html).
* [Graphillion](https://github.com/takemaru/graphillion).
* [CMake](https://cmake.org).

### Compile
```bash
mkdir build; cd build; cmake ..; make -j4;
```

### Test
```bash
./hmc_test <path_to_script_directory>
```

The <path_to_script_directory> is the path to the script directory of the project. If the project is compiled according to the previous section, it should be "../script"

### Run
```bash
./hmc_main <map_filename> <graph_hopper_script> <tmp_dir> <thread_num> <sdd_filename> <vtree_filename>
```
* map_filename: a json file representing the binary hierarchical map. Examples are given below.
* graph_hopper_script: the path to the python script that generates a logical circuit representing the simple path constraint. By default, it should be ../script/compile_graph.py.
* tmp_dir: a directory that stores temp files.
* thread_num: the number of parallel threads to compile the binary hierarchical map.
* sdd_filename: the filename for the resulting sdd circuit.
* vtree_filename: the filename for the vtree of the resulting sdd circuit.
..* To train a PSDD, which represents a joint distribution over all possible paths, one needs to provide the learner with both sdd_file and vtre_file. See [PSDD](https://github.com/hahaXD/psdd), [PyPSDD](https://github.com/art-ai/pypsdd).


### Map File Format
```json
 { "edge_size": 3,      
   "edges": [ {"__Edge__": true,  "x": 3, "name": "0", "y":0 }, {"__Edge__": true,  "x": 0, "name": "1", "y":1 } ,{"__Edge__": true,  "x": 1, "name": "2", "y":2 }], 
   "node_size": 4,
     "clusters": {       
         "left": { "nodes":[0, 3]},  
         "right": {"nodes":[1, 2]},
         "root": {  "sub_clusters": [ "left",  "right"] }
      }
  }
```
* edge_size: Number of edges in your graph
* edges: The edges in your graph. Field x and y are the node indexes that the edge connects.
* node_size: Number of nodes in your graph.
* clusters: A dictionary of hierarhical clusters, with the format <key:cluster_name>:<value: cluster_obj>.
..* For leaf clusters, the corresponding cluster_obj has a field called "nodes" that contains a list of node indexes in the leaf cluster.
..* For internal clusters, the corresponding cluster_obj has a field called "sub_clusters" that contains the name of the child clusters.
