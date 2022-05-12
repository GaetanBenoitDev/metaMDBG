#ifndef ABPOA_GRAPH_H
#define ABPOA_GRAPH_H

#include <stdint.h>
#include "abpoa.h"
#include "utils.h"

//#define CIGAR_STR "MIDNSHP=XB"
//#define ABPOA_GRAPH_CIGAR_STR "=XIDNSH"
//#define ABPOA_GRAPH_CEQUAL 0
//#define ABPOA_GRAPH_CMISMATCH 1
//#define ABPOA_GRAPH_CINS 2
//#define ABPOA_GRAPH_CDEL 3
//#define ABPOA_GRAPH_CREF_SKIP 4
//#define ABPOA_GRAPH_CCLIP 5

#ifdef __cplusplus
extern "C" {
#endif

int abpoa_get_aligned_id(abpoa_graph_t *abg, int node_id, uint8_t base);
void abpoa_add_graph_aligned_node(abpoa_graph_t *abg, int node_id, int aligned_id);
void abpoa_set_msa_rank(abpoa_graph_t *abg, int src_id, int sink_id);
abpoa_graph_t *abpoa_init_graph(void);
void abpoa_free_graph(abpoa_graph_t *graph);

static inline int abpoa_graph_node_id_to_index(abpoa_graph_t *graph, int node_id) {
    if (node_id < 0 || node_id >= graph->node_n) err_fatal(__func__, "Wrong node id: %d\n", node_id);
    return graph->node_id_to_index[node_id];
}

static inline int abpoa_graph_node_id_to_max_pos_right(abpoa_graph_t *graph, int node_id) {
    if (node_id < 0 || node_id >= graph->node_n) err_fatal(__func__, "Wrong node id: %d\n", node_id);
    return graph->node_id_to_max_pos_right[node_id];
}

static inline int abpoa_graph_node_id_to_max_pos_left(abpoa_graph_t *graph, int node_id) {
    if (node_id < 0 || node_id >= graph->node_n) err_fatal(__func__, "Wrong node id: %d\n", node_id);
    return graph->node_id_to_max_pos_left[node_id];
}

static inline int abpoa_graph_node_id_to_max_remain(abpoa_graph_t *graph, int node_id) {
    if (node_id < 0 || node_id >= graph->node_n) err_fatal(__func__, "Wrong node id: %d\n", node_id);
    return graph->node_id_to_max_remain[node_id];
}

static inline int abpoa_graph_index_to_node_id(abpoa_graph_t *graph, int index_i) {
    if (index_i < 0 || index_i >= graph->node_n) err_fatal(__func__, "Wrong index: %d\n", index_i);
    return graph->index_to_node_id[index_i];
}

static inline int abpoa_graph_node_id_to_msa_rank(abpoa_graph_t *graph, int node_id) {
    if (node_id < 0 || node_id >= graph->node_n) err_fatal(__func__, "Wrong node id: %d\n", node_id);
    return graph->node_id_to_msa_rank[node_id];
}

#ifdef __cplusplus
}
#endif

#endif
