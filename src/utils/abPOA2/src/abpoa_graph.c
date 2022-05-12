#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "abpoa_align.h"
#include "abpoa_seq.h"
#include "simd_abpoa_align.h"
#include "kdq.h"

KDQ_INIT(int)
#define kdq_int_t kdq_t(int)

abpoa_node_t *abpoa_init_node(int n) {
    abpoa_node_t *node = (abpoa_node_t*)_err_calloc(n, sizeof(abpoa_node_t));
    return node;
}

void abpoa_set_graph_node(abpoa_graph_t *abg, int node_i) {
    abg->node[node_i].node_id = node_i;
    abg->node[node_i].in_edge_n = 0; abg->node[node_i].in_edge_m = 0;
    abg->node[node_i].out_edge_n = 0; abg->node[node_i].out_edge_m = 0;
    abg->node[node_i].aligned_node_n = 0; abg->node[node_i].aligned_node_m = 0;
    abg->node[node_i].read_ids_n = 0;
}

void abpoa_free_node(abpoa_node_t *node, int n) {
    int i, j;
    for (i = 0; i < n; ++i) {
        if (node[i].in_edge_m > 0) free(node[i].in_id);
        if (node[i].out_edge_m > 0) {
            free(node[i].out_id); free(node[i].out_weight);
            if (node[i].read_ids_n > 0) {
                for (j = 0; j < node[i].out_edge_m; ++j) {
                    free(node[i].read_ids[j]);
                } 
            } free(node[i].read_ids);
        }
        if (node[i].aligned_node_m > 0) free(node[i].aligned_node_id);
    }
    free(node);
}

// 0: in_edge, 1: out_edge
abpoa_graph_t *abpoa_realloc_graph_edge(abpoa_graph_t *abg, int io, int id, int use_read_ids) {
    if (io == 0) {
        _uni_realloc(abg->node[id].in_id, abg->node[id].in_edge_n, abg->node[id].in_edge_m, int);
    } else {
        int edge_m = abg->node[id].out_edge_m;
        if (edge_m <= 0) {
            abg->node[id].out_edge_m = MAX_OF_TWO(abg->node[id].out_edge_n, 1);
            abg->node[id].out_id = (int*)_err_malloc(abg->node[id].out_edge_m * sizeof(int));
            abg->node[id].out_weight = (int*)_err_malloc(abg->node[id].out_edge_m * sizeof(int));
            abg->node[id].read_ids = (uint64_t**)_err_malloc(abg->node[id].out_edge_m * sizeof(uint64_t*));
            if (abg->node[id].read_ids_n > 0) {
                int i;
                for (i = 0; i < abg->node[id].out_edge_m; ++i) {
                    abg->node[id].read_ids[i] = (uint64_t*)_err_calloc(abg->node[id].read_ids_n, sizeof(uint64_t));
                }
            }
        } else if (abg->node[id].out_edge_n >= edge_m) {
            abg->node[id].out_edge_m = abg->node[id].out_edge_n+1; kroundup32(abg->node[id].out_edge_m);
            abg->node[id].out_id = (int*)_err_realloc(abg->node[id].out_id, abg->node[id].out_edge_m * sizeof(int));
            abg->node[id].out_weight = (int*)_err_realloc(abg->node[id].out_weight, abg->node[id].out_edge_m * sizeof(int));
            abg->node[id].read_ids = (uint64_t**)_err_realloc(abg->node[id].read_ids, abg->node[id].out_edge_m * sizeof(uint64_t*));
            if (abg->node[id].read_ids_n > 0) {
                int i;
                for (i = edge_m; i < abg->node[id].out_edge_m; ++i) {
                    abg->node[id].read_ids[i] = (uint64_t*)_err_calloc(abg->node[id].read_ids_n, sizeof(uint64_t));
                }
            }
        }
    }
    return abg;
}

abpoa_graph_t *abpoa_realloc_graph_node(abpoa_graph_t *abg) {
    if (abg->node_m <= 0) {
        abg->node_m = 1;
        abg->node = (abpoa_node_t*)_err_calloc(1, sizeof(abpoa_node_t));
    }
    if (abg->node_n == abg->node_m) {
        int i;
        abg->node_m <<= 1;
        abg->node = (abpoa_node_t*)_err_realloc(abg->node, abg->node_m * sizeof(abpoa_node_t));
        for (i = abg->node_m >> 1; i < abg->node_m; ++i) {
            abpoa_set_graph_node(abg, i);
        }
    }
    return abg;
}

abpoa_graph_t *abpoa_init_graph(void) {
    abpoa_graph_t *abg = (abpoa_graph_t*)_err_malloc(sizeof(abpoa_graph_t));
    abg->node_n = 2, abg->node_m = 2, abg->index_rank_m = 0;
    abg->node = abpoa_init_node(2);
    abg->node[0].node_id = 0; abg->node[1].node_id = 1;
    abg->node[0].read_ids_n = 0; abg->node[1].read_ids_n = 0;
    abg->is_topological_sorted = abg->is_called_cons = 0;
    abg->node_id_to_index = NULL; abg->index_to_node_id = NULL; abg->node_id_to_msa_rank = NULL;
    abg->node_id_to_max_pos_left = NULL; abg->node_id_to_max_pos_right = NULL; abg->node_id_to_max_remain = NULL;
    return abg;
}

void abpoa_free_graph(abpoa_graph_t *abg) {
    if (abg->node_m > 0) abpoa_free_node(abg->node, abg->node_m);

    if (abg->node_n > 0) {
        free(abg->index_to_node_id);
        free(abg->node_id_to_index);
        if (abg->node_id_to_msa_rank) free(abg->node_id_to_msa_rank);

        if (abg->node_id_to_max_pos_left) free(abg->node_id_to_max_pos_left);
        if (abg->node_id_to_max_pos_right) free(abg->node_id_to_max_pos_right);
        if (abg->node_id_to_max_remain) free(abg->node_id_to_max_remain);
    }
    free(abg);
}

abpoa_cons_t *abpoa_init_cons(void) {
    abpoa_cons_t *abc = (abpoa_cons_t*)_err_malloc(sizeof(abpoa_cons_t));
    abc->n_cons = 0; abc->msa_len = 0;
    abc->clu_n_seq = NULL;
    abc->cons_len = NULL;
    abc->cons_node_ids = NULL;
    abc->cons_base = NULL;
    abc->msa_base = NULL;
    abc->cons_cov = NULL;
    abc->clu_read_ids = NULL;
    abc->cons_phred_score = NULL;
    return abc;
}

void abpoa_free_cons(abpoa_cons_t *abc) {
    int i;
    if (abc->n_cons > 0) {
        if (abc->clu_n_seq != NULL) free(abc->clu_n_seq);
        if (abc->cons_len != NULL) free(abc->cons_len);
        if (abc->cons_node_ids != NULL) {
            for (i = 0; i < abc->n_cons; ++i) free(abc->cons_node_ids[i]); free(abc->cons_node_ids);
        }
        if (abc->cons_base != NULL) {
            for (i = 0; i < abc->n_cons; ++i) free(abc->cons_base[i]); free(abc->cons_base);
        }
        if (abc->cons_cov != NULL) {
            for (i = 0; i < abc->n_cons; ++i) free(abc->cons_cov[i]); free(abc->cons_cov);
        }
        if (abc->clu_read_ids != NULL) {
            for (i = 0; i < abc->n_cons; ++i) free(abc->clu_read_ids[i]); free(abc->clu_read_ids);
        }
        if (abc->cons_phred_score != NULL) {
            for (i = 0; i < abc->n_cons; ++i) free(abc->cons_phred_score[i]); free(abc->cons_phred_score);
        }
    }
    if (abc->msa_len > 0) {
        if (abc->msa_base != NULL) {
            for (i = 0; i < abc->n_seq+abc->n_cons; ++i) free(abc->msa_base[i]);
            free(abc->msa_base);
        }
    }
    free(abc);
}

abpoa_t *abpoa_init(void) {
    abpoa_t *ab = (abpoa_t*)_err_malloc(sizeof(abpoa_t));
    ab->abg = abpoa_init_graph();
    ab->abs = abpoa_init_seq();
    ab->abm = abpoa_init_simd_matrix();
    ab->abc = abpoa_init_cons();
    return ab;
}

void abpoa_free(abpoa_t *ab) {
    abpoa_free_graph(ab->abg);
    abpoa_free_seq(ab->abs);
    abpoa_free_simd_matrix(ab->abm);
    abpoa_free_cons(ab->abc);
    free(ab);
}

void abpoa_BFS_set_node_index(abpoa_graph_t *abg, int src_id, int sink_id) {
    int *id, cur_id, out_id, aligned_id;
    int index = 0, q_size, new_q_size;

    int *in_degree = (int*)_err_malloc(abg->node_n * sizeof(int));
    int i, j;
    for (i = 0; i < abg->node_n; ++i) in_degree[i] = abg->node[i].in_edge_n;

    kdq_int_t *q = kdq_init_int();

    // Breadth-First-Search
    kdq_push_int(q, src_id); q_size = 1; new_q_size = 0; // node[q.id].in_degree equals 0
    while (q_size > 0) {
        if ((id = kdq_shift_int(q)) == 0) err_fatal_simple("Error in queue.");
        cur_id = *id;
        abg->index_to_node_id[index] = cur_id;
        abg->node_id_to_index[cur_id] = index++;

        if (cur_id == sink_id) {
            kdq_destroy_int(q); free(in_degree);
            return;
        }
        for (i = 0; i < abg->node[cur_id].out_edge_n; ++i) {
            out_id = abg->node[cur_id].out_id[i];
            if (--in_degree[out_id] == 0) {
                for (j = 0; j < abg->node[out_id].aligned_node_n; ++j) {
                    aligned_id = abg->node[out_id].aligned_node_id[j];
                    if (in_degree[aligned_id] != 0) goto next_out_node;
                }
                kdq_push_int(q, out_id);
                ++new_q_size;
                for (j = 0; j < abg->node[out_id].aligned_node_n; ++j) {
                    aligned_id = abg->node[out_id].aligned_node_id[j];
                    kdq_push_int(q, aligned_id);
                    ++new_q_size;
                }
            }
next_out_node:;
        }
        if (--q_size == 0) {
            q_size = new_q_size;
            new_q_size = 0;
        }
    }
    err_fatal_simple("Failed to set node index.");
}

void abpoa_BFS_set_node_remain(abpoa_graph_t *abg, int src_id, int sink_id) {
    int *id, cur_id, i, out_id, in_id;

    int *out_degree = (int*)_err_malloc(abg->node_n * sizeof(int));
    for (i = 0; i < abg->node_n; ++i) {
        out_degree[i] = abg->node[i].out_edge_n;
        abg->node_id_to_max_remain[i] = 0;
    }

    kdq_int_t *q = kdq_init_int();

    // Breadth-First-Search
    kdq_push_int(q, sink_id); // node[q.id].in_degree equals 0
    abg->node_id_to_max_remain[sink_id] = -1; // XXX not 0
    while ((id = kdq_shift_int(q)) != 0) {
        cur_id = *id;

        // all out_id of cur_id have beed visited
        // max weight out_id
        if (cur_id != sink_id) {
            int max_w=-1, max_id=-1;
            for (i = 0; i < abg->node[cur_id].out_edge_n; ++i) {
                out_id = abg->node[cur_id].out_id[i];
                if (abg->node[cur_id].out_weight[i] > max_w) {
                    max_w = abg->node[cur_id].out_weight[i];
                    max_id = out_id;
                }
            }
            abg->node_id_to_max_remain[cur_id] = abg->node_id_to_max_remain[max_id] + 1;
            // fprintf(stderr, "%d -> %d\n", abg->node_id_to_index[cur_id], abg->node_id_to_max_remain[cur_id]);
        }
        if (cur_id == src_id) {
            kdq_destroy_int(q); free(out_degree);
            return;
        }
        for (i = 0; i < abg->node[cur_id].in_edge_n; ++i) {
            in_id = abg->node[cur_id].in_id[i];
            if (--out_degree[in_id] == 0) kdq_push_int(q, in_id);
        }
    }
    err_fatal_simple("Failed to set node remain.");
}

// 1. index_to_node_id
// 2. node_id_to_index
// 3. node_id_to_rank
void abpoa_topological_sort(abpoa_graph_t *abg, abpoa_para_t *abpt) {
    if (abg->node_n <= 0) {
        err_func_format_printf(__func__, "Empty graph.\n");
        return;
    }
    int node_n = abg->node_n;
    if (node_n > abg->index_rank_m) {
        abg->index_rank_m = node_n; kroundup32(abg->index_rank_m);
        // fprintf(stderr, "node_n: %d, index_rank_m: %d\n", node_n, abg->index_rank_m);
        abg->index_to_node_id = (int*)_err_realloc(abg->index_to_node_id, abg->index_rank_m * sizeof(int));
        abg->node_id_to_index = (int*)_err_realloc(abg->node_id_to_index, abg->index_rank_m * sizeof(int));
        if (abpt->out_msa || abpt->max_n_cons > 1) 
            abg->node_id_to_msa_rank = (int*)_err_realloc(abg->node_id_to_msa_rank, abg->index_rank_m * sizeof(int));
        if (abpt->wb >= 0) {
            abg->node_id_to_max_pos_left = (int*)_err_realloc(abg->node_id_to_max_pos_left, abg->index_rank_m * sizeof(int));
            abg->node_id_to_max_pos_right = (int*)_err_realloc(abg->node_id_to_max_pos_right, abg->index_rank_m * sizeof(int));
            abg->node_id_to_max_remain = (int*)_err_realloc(abg->node_id_to_max_remain, abg->index_rank_m * sizeof(int));
        } else if (abpt->zdrop > 0) {
            abg->node_id_to_max_remain = (int*)_err_realloc(abg->node_id_to_max_remain, abg->index_rank_m * sizeof(int));
        }
    }
    // start from ABPOA_SRC_NODE_ID to ABPOA_SINK_NODE_ID
    abpoa_BFS_set_node_index(abg, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID);
    // init min/max rank
    if (abpt->wb >= 0) {
        int i;
        for (i = 0; i < node_n; ++i) {
            abg->node_id_to_max_pos_right[i] = 0;
            abg->node_id_to_max_pos_left[i] = node_n;
        }
        abpoa_BFS_set_node_remain(abg, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID);
    } else if (abpt->zdrop > 0)
        abpoa_BFS_set_node_remain(abg, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID);
    abg->is_topological_sorted = 1;
}

void abpoa_DFS_set_msa_rank(abpoa_graph_t *abg, int src_id, int sink_id, int *in_degree) {
    // fprintf(stderr, "node_n: %d, m: %d\n", abg->node_n, abg->index_rank_m);
    if (abg->node_n > abg->index_rank_m) {
        int m = abg->node_n; kroundup32(m);
        abg->node_id_to_msa_rank = (int*)_err_realloc(abg->node_id_to_msa_rank, m * sizeof(int));
    }
    int *id, cur_id, i, j, out_id, aligned_id;
    int msa_rank = 0;
    kdq_int_t *q = kdq_init_int();

    // Depth-First-Search
    kdq_push_int(q, src_id); // node[q.id].in_degree equals 0
    abg->node_id_to_msa_rank[src_id] = -1;
    // printf("tot_node_n: %d, node_m: %d\n", abg->node_n, abg->node_m);

    while((id = kdq_pop_int(q)) != 0) {
        cur_id = *id;
        if (abg->node_id_to_msa_rank[cur_id] < 0) {
            abg->node_id_to_msa_rank[cur_id] = msa_rank;
            for (i = 0; i < abg->node[cur_id].aligned_node_n; ++i) {
                aligned_id = abg->node[cur_id].aligned_node_id[i];
                abg->node_id_to_msa_rank[aligned_id] = msa_rank;
            }
            msa_rank++;
        }

        if (cur_id == sink_id) {
            kdq_destroy_int(q);
            abg->is_set_msa_rank = 1;
            return;
        }
        for (i = 0; i < abg->node[cur_id].out_edge_n; ++i) {
            out_id = abg->node[cur_id].out_id[i];
            if (--in_degree[out_id] == 0) {
                for (j = 0; j < abg->node[out_id].aligned_node_n; ++j) {
                    aligned_id = abg->node[out_id].aligned_node_id[j];
                    if (in_degree[aligned_id] != 0) goto next_out_node;
                }
                kdq_push_int(q, out_id);
                abg->node_id_to_msa_rank[out_id] = -1;
                for (j = 0; j < abg->node[out_id].aligned_node_n; ++j) {
                    aligned_id = abg->node[out_id].aligned_node_id[j];
                    kdq_push_int(q, aligned_id);
                    // printf("aln_id: %d\n", aligned_id);
                    abg->node_id_to_msa_rank[aligned_id] = -1;
                }
            }
next_out_node:;
        }
    }
    err_fatal_simple("Error in set_msa_rank.\n");
}

void abpoa_set_msa_rank(abpoa_graph_t *abg, int src_id, int sink_id) {
    if (abg->is_set_msa_rank == 0) {
        int i, *in_degree = (int*)_err_malloc(abg->node_n * sizeof(int));
        for (i = 0; i < abg->node_n; ++i) in_degree[i] = abg->node[i].in_edge_n;
        abpoa_DFS_set_msa_rank(abg, src_id, sink_id, in_degree);
        free(in_degree);
    }
}

int abpoa_get_aligned_id(abpoa_graph_t *abg, int node_id, uint8_t base) {
    int i, aln_id;
    abpoa_node_t *node = abg->node;
    for (i = 0; i < node[node_id].aligned_node_n; ++i) {
        aln_id = node[node_id].aligned_node_id[i];
        if (node[aln_id].base == base)
            return aln_id;
    }
    return -1;
}

void abpoa_add_graph_aligned_node1(abpoa_node_t *node, int aligned_id) {
    _uni_realloc(node->aligned_node_id, node->aligned_node_n, node->aligned_node_m, int);
    node->aligned_node_id[node->aligned_node_n++] = aligned_id;
}

void abpoa_add_graph_aligned_node(abpoa_graph_t *abg, int node_id, int aligned_id) {
    int i; abpoa_node_t *node = abg->node;
    for (i = 0; i < node[node_id].aligned_node_n; ++i) {
        abpoa_add_graph_aligned_node1(node + node[node_id].aligned_node_id[i], aligned_id);
        abpoa_add_graph_aligned_node1(node + aligned_id, node[node_id].aligned_node_id[i]);
    }
    abpoa_add_graph_aligned_node1(abg->node + node_id, aligned_id);
    abpoa_add_graph_aligned_node1(abg->node + aligned_id, node_id);
}

//TODO
void abpoa_set_read_id(uint64_t *read_ids, int read_id) {
    int n = read_id / 64;
    // if ((n = read_id / 64) >= read_ids_n) err_fatal("abpoa_set_read_id", "Unexpected read id: %d\n", read_id);
    uint64_t one = 1; int b = read_id & 0x3f;
    read_ids[n] |= (one << b);
}

int abpoa_add_graph_node(abpoa_graph_t *abg, uint8_t base) {
    int node_id = abg->node_n;
    abpoa_realloc_graph_node(abg);
    // add node
    abg->node[node_id].base = base;
    ++abg->node_n;
    return node_id;
}

int abpoa_add_graph_edge(abpoa_graph_t *abg, int from_id, int to_id, int check_edge, int w, uint8_t add_read_id, int read_id, int read_ids_n) {
    int ret = 1;
    if (from_id < 0 || from_id >= abg->node_n || to_id < 0 || to_id >= abg->node_n) err_fatal(__func__, "node_n: %d\tfrom_id: %d\tto_id: %d.", abg->node_n, from_id, to_id);
    int out_edge_n = abg->node[from_id].out_edge_n;
    int edge_exist = 0;
    int out_edge_i = -1;
    if (check_edge) {
        int i;
        for (i = 0; i < out_edge_n; ++i) {
            if (abg->node[from_id].out_id[i] == to_id) { // edge exists
                abg->node[from_id].out_weight[i] += w; // update weight on existing edge
                // update label id
                edge_exist = 1;
                out_edge_i = i;
                break;
            }
        }
    }

    // add edge
    if (edge_exist == 0) {
        /// in edge
        abpoa_realloc_graph_edge(abg, 0, to_id, 0);
        abg->node[to_id].in_id[abg->node[to_id].in_edge_n] = from_id;
        ++abg->node[to_id].in_edge_n;
        /// out edge
        abpoa_realloc_graph_edge(abg, 1, from_id, add_read_id);
        abg->node[from_id].out_id[out_edge_n] = to_id;
        abg->node[from_id].out_weight[out_edge_n] = w; // initial weight for new edge
        out_edge_i = out_edge_n;
        ++abg->node[from_id].out_edge_n;
    }
    
    // add read_id to out edge
    if (add_read_id) {
        if (out_edge_i < 0) err_fatal_simple("No edge found.");
        if (read_ids_n <= 0) err_fatal(__func__, "Unexpected read_ids_n: %d.", read_ids_n);
        int i, j;
        abpoa_node_t *from_node = abg->node + from_id;
        if (from_node->read_ids_n == 0) {
            for (i = 0; i < from_node->out_edge_m; ++i) {
                from_node->read_ids[i] = (uint64_t*)_err_calloc(read_ids_n, sizeof(uint64_t));
            }
            from_node->read_ids_n = read_ids_n;
        } else if (from_node->read_ids_n < read_ids_n) {
            // reallocate from_node->read_ids
            for (i = 0; i < from_node->out_edge_m; ++i) {
                from_node->read_ids[i] = (uint64_t*)_err_realloc(from_node->read_ids[i], read_ids_n * sizeof(uint64_t));
                for (j = from_node->read_ids_n; j < read_ids_n; ++j) from_node->read_ids[i][j] = 0;
            }
            from_node->read_ids_n = read_ids_n;
        }
        abpoa_set_read_id(from_node->read_ids[out_edge_i], read_id);
    }
    return ret;
}

void abpoa_add_graph_sequence(abpoa_graph_t *abg, uint8_t *seq, int seq_l, int *qpos_to_node_id, int start, int end, uint8_t add_read_id, int read_id, int read_ids_n) {
    if (start >= seq_l || end <= start) err_fatal(__func__, "seq_l: %d\tstart: %d\tend: %d.", seq_l, start, end);
    if (end > seq_l) end = seq_l;

    int i, last_node_id, cur_node_id;
    last_node_id = ABPOA_SRC_NODE_ID;
    for (i = start; i < end; ++i) {
        cur_node_id = abpoa_add_graph_node(abg, seq[i]);
        if (qpos_to_node_id) qpos_to_node_id[i] = cur_node_id;
        abpoa_add_graph_edge(abg, last_node_id, cur_node_id, 0, 1, add_read_id, read_id, read_ids_n);
        last_node_id = cur_node_id;
    }

    abpoa_add_graph_edge(abg, last_node_id, ABPOA_SINK_NODE_ID, 0, 1, add_read_id, read_id, read_ids_n);
    abg->is_called_cons = abg->is_set_msa_rank = abg->is_topological_sorted = 0;
    // abpoa_topological_sort(abg, abpt);
}

int is_full_upstream_subgraph(abpoa_graph_t *abg, int up_index, int down_index) {
    int i, j, id, in_id;
    for (i = up_index+1; i <= down_index; ++i) {
        id = abg->index_to_node_id[i];
        for (j = 0; j < abg->node[id].in_edge_n; ++j) {
            in_id = abg->node[id].in_id[j];
            if (abg->node_id_to_index[in_id] < up_index) return 0;
        }
    }
    return 1;
}

int abpoa_upstream_index(abpoa_graph_t *abg, int beg_index, int end_index) {
    int min_index, in_index, i, j, node_id, in_id;

    while (1) {
        min_index = beg_index;
        for (i = beg_index; i <= end_index; ++i) {
            node_id = abg->index_to_node_id[i];
            for (j = 0; j < abg->node[node_id].in_edge_n; ++j) {
                in_id = abg->node[node_id].in_id[j];
                in_index = abg->node_id_to_index[in_id];
                min_index = MIN_OF_TWO(min_index, in_index);
            }
        }
        if (is_full_upstream_subgraph(abg, min_index, beg_index)) {
            return min_index;
        } else {
            end_index = beg_index;
            beg_index = min_index; 
        }
    }
}

int is_full_downstream_subgraph(abpoa_graph_t *abg, int up_index, int down_index) {
    int i, j, id, out_id;
    for (i = up_index; i < down_index; ++i) {
        id = abg->index_to_node_id[i];
        for (j = 0; j < abg->node[id].out_edge_n; ++j) {
            out_id = abg->node[id].out_id[j];
            if (abg->node_id_to_index[out_id] > down_index) return 0;
        }
    }
    return 1;
}

int abpoa_downstream_index(abpoa_graph_t *abg, int beg_index, int end_index) {
    int max_index, out_index, i, j, node_id, out_id;

    while (1) {
        max_index = end_index;
        for (i = beg_index; i <= end_index; ++i) {
            node_id = abg->index_to_node_id[i];
            for (j = 0; j < abg->node[node_id].out_edge_n; ++j) {
                out_id = abg->node[node_id].out_id[j];
                out_index = abg->node_id_to_index[out_id];
                max_index = MAX_OF_TWO(max_index, out_index);
            }
        }
        if (is_full_upstream_subgraph(abg, end_index, max_index)) {
            return max_index;
        } else {
            beg_index = end_index;
            end_index = max_index;
        }
    }
}

//   exc_beg | inc_beg ... inc_end | exc_end
void abpoa_subgraph_nodes(abpoa_t *ab, abpoa_para_t *abpt, int exc_beg0, int exc_end0, int *exc_beg, int *exc_end) {
    abpoa_graph_t *abg = ab->abg;
    if (ab->abg->is_topological_sorted == 0) abpoa_topological_sort(abg, abpt);
    int inc_beg_index = abg->node_id_to_index[exc_beg0], inc_end_index = abg->node_id_to_index[exc_end0];

    int exc_beg_index = abpoa_upstream_index(abg, inc_beg_index, inc_end_index);
    int exc_end_index = abpoa_downstream_index(abg, inc_beg_index, inc_end_index);

    if (exc_beg_index < 0 || exc_end_index >= abg->node_n)
        err_fatal_simple("Error in subgraph_nodes");
    *exc_beg = abg->index_to_node_id[exc_beg_index];
    *exc_end = abg->index_to_node_id[exc_end_index];
}

// fusion stratergy :
// 1. Match: merge to one node
// 2. Mismatch: check if B is identical to A' aligned nodes, then merge to node; if not, add node
// 3. Insertion: add node
// 4. Deletion: nothing
// 5. Clipping: add node
// 6. For all first/last node, link to virtual start/end node

// inc_both_ends: set as 1 to add weight for edge between beg_node_id/end_node_id and internal node
int abpoa_add_subgraph_alignment(abpoa_t *ab, abpoa_para_t *abpt, int beg_node_id, int end_node_id, uint8_t *seq, int seq_l, int *qpos_to_node_id, abpoa_res_t res, int read_id, int tot_read_n, int inc_both_ends) {
    abpoa_graph_t *abg = ab->abg;
    int n_cigar = res.n_cigar; abpoa_cigar_t *abpoa_cigar = res.graph_cigar;
    int read_ids_n = 1 + ((tot_read_n-1) >> 6);
    uint8_t add_read_id = abpt->use_read_ids, add, w=1;
    if (abg->node_n == 2) { // empty graph
        abpoa_add_graph_sequence(abg, seq, seq_l, qpos_to_node_id, 0, seq_l, add_read_id, read_id, read_ids_n);
        return 0;
    } else {
        if (abg->node_n < 2) {
            err_fatal(__func__, "Graph node: %d.", abg->node_n);
        } else if (n_cigar == 0) {
            return 0;
            //err_fatal(__func__, "Empty graph cigar.");
        }
    }
    // normal graph, normal graph_cigar
    int i, j; int op, len;
    int node_id, query_id=-1, last_new = 0, last_id = beg_node_id, new_id, aligned_id;

    for (i = 0; i < n_cigar; ++i) {
        op = abpoa_cigar[i] & 0xf;
        if (op == ABPOA_CMATCH) {
            node_id = (abpoa_cigar[i] >> 34) & 0x3fffffff;
            query_id++; // = (abpoa_cigar[i] >> 4) & 0x3fffffff;
            if (abg->node[node_id].base != seq[query_id]) { // mismatch
                // check if query base is identical to node_id's aligned node
                if ((aligned_id = abpoa_get_aligned_id(abg, node_id, seq[query_id])) != -1) {
                    if (last_id != beg_node_id || inc_both_ends) add = 1; else add = 0;
                    abpoa_add_graph_edge(abg, last_id, aligned_id, 1-last_new, w, add_read_id&add, read_id, read_ids_n);
                    last_id = aligned_id; last_new = 0;
                } else {
                    new_id = abpoa_add_graph_node(abg, seq[query_id]);
                    if (last_id != beg_node_id || inc_both_ends) add = 1; else add = 0;
                    abpoa_add_graph_edge(abg, last_id, new_id, 0, w, add_read_id&add, read_id, read_ids_n);
                    last_id = new_id; last_new = 1;
                    // add new_id to node_id's aligned node
                    abpoa_add_graph_aligned_node(abg, node_id, new_id);
                }
            } else { // match
                if (last_id != beg_node_id || inc_both_ends) add = 1; else add = 0;
                abpoa_add_graph_edge(abg, last_id, node_id, 1-last_new, w, add_read_id&add, read_id, read_ids_n);
                last_id = node_id; last_new = 0;
            }
            if (qpos_to_node_id) qpos_to_node_id[query_id] = last_id;
        } else if (op == ABPOA_CINS || op == ABPOA_CSOFT_CLIP || op == ABPOA_CHARD_CLIP) {
            len = (abpoa_cigar[i] >> 4) & 0x3fffffff;
            query_id+=len; // = (abpoa_cigar[i] >> 34) & 0x3fffffff;
            for (j = len-1; j >= 0; --j) { // XXX use dynamic id, instead of static query_id
                new_id = abpoa_add_graph_node(abg, seq[query_id-j]);
                if (last_id != beg_node_id || inc_both_ends) add = 1; else add = 0;
                abpoa_add_graph_edge(abg, last_id, new_id, 0, w, add_read_id&add, read_id, read_ids_n);
                last_id = new_id; last_new = 1;
                if (qpos_to_node_id) qpos_to_node_id[query_id-j] = last_id;
            }
        } else if (op == ABPOA_CDEL) {
            // nothing;
            continue;
        }
    } 
    // if (inc_both_ends) add = 1; else add = 0; XXX end_node_id is always excluded when adding weight
    // abpoa_add_graph_edge(abg, last_id, end_node_id, 1-last_new, w, add_read_id&add, read_id, read_ids_n);
    abpoa_add_graph_edge(abg, last_id, end_node_id, 1-last_new, w, add_read_id, read_id, read_ids_n);
    abg->is_called_cons = abg->is_topological_sorted = 0;
    // abpoa_topological_sort(abg, abpt);
    return 0;
}

int abpoa_add_graph_alignment(abpoa_t *ab, abpoa_para_t *abpt, uint8_t *seq, int seq_l, int *qpos_to_node_id, abpoa_res_t res, int read_id, int tot_read_n, int inc_both_ends) {
    return abpoa_add_subgraph_alignment(ab, abpt, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, seq, seq_l, qpos_to_node_id, res, read_id, tot_read_n, inc_both_ends);
}

// reset allocated memery everytime init the graph
// * node
// * index_to_node_id/node_id_to_index/node_id_to_max_remain, max_pos_left/right
void abpoa_reset(abpoa_t *ab, abpoa_para_t *abpt, int qlen) {
    abpoa_graph_t *abg = ab->abg;
    int i, j, k, node_m;
    abg->is_topological_sorted = abg->is_called_cons = 0;
    for (i = 0; i < abg->node_n; ++i) {
        for (j = 0; j < abg->node[i].out_edge_n; ++j) {
            for (k = 0; k < abg->node[i].read_ids_n; ++k) abg->node[i].read_ids[j][k] = 0;
        }
        abg->node[i].in_edge_n = abg->node[i].out_edge_n = abg->node[i].aligned_node_n = 0;
            
    }
    abg->node_n = 2;
    if (qlen+2 > abg->node_m) {
        node_m = qlen+2; kroundup32(node_m);
        abg->node = (abpoa_node_t*)_err_realloc(abg->node, node_m * sizeof(abpoa_node_t));
        for (i = abg->node_m; i < node_m; ++i) 
            abpoa_set_graph_node(abg, i);
        abg->node_m = abg->index_rank_m = node_m;
        abg->index_to_node_id = (int*)_err_realloc(abg->index_to_node_id, node_m * sizeof(int));
        abg->node_id_to_index = (int*)_err_realloc(abg->node_id_to_index, node_m * sizeof(int));
        if (abpt->out_msa || abpt->max_n_cons > 1) 
            abg->node_id_to_msa_rank = (int*)_err_realloc(abg->node_id_to_msa_rank, node_m * sizeof(int));
        if (abpt->wb >= 0) {
            abg->node_id_to_max_pos_left = (int*)_err_realloc(abg->node_id_to_max_pos_left, node_m * sizeof(int));
            abg->node_id_to_max_pos_right = (int*)_err_realloc(abg->node_id_to_max_pos_right, node_m * sizeof(int));
            abg->node_id_to_max_remain = (int*)_err_realloc(abg->node_id_to_max_remain, node_m * sizeof(int));
        } else if (abpt->zdrop > 0) {
            abg->node_id_to_max_remain = (int*)_err_realloc(abg->node_id_to_max_remain, node_m * sizeof(int));
        }
    }
    // fprintf(stderr, "qlen: %d, node_n: %d, node_m: %d\n", qlen, abg->node_n, abg->node_m);
    // reset abs
    ab->abs->n_seq = 0;
    // reset cons
    abpoa_cons_t *abc = ab->abc;
    if (abc->n_cons > 0) {
        if (abc->clu_n_seq != NULL) free(abc->clu_n_seq);
        if (abc->cons_len != NULL) free(abc->cons_len);
        if (abc->cons_node_ids != NULL) {
            for (i = 0; i < abc->n_cons; ++i) free(abc->cons_node_ids[i]); free(abc->cons_node_ids);
        }
        if (abc->cons_base != NULL) {
            for (i = 0; i < abc->n_cons; ++i) free(abc->cons_base[i]); free(abc->cons_base);
        }
        if (abc->cons_cov != NULL) {
            for (i = 0; i < abc->n_cons; ++i) free(abc->cons_cov[i]); free(abc->cons_cov);
        }
        if (abc->clu_read_ids != NULL) {
            for (i = 0; i < abc->n_cons; ++i) free(abc->clu_read_ids[i]); free(abc->clu_read_ids);
        }
        if (abc->cons_phred_score != NULL) {
            for (i = 0; i < abc->n_cons; ++i) free(abc->cons_phred_score[i]); free(abc->cons_phred_score);
        }
    }
    if (abc->msa_len > 0) {
        if (abc->msa_base != NULL) {
            for (i = 0; i < abc->n_seq+abc->n_cons; ++i) free(abc->msa_base[i]);
            free(abc->msa_base);
        }
    }
    abc->n_seq = abc->n_cons = abc->msa_len = 0;
}
