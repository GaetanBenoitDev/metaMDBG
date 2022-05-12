#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "abpoa.h"
#include "abpoa_graph.h"
#include "utils.h"
#include "abpoa_seq.h"
#include "kdq.h"

extern char ab_char256_table[256];
char ab_LogTable65536[65536];
char ab_bit_table16[65536];

#define NAT_E 2.718281828459045
static const char ab_LogTable256[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
    -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
    LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
    LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2_32(uint32_t v)
{
    uint32_t t, tt;
    if ((tt = v>>16)) return (t = tt>>8) ? 24 + ab_LogTable256[t] : 16 + ab_LogTable256[tt];
    return (t = v>>8) ? 8 + ab_LogTable256[t] : ab_LogTable256[v];
}

void set_65536_table(void) {
    int i;
    for (i = 0; i < 65536; ++i) {
        ab_LogTable65536[i] = ilog2_32(i);
    }
}

void set_bit_table16(void) {
    int i; ab_bit_table16[0] = 0;
    for (i = 0; i != 65536; ++i) ab_bit_table16[i] = (i&1) + ab_bit_table16[i>>1];
}

#define get_bit_cnt4(table, b) (table[(b)&0xffff] + table[(b)>>16&0xffff] + table[(b)>>32&0xffff] + table[(b)>>48&0xffff])

static inline int ilog2_64(uint64_t v) {
    uint64_t t, tt;
    if ((tt = v >> 32)) return (t = tt >> 16) ? 48 + ab_LogTable65536[t] : 32 + ab_LogTable65536[tt];
    return (t = v>>16) ? 16 + ab_LogTable65536[t] : ab_LogTable65536[v];
}

KDQ_INIT(int)
#define kdq_int_t kdq_t(int)

static inline int get_read_cnt(uint64_t *read_ids, int read_ids_n) {
    int i, c;
    for (i = c =0; i < read_ids_n; ++i) {
        c += get_bit_cnt4(ab_bit_table16, read_ids[i]);
    }
    return c;
}

abpoa_cons_t *abpoa_allocate_rc_msa(abpoa_cons_t *abc, int msa_len, int n_seq, int n_cons) {
    int i;
    abc->n_seq = n_seq; abc->msa_len = msa_len;
    abc->msa_base = (uint8_t**)_err_malloc((n_seq+n_cons) * sizeof(uint8_t*));
    for (i = 0; i < n_seq+n_cons; ++i) {
        abc->msa_base[i] = (uint8_t*)_err_malloc(msa_len * sizeof(uint8_t));
    }
    return abc;
}

void abpoa_output_rc_msa(abpoa_t *ab, abpoa_para_t *abpt, FILE *out_fp) {
    if (out_fp == NULL) return;
    int i, j;
    abpoa_seq_t *abs = ab->abs; abpoa_cons_t *abc = ab->abc;
    if (abc->msa_len <= 0) return;
    for (i = 0; i < abs->n_seq; ++i) {
        if (abs->name[i].l > 0) {
            if (abs->is_rc[i]) fprintf(out_fp, ">%s_reverse_complement\n", abs->name[i].s);
            else fprintf(out_fp, ">%s\n", abs->name[i].s);
        } else {
            fprintf(out_fp, ">Seq_%d\n", i+1);
        }
        for (j = 0; j < abc->msa_len; ++j) fprintf(out_fp, "%c", ab_char256_table[abc->msa_base[i][j]]);
        fprintf(out_fp, "\n");
    }
    if (abpt->out_cons) { // RC-MSA for consensus sequence
        int cons_i;
        for (cons_i = 0; cons_i < abc->n_cons; cons_i++) {
            fprintf(out_fp, ">Consensus_sequence");
            if (abc->n_cons > 1) {
                fprintf(out_fp, "_%d ", cons_i+1);
                for (j = 0; j < abc->clu_n_seq[cons_i]; ++j) { // cluter read_id
                    if (j != 0) fprintf(out_fp, ",");
                    fprintf(out_fp, "%d", abc->clu_read_ids[cons_i][j]);
                }
            }
            fprintf(out_fp, "\n");
            for (i = 0; i < abc->msa_len; ++i) fprintf(out_fp, "%c", ab_char256_table[abc->msa_base[abc->n_seq+cons_i][i]]);
            fprintf(out_fp, "\n");
        }
    }
}

void abpoa_set_msa_seq(abpoa_node_t node, int rank, uint8_t **msa_base) {
    int i, j, b, read_id; uint8_t base = node.base;
    uint64_t num, tmp;

    b = 0;
    for (i = 0; i < node.read_ids_n; ++i) {
        for (j = 0; j < node.out_edge_n; ++j) {
            num = node.read_ids[j][i];
            while (num) {
                tmp = num & -num;
                read_id = ilog2_64(tmp);
                msa_base[b+read_id][rank-1] = base;
                num ^= tmp;
            }
        }
        b += 64;
    }
}

// only generate rc-msa, output in separated func
void abpoa_generate_rc_msa(abpoa_t *ab, abpoa_para_t *abpt) {
    abpoa_graph_t *abg = ab->abg;
    if (abg->node_n <= 2) return;
    abpoa_set_msa_rank(abg, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID);
    if (abpt->out_cons) abpoa_generate_consensus(ab, abpt);

    abpoa_seq_t *abs = ab->abs; abpoa_cons_t *abc = ab->abc;
    int i, j, aligned_id, n_seq = abs->n_seq;
    int msa_len = abg->node_id_to_msa_rank[ABPOA_SINK_NODE_ID]-1;

    abpoa_allocate_rc_msa(abc, msa_len, n_seq, abc->n_cons);
    for (i = 0; i < n_seq; ++i) {
        for (j = 0; j < abc->msa_len; ++j) 
            abc->msa_base[i][j] = abpt->m;
    }

    int rank;
    // if (out_fp && abpt->out_msa_header == 0) fprintf(out_fp, ">Multiple_sequence_alignment\n");
    for (i = 2; i < abg->node_n; ++i) {
        // get msa rank
        rank = abpoa_graph_node_id_to_msa_rank(abg, i);
        for (j = 0; j < abg->node[i].aligned_node_n; ++j) {
            aligned_id = abg->node[i].aligned_node_id[j];
            rank = MAX_OF_TWO(rank, abpoa_graph_node_id_to_msa_rank(abg, aligned_id));
        }
        // assign seq
        abpoa_set_msa_seq(abg->node[i], rank, abc->msa_base);
    }
    if (abpt->out_cons) {
        int cons_i, cur_id;
        for (cons_i = 0; cons_i < abc->n_cons; cons_i++) {
            for (i = 0; i < msa_len; ++i) abc->msa_base[n_seq+cons_i][i] = abpt->m;
            for (i = 0; i < abc->cons_len[cons_i]; ++i) {
                cur_id = abc->cons_node_ids[cons_i][i];
                rank = abpoa_graph_node_id_to_msa_rank(abg, cur_id);
                for (j = 0; j < abg->node[cur_id].aligned_node_n; ++j) {
                    aligned_id = abg->node[cur_id].aligned_node_id[j];
                    rank = MAX_OF_TWO(rank, abpoa_graph_node_id_to_msa_rank(abg, aligned_id));
                }
                abc->msa_base[n_seq+cons_i][rank-1] = abc->cons_base[cons_i][i];
            }
        }
    }
}

// generate & output gfa
void abpoa_generate_gfa(abpoa_t *ab, abpoa_para_t *abpt, FILE *out_fp) {
    if (out_fp == NULL) return;
    abpoa_seq_t *abs = ab->abs; abpoa_graph_t *abg = ab->abg;
    if (abg->node_n <= 2) return;

    // traverse graph 
    int *in_degree = (int*)_err_malloc(abg->node_n * sizeof(int));
    int n_seq = abs->n_seq;
    int **read_paths = (int**)_err_malloc(n_seq * sizeof(int*)), *read_path_i = (int*)_err_calloc(n_seq, sizeof(int));
    int i, j, cur_id, pre_id, out_id, *id;
    for (i = 0; i < abg->node_n; ++i) in_degree[i] = abg->node[i].in_edge_n;
    for (i = 0; i < n_seq; ++i) read_paths[i] = (int*)_err_malloc(abg->node_n * sizeof(int));

    // output comment and header
    int nl = 0;
    for (i = 2; i < abg->node_n; ++i) nl += abg->node[i].in_edge_n;
    fprintf(out_fp, "H\tVN:Z:1.0\tNS:i:%d\tNL:i:%d\tNP:i:%d\n", abg->node_n-2, nl - abg->node[ABPOA_SRC_NODE_ID].out_edge_n, n_seq + abpt->out_cons);

    kdq_int_t *q = kdq_init_int();

    // Breadth-First-Search
    kdq_push_int(q, ABPOA_SRC_NODE_ID); 
    while ((id = kdq_shift_int(q)) != 0) {
        cur_id = *id;
        if (cur_id == ABPOA_SINK_NODE_ID) {
            kdq_destroy_int(q);
            break;
        } else {
            if (cur_id != ABPOA_SRC_NODE_ID) {
                // output node
                fprintf(out_fp, "S\t%d\t%c\n", cur_id-1, ab_char256_table[abg->node[cur_id].base]);
                // output all links based pre_ids
                for (i = 0; i < abg->node[cur_id].in_edge_n; ++i) {
                    pre_id = abg->node[cur_id].in_id[i];
                    if (pre_id != ABPOA_SRC_NODE_ID)
                        fprintf(out_fp, "L\t%d\t+\t%d\t+\t0M\n", pre_id-1, cur_id-1);
                }
                // add node id to read path
                int b, read_id; uint64_t num, tmp;
                b = 0;
                for (i = 0; i < abg->node[cur_id].read_ids_n; ++i) {
                    for (j = 0; j < abg->node[cur_id].out_edge_n; ++j) {
                        num = abg->node[cur_id].read_ids[j][i];
                        while (num) {
                            tmp = num & -num;
                            read_id = ilog2_64(tmp);
                            read_paths[b+read_id][read_path_i[b+read_id]++] = cur_id-1;
                            num ^= tmp;
                        }
                    }
                    b += 64;
                }
            }
            for (i = 0; i < abg->node[cur_id].out_edge_n; ++i) {
                out_id = abg->node[cur_id].out_id[i];
                if (--in_degree[out_id] == 0) {
                    kdq_push_int(q, out_id);
                }
            }
        }
    }
    // output read paths
    for (i = 0; i < n_seq; ++i) {
        if (abs->name[i].l > 0) fprintf(out_fp, "P\t%s\t", abs->name[i].s);
        else fprintf(out_fp, "P\t%d\t", i+1);
        if (abs->is_rc[i]) {
            for (j = read_path_i[i]-1; j >= 0; --j) {
                fprintf(out_fp, "%d-", read_paths[i][j]);
                if (j != 0) fprintf(out_fp, ",");
                else fprintf(out_fp, "\t*\n");
            }
        } else {
            for (j = 0; j < read_path_i[i]; ++j) {
                fprintf(out_fp, "%d+", read_paths[i][j]);
                if (j != read_path_i[i]-1) fprintf(out_fp, ",");
                else fprintf(out_fp, "\t*\n");
            }
        }
    }
    if (abpt->out_cons) {
        abpoa_generate_consensus(ab, abpt);
        abpoa_cons_t *abc = ab->abc;
        int cons_i;
        for (cons_i = 0; cons_i < abc->n_cons; ++cons_i) {
            fprintf(out_fp, "P\tConsensus_sequence");
            if (abc->n_cons > 1) fprintf(out_fp, "_%d", cons_i+1);
            fprintf(out_fp, "\t");
            for (i = 0; i < abc->cons_len[cons_i]; ++i) {
                cur_id = abc->cons_node_ids[cons_i][i];
                fprintf(out_fp, "%d+", cur_id-1);
                if (i != abc->cons_len[cons_i]-1) fprintf(out_fp, ",");
                else fprintf(out_fp, "\t*\n"); 
            }

        }
    }
    free(in_degree);
    for (i = 0; i < n_seq; ++i) free(read_paths[i]); 
    free(read_paths); free(read_path_i);
}

int abpoa_cons_phred_score(int n_cov, int n_seq) {
    if (n_cov > n_seq) err_fatal(__func__, "Error: unexpected n_cov/n_seq (%d/%d).", n_cov, n_seq);
    double x, p;
    x = 13.8 * (1.25 * n_cov / n_seq - 0.25);
    p = 1 - 1.0 / (1.0 + pow(NAT_E, -1 * x));
    return (33 + (int)(-10 * log10(p) + 0.499));
}

int get_read_ids_clu_weight(uint64_t *cur_read_ids, int read_ids_n, uint64_t *clu_read_ids) {
    int w = 0, i; uint64_t b;
    for (i = 0; i < read_ids_n; ++i) {
        b = cur_read_ids[i] & clu_read_ids[i];
        w += get_bit_cnt4(ab_bit_table16, b);
    }
    return w;
}

int abpoa_consensus_cov(abpoa_graph_t *abg, int id, uint64_t *clu_read_ids) {
    int i, j, in_id, left_w, right_w;
    // for each id: get max{left_weigth, right_weight}
    left_w = right_w = 0;
    for (i = 0; i < abg->node[id].in_edge_n; ++i) {
        in_id = abg->node[id].in_id[i];
        for (j = 0; j < abg->node[in_id].out_edge_n; ++j) {
            if (abg->node[in_id].out_id[j] == id) {
                left_w += get_read_ids_clu_weight(abg->node[in_id].read_ids[j], abg->node[in_id].read_ids_n, clu_read_ids);
                break;
            }
        }
    }
    for (i = 0; i < abg->node[id].out_edge_n; ++i) {
        right_w += get_read_ids_clu_weight(abg->node[id].read_ids[i], abg->node[id].read_ids_n, clu_read_ids);
    }
    return MAX_OF_TWO(left_w, right_w);
}

void abpoa_set_hb_cons(abpoa_graph_t *abg, int **max_out_id, int n_cons, uint64_t **clu_read_ids, int src_id, int sink_id, abpoa_cons_t *abc) {
    abc->n_cons = n_cons;
    int i, j, cur_id;
    for (i = 0; i < n_cons; ++i) {
        cur_id = max_out_id[i][src_id];
        j = 0;
        while (cur_id != sink_id) {
            abc->cons_node_ids[i][j] = cur_id;
            abc->cons_base[i][j] = abg->node[cur_id].base;
            abc->cons_cov[i][j] = abpoa_consensus_cov(abg, cur_id, clu_read_ids[i]);
            abc->cons_phred_score[i][j] = abpoa_cons_phred_score(abc->cons_cov[i][j], abc->clu_n_seq[i]);
            ++j;
            cur_id = max_out_id[i][cur_id];
        }
        abc->cons_len[i] = j;
    }
}

int abpoa_consensus_cov1(abpoa_graph_t *abg, int id) {
    int i, j, in_id, left_w, right_w;
    // for each id: get max{left_weigth, right_weight}
    left_w = right_w = 0;
    for (i = 0; i < abg->node[id].in_edge_n; ++i) {
        in_id = abg->node[id].in_id[i];
        for (j = 0; j < abg->node[in_id].out_edge_n; ++j) {
            if (abg->node[in_id].out_id[j] == id) {
                left_w += abg->node[in_id].out_weight[j];
                break;
            }
        }
    }
    for (i = 0; i < abg->node[id].out_edge_n; ++i) {
        right_w += abg->node[id].out_weight[i];
    }
    return MAX_OF_TWO(left_w, right_w);
}

void abpoa_set_hb_cons1(abpoa_graph_t *abg, int *max_out_id, int src_id, int sink_id, abpoa_cons_t *abc) {
    int i = 0, cur_id; 
    abc->n_cons = 1;
    cur_id = max_out_id[src_id];
    while (cur_id != sink_id) {
        abc->cons_node_ids[0][i] = cur_id;
        abc->cons_base[0][i] = abg->node[cur_id].base;
        abc->cons_cov[0][i] = abpoa_consensus_cov1(abg, cur_id);
        abc->cons_phred_score[0][i] = abpoa_cons_phred_score(abc->cons_cov[0][i], abc->n_seq);
        cur_id = max_out_id[cur_id];
        ++i;
    }
    abc->cons_len[0] = i;
}

// heaviest_bundling
// 1. argmax{cur->weight}
// 2. argmax{out_node->weight}
void abpoa_heaviest_bundling(abpoa_graph_t *abg, int src_id, int sink_id, int *out_degree, abpoa_cons_t *abc) {
    int *id, i, cur_id, in_id, out_id, max_id; int max_w, out_w;
    int *score = (int*)_err_malloc(abg->node_n * sizeof(int));
    int *max_out_id = (int*)_err_malloc(abg->node_n * sizeof(int));

    kdq_int_t *q = kdq_init_int();
    kdq_push_int(q, sink_id);
    // reverse Breadth-First-Search
    while ((id = kdq_shift_int(q)) != 0) {
        cur_id = *id;
        if (cur_id == sink_id) {
            max_out_id[cur_id] = -1;
            score[cur_id] = 0;
        } else {
            max_id = -1;
            if (cur_id == src_id) {
                int path_score = -1, path_max_w = -1;
                for (i = 0; i < abg->node[cur_id].out_edge_n; ++i) {
                    out_id = abg->node[cur_id].out_id[i];
                    out_w = abg->node[cur_id].out_weight[i];
                    if (out_w > path_max_w || (out_w == path_max_w && score[out_id] > path_score)) {
                        max_id = out_id;
                        path_score = score[out_id];
                        path_max_w = out_w;
                    }
                }
                max_out_id[cur_id] = max_id;
                kdq_destroy_int(q);
                break;
            } else {
                max_w = INT32_MIN;
                for (i = 0; i < abg->node[cur_id].out_edge_n; ++i) {
                    out_id = abg->node[cur_id].out_id[i];
                    out_w = abg->node[cur_id].out_weight[i];
                    if (max_w < out_w) {
                        max_w = out_w; max_id = out_id;
                    } else if (max_w == out_w && score[max_id] <= score[out_id]) {
                        max_id = out_id;
                    }
                }
                score[cur_id] = max_w + score[max_id];
                max_out_id[cur_id] = max_id;
            }
        }
        for (i = 0; i < abg->node[cur_id].in_edge_n; ++i) {
            in_id = abg->node[cur_id].in_id[i];
            if (--out_degree[in_id] == 0) kdq_push_int(q, in_id);
        }
    }
    abc->clu_n_seq[0] = abc->n_seq;
    // set cons read ids
    for (i = 0; i < abc->n_seq; ++i) abc->clu_read_ids[0][i] = i;
    abpoa_set_hb_cons1(abg, max_out_id, src_id, sink_id, abc);
    free(score); free(max_out_id);
}

void set_clu_read_ids(abpoa_cons_t *abc, uint64_t **read_ids, int cons_i, int n_seq) {
    int n, i, j; uint64_t b, one = 1;
    for (i = n = 0; i < n_seq; ++i) {
        j = i / 64; b = i & 0x3f;
        if (read_ids[cons_i][j] & (one << b)) {
            abc->clu_read_ids[cons_i][n++] = i;
        }
    }
    if (n != abc->clu_n_seq[cons_i])
        err_fatal(__func__, "Error in set cluster read ids. (%d, %d)", n, abc->clu_n_seq[cons_i]);
}

void abpoa_multip_heaviest_bundling(abpoa_graph_t *abg, int src_id, int sink_id, int *out_degree, int n_clu, int read_ids_n, uint64_t **clu_read_ids, abpoa_cons_t *abc) {
    int *id, cons_i, i, cur_id, in_id, out_id, max_id; int max_w, out_w;

    int *_out_degree = (int*)_err_malloc(abg->node_n * sizeof(int));
    int *score = (int*)_err_malloc(abg->node_n * sizeof(int));
    int **max_out_id = (int**)_err_malloc(n_clu * sizeof(int*));
    for (i = 0; i < n_clu; ++i) max_out_id[i] = (int*)_err_malloc(abg->node_n * sizeof(int));

    for (cons_i = 0; cons_i < n_clu; cons_i++) {
        for (i = 0; i < abg->node_n; ++i) _out_degree[i] = out_degree[i];
        abc->clu_n_seq[cons_i] = get_read_cnt(clu_read_ids[cons_i], read_ids_n);
        set_clu_read_ids(abc, clu_read_ids, cons_i, abc->n_seq);
        kdq_int_t *q = kdq_init_int();
        kdq_push_int(q, sink_id);
        // reverse Breadth-First-Search
        while ((id = kdq_shift_int(q)) != 0) {
            cur_id = *id;
            if (cur_id == sink_id) {
                max_out_id[cons_i][cur_id] = -1;
                score[cur_id] = 0;
            } else {
                max_id = -1;
                if (cur_id == src_id) {
                    int path_score = -1, path_max_w = -1;
                    for (i = 0; i < abg->node[cur_id].out_edge_n; ++i) {
                        out_id = abg->node[cur_id].out_id[i];
                        out_w = get_read_ids_clu_weight(abg->node[cur_id].read_ids[i], abg->node[cur_id].read_ids_n, clu_read_ids[cons_i]);
                        out_w = abg->node[cur_id].out_weight[i];
                        if (out_w > path_max_w || (out_w == path_max_w && score[out_id] > path_score)) {
                            max_id = out_id;
                            path_score = score[out_id];
                            path_max_w = out_w;
                        }
                    }
                    max_out_id[cons_i][cur_id] = max_id;
                    kdq_destroy_int(q);
                    break;
                } else {
                    max_w = INT32_MIN;
                    for (i = 0; i < abg->node[cur_id].out_edge_n; ++i) {
                        out_id = abg->node[cur_id].out_id[i];
                        out_w = get_read_ids_clu_weight(abg->node[cur_id].read_ids[i], abg->node[cur_id].read_ids_n, clu_read_ids[cons_i]);
                        if (max_w < out_w) {
                            max_w = out_w;
                            max_id = out_id;
                        } else if (max_w == out_w && score[max_id] <= score[out_id]) {
                            max_id = out_id;
                        }
                    }
                    score[cur_id] = max_w + score[max_id];
                    max_out_id[cons_i][cur_id] = max_id;
                }
            }
            for (i = 0; i < abg->node[cur_id].in_edge_n; ++i) {
                in_id = abg->node[cur_id].in_id[i];
                if (--_out_degree[in_id] == 0) 
                    kdq_push_int(q, in_id);
            }
        }
    }
    abpoa_set_hb_cons(abg, max_out_id, n_clu, clu_read_ids, src_id, sink_id, abc);

    free(score); free(_out_degree);
    for (i = 0; i < n_clu; ++i) free(max_out_id[i]); free(max_out_id);
}

void abpoa_output_fx_consensus(abpoa_t *ab, abpoa_para_t *abpt, FILE *out_fp) {
    if (out_fp == NULL) return;
    int cons_i, j;
    abpoa_cons_t *abc = ab->abc;
    for (cons_i = 0; cons_i < abc->n_cons; ++cons_i) {
        if (abpt->out_fq) fprintf(out_fp, "@Consensus_sequence");
        else fprintf(out_fp, ">Consensus_sequence");
        if (abc->n_cons > 1) {
            fprintf(out_fp, "_%d ", cons_i+1); // cons_id
            for (j = 0; j < abc->clu_n_seq[cons_i]; ++j) { // cluter read_id
                if (j != 0) fprintf(out_fp, ",");
                fprintf(out_fp, "%d", abc->clu_read_ids[cons_i][j]);
            }
        }
        fprintf(out_fp, "\n");
        for (j = 0; j < abc->cons_len[cons_i]; ++j) {
            fprintf(out_fp, "%c", ab_char256_table[abc->cons_base[cons_i][j]]);
        } fprintf(out_fp, "\n");
        if (abpt->out_fq) {
            fprintf(out_fp, "+Consensus_sequence");
            if (abc->n_cons > 1) {
                fprintf(out_fp, "_%d ", cons_i+1); // cons_id
                for (j = 0; j < abc->clu_n_seq[cons_i]; ++j) { // cluter read_id
                    if (j != 0) fprintf(out_fp, ",");
                    fprintf(out_fp, "%d", abc->clu_read_ids[cons_i][j]);
                }
            }
            fprintf(out_fp, "\n");
            for (j = 0; j < abc->cons_len[cons_i]; ++j) {
                fprintf(out_fp, "%c", abc->cons_phred_score[cons_i][j]);
            } fprintf(out_fp, "\n");
        }
    }
}

abpoa_cons_t *abpoa_allocate_cons(abpoa_cons_t *abc, int n_node, int n_seq, int n_cons) {
    int i;
    abc->n_cons = n_cons, abc->n_seq = n_seq;
    abc->clu_n_seq = (int*)_err_calloc(n_cons, sizeof(int));
    abc->cons_len = (int*)_err_calloc(n_cons, sizeof(int));
    abc->cons_node_ids = (int**)_err_malloc(n_cons * sizeof(int*));
    abc->cons_base = (uint8_t**)_err_malloc(n_cons * sizeof(uint8_t*));
    abc->cons_cov = (int**)_err_malloc(n_cons * sizeof(int*));
    abc->clu_read_ids = (int**)_err_malloc(n_cons * sizeof(int*));
    abc->cons_phred_score = (int**)_err_malloc(n_cons * sizeof(int*));
    for (i = 0; i < n_cons; ++i) {
        abc->cons_node_ids[i] = (int*)_err_malloc(n_node * sizeof(int));
        abc->cons_base[i] = (uint8_t*)_err_malloc(n_node * sizeof(uint8_t));
        abc->cons_cov[i] = (int*)_err_malloc(n_node * sizeof(int));
        abc->clu_read_ids[i] = (int*)_err_malloc(n_seq * sizeof(int));
        abc->cons_phred_score[i] = (int*)_err_malloc(n_node * sizeof(int));
    }
    return abc;
}

int abpoa_check_iden_read_ids(int **rc_weight, uint64_t ***read_ids, int m, int read_ids_n, int pos1, int pos2) {
    int i, j, k, iden = 1;
    uint8_t *map = (uint8_t*)_err_calloc(m, sizeof(uint8_t));

    for (i = 0; i < m ; ++i) {
        if (rc_weight[pos1][i] == 0) continue;
        int found_iden = 0;
        for (j = 0; j < m; ++j) { // find from 0~m that is identical to i'th read_ids
            if (map[j] == 1 || rc_weight[pos1][i] != rc_weight[pos2][j]) continue;
            // iden rc_weight
            int diff = 0;
            for (k = 0; k < read_ids_n; ++k) {
                if (read_ids[pos1][i][k] != read_ids[pos2][j][k]) {
                    diff = 1; break;
                }
            }
            if (diff == 0) { // i is identical to j
                found_iden = 1; 
                map[j] = 1;
                break;
            }
        }
        if (found_iden == 0) { // no iden for i'th base
            iden = 0; break;
        }
    }
    free(map);
    return iden;
}

// return: 1 if redundent else 0
int check_redundent_hap(int **clu_haps, int *clu_size, uint64_t **clu_read_ids, int n_clu, int new_clu_i, int n_het_pos, int read_id_i, uint64_t read_id) {
    int i, j, redundent = 0;
    for (i = n_clu-1; i >= 0; --i) {
        int iden = 1;
        for (j = 0; j < n_het_pos; ++j) {
            if (clu_haps[i][j] != clu_haps[new_clu_i][j]) {
                iden = 0; break;
            }
        }
        if (iden == 1) {
            clu_size[i] += 1;
            clu_read_ids[i][read_id_i] |= read_id;
            redundent = 1; break;
        }
    }
    if (redundent == 0) {
        clu_size[new_clu_i] += 1;
        clu_read_ids[new_clu_i][read_id_i] |= read_id;
    }
    return redundent;
}

int reassign_hap_by_min_w(int **clu_haps, int *clu_size, uint64_t **clu_read_ids, int read_ids_n, int n_clu, int min_w, int n_het_pos) {
    int i, j, k, n_reassign = 0;
    for (i = 0; i < n_clu; ++i) {
        if (clu_size[i] >= min_w || clu_size[i] == 0) continue;
        int reassign_i = -1, max_iden_pos = 0;
        for (j = 0; j < n_clu; ++j) {
            int n_iden_pos = 0;
            if (clu_size[j] < min_w) continue;
            // i < min_w, j >= min_w
            for (k = 0; k < n_het_pos; ++k) {
                if (clu_haps[i][k] == clu_haps[j][k]) n_iden_pos++;
            }
            if (n_iden_pos > max_iden_pos) {
                max_iden_pos = n_iden_pos;
                reassign_i = j;
            }
        }
        if (reassign_i >= 0) {
            for (j = 0; j < read_ids_n; ++j) {
                clu_read_ids[reassign_i][j] |= clu_read_ids[i][j];
                clu_read_ids[i][j] = 0;
            }
            clu_size[reassign_i] += clu_size[i];
            clu_size[i] = 0;
            n_reassign += 1;
        }
    }
    return n_clu - n_reassign;
}

int reassign_max_n_hap1(int **clu_haps, int *clu_size, uint64_t **clu_read_ids, int read_ids_n, int n_clu, int *clu_poss, int max_n_cons, int n_het_pos) {
    int i, j, k, n_reassign = 0;
    for (i = 0; i < n_clu; ++i) {
        int is_clu = 0;
        if (clu_size[i] == 0) continue;
        for (j = 0; j < max_n_cons; ++j) {
            if (i == clu_poss[j]) {
                is_clu = 1;
                break;
            }
        }
        if (is_clu) continue;

        int reassign_i = -1, max_iden_pos = 0;
        for (j = 0; j < max_n_cons; ++j) {
            int clu_i = clu_poss[j], n_iden_pos = 0;
            // i < min_w, clu_i >= min_w
            for (k = 0; k < n_het_pos; ++k) {
                if (clu_haps[i][k] == clu_haps[clu_i][k]) n_iden_pos++;
            }
            if (n_iden_pos > max_iden_pos) {
                max_iden_pos = n_iden_pos;
                reassign_i = clu_i;
            }
        }
        if (reassign_i >= 0) {
            for (j = 0; j < read_ids_n; ++j) {
                clu_read_ids[reassign_i][j] |= clu_read_ids[i][j];
                clu_read_ids[i][j] = 0;
            }
            clu_size[reassign_i] += clu_size[i];
            clu_size[i] = 0;
            n_reassign += 1;
        } else {
            clu_size[i] = 0;
        }
    }
    return n_clu - n_reassign;
}

typedef struct {
    int size, pos;
} clu_hap_tuple_t;

// descending order
int tup_cmpfunc (const void * a, const void * b) {
    return -(((clu_hap_tuple_t*)a)->size - ((clu_hap_tuple_t*)b)->size);
}

int reassign_max_n_hap(int **clu_haps, int *clu_size, uint64_t **clu_read_ids, int read_ids_n, int n_clu, int min_w, int n_het_pos, int max_n_cons) {
    int i;
    clu_hap_tuple_t *tup = (clu_hap_tuple_t*)_err_malloc(n_clu * sizeof(clu_hap_tuple_t));
    int *clu_poss = (int*)_err_malloc(max_n_cons * sizeof(int));

    while (n_clu > max_n_cons) {
        for (i = 0; i < n_clu; ++i) {
            tup[i].size = clu_size[i];
            tup[i].pos = i;
        }
        qsort(tup, n_clu, sizeof(clu_hap_tuple_t), tup_cmpfunc);
        // new min_w
        for (i = 0; i < max_n_cons; ++i) clu_poss[i] = tup[i].pos;
        int new_n_clu = reassign_max_n_hap1(clu_haps, clu_size, clu_read_ids, read_ids_n, n_clu, clu_poss, max_n_cons, n_het_pos);
        if (new_n_clu == n_clu) { // no further reassignment, but still have more than _max_n_cons_ clus
            err_func_printf(__func__, "%d small clusters of sequences remain un-assigned.", n_clu-max_n_cons);
            break;
        }
        n_clu = new_n_clu;
    }
    free(tup); free(clu_poss);
    return n_clu;
}

int reassign_hap(int **clu_haps, int *clu_size, uint64_t **clu_read_ids, int read_ids_n, int n_clu, int min_w, int max_n_cons, int n_het_pos) {
    // assign haplotype with reads < min_w to haplotype with reads >= min_w
    int new_n_clu = reassign_hap_by_min_w(clu_haps, clu_size, clu_read_ids, read_ids_n, n_clu, min_w, n_het_pos);
    if (new_n_clu > max_n_cons) // keep at most _max_n_cons_
        new_n_clu = reassign_max_n_hap(clu_haps, clu_size, clu_read_ids, read_ids_n, n_clu, min_w, n_het_pos, max_n_cons);
    // move max_n_cons to the front
    int i, j, pos_i;
    for (i = pos_i = 0; i < n_clu; ++i) {
        if (clu_size[i] == 0) continue;
        if (i == pos_i) {
            pos_i++; continue;
        }
        // move i to pos_i
        for (j = 0; j < read_ids_n; ++j) {
            clu_read_ids[pos_i][j] = clu_read_ids[i][j];
            clu_size[pos_i] = clu_size[i];
        }
        pos_i++;
    }
    if (pos_i > max_n_cons) err_fatal_core(__func__, "Error: collected %d clusters.", pos_i);
    return pos_i;
}

// collect minimized set of het bases
int abpoa_set_het_row_column_ids_weight(abpoa_graph_t *abg, uint64_t ***read_ids, int *het_poss, int **rc_weight, int msa_l, int n_seq, int m, int min_w, int read_ids_n) {
    int i, j, k, n, rank;
    uint64_t b, one = 1, *whole_read_ids = (uint64_t*)_err_calloc(read_ids_n, sizeof(uint64_t));
    for (i = 0; i < n_seq; ++i) {
        j = i / 64; b = i & 0x3f;
        whole_read_ids[j] |= (one << b);
    }
    for (i = 0; i < msa_l; ++i) {
        for (j = 0; j < read_ids_n; ++j) {
            read_ids[i][m-1][j] = whole_read_ids[j];
        }
    } free(whole_read_ids);

    uint8_t *node_map = (uint8_t*)_err_calloc(abg->node_n, sizeof(uint8_t));
    int *n_branch = (int*)_err_calloc(msa_l, sizeof(int)), n_het_pos = 0;
    for (i = 2; i < abg->node_n; ++i) {
        if (abg->node[i].out_edge_n < 2) continue;

        for (j = 0; j < abg->node[i].out_edge_n; ++j) {
            int out_id = abg->node[i].out_id[j];
            if (node_map[out_id]) continue;
            else node_map[out_id] = 1;
            int sum_out_w = 0;
            for (k = 0; k < abg->node[out_id].out_edge_n; ++k)
                sum_out_w += abg->node[out_id].out_weight[k];
            if (sum_out_w < min_w || sum_out_w > n_seq-min_w) continue;
            rank = abpoa_graph_node_id_to_msa_rank(abg, out_id); 
            n_branch[rank-1] += 1;
            // assign seq
            for (n = 0; n < abg->node[out_id].out_edge_n; ++n) {
                for (k = 0; k < abg->node[out_id].read_ids_n; ++k) {
                    b = abg->node[out_id].read_ids[n][k];
                    rc_weight[rank-1][abg->node[out_id].base] += get_bit_cnt4(ab_bit_table16, b);
                    read_ids[rank-1][abg->node[out_id].base][k] |= b;
                    read_ids[rank-1][m-1][k] ^= b;
                }
            }
            rc_weight[rank-1][m-1] -= rc_weight[rank-1][abg->node[out_id].base];
        }
    }
    for (rank = 0; rank < msa_l; ++rank) {
        if (rc_weight[rank][m-1] >= min_w && rc_weight[rank][m-1] <= n_seq-min_w) n_branch[rank]++;
        if (n_branch[rank] > 1) {
            // filter out identical read_ids
            int iden = 0;
            for (i = n_het_pos-1; i >= 0; i--) {
                int het_pos = het_poss[i];
                // remove het bases that share the identical read groups
                iden = abpoa_check_iden_read_ids(rc_weight, read_ids, m, read_ids_n, rank, het_pos);
                if (iden == 1) break;
            }
            if (iden == 1) continue;

            het_poss[n_het_pos++] = rank;
#ifdef __DEBUG__
            fprintf(stderr, "%d\t", rank);
            for (j = 0; j < m; ++j) {
                fprintf(stderr, "%c: %d\t", "ACGT-"[j], rc_weight[rank][j]);
            } fprintf(stderr, "\n");
#endif
        }
    }
    free(n_branch); free(node_map);
    return n_het_pos;
}

// group read into clusters based on all het bases
// initial cluster size could be > max_n_cons
int abpoa_collect_clu_hap_read_ids(int *het_poss, int n_het_pos, uint64_t ***read_ids, int read_ids_n, int n_seq, int m, int min_w, int max_n_cons, uint64_t ***clu_read_ids, int *_m_clu) {
    if (n_het_pos == 0) return 1;
    int i, j, k, n_clu = 0, m_clu = 2;
    int **clu_haps = (int**)_err_malloc(2 * sizeof(int*));
    int *clu_size = (int*)_err_calloc(2, sizeof(int));
    *clu_read_ids = (uint64_t**)_err_malloc(2 * sizeof(uint64_t**));
    for (i = 0; i < 2; ++i) {
        clu_haps[i] = (int*)_err_calloc(n_het_pos, sizeof(int));
        (*clu_read_ids)[i] = (uint64_t*)_err_calloc(read_ids_n, sizeof(uint64_t));
    }
    
    for (i = 0; i < n_seq; ++i) { // collect haplotype for each sequence
        int read_id_i = i / 64; uint64_t read_id = 1ULL << (i & 0x3f);
        for (j = 0; j < n_het_pos; ++j) {
            int het_pos = het_poss[j];
            for (k = 0; k < m; ++k) {
                if (read_ids[het_pos][k][read_id_i] & read_id) {
                    clu_haps[n_clu][j] = k; break;
                }
            }
        }
        if (check_redundent_hap(clu_haps, clu_size, *clu_read_ids, n_clu, n_clu, n_het_pos, read_id_i, read_id) == 0) {
            if (++n_clu == m_clu) {
                m_clu <<= 1;
                clu_haps = (int**)_err_realloc(clu_haps, m_clu * sizeof(int*));
                clu_size = (int*)_err_realloc(clu_size, m_clu * sizeof(int));
                (*clu_read_ids) = (uint64_t**)_err_realloc(*clu_read_ids, m_clu * sizeof(uint64_t**));
                for (j = n_clu; j < m_clu; ++j) {
                    clu_haps[j] = (int*)_err_calloc(n_het_pos, sizeof(int));
                    clu_size[j] = 0;
                    (*clu_read_ids)[j] = (uint64_t*)_err_calloc(read_ids_n, sizeof(uint64_t)); // mem may lost
                }
            }
        }
    }
    if (n_clu < 2) err_fatal(__func__, "# haplotypes: %d\n", n_clu);
#ifdef __DEBUG__
    fprintf(stderr, "n_clu: %d\n", n_clu);
    for (i = 0; i < n_clu; ++i) {
        for (j = 0; j < n_het_pos; ++j) {
            fprintf(stderr, "%d\t", clu_haps[i][j]);
        }
        fprintf(stderr, "\tsize: %d\n", clu_size[i]);
    }
#endif

    // assign haplotype with reads < min_w to haplotype with reads >= min_w
    // keep at most _max_n_cons_ haps and read ids, weight need to >= min_w
    n_clu = reassign_hap(clu_haps, clu_size, *clu_read_ids, read_ids_n, n_clu, min_w, max_n_cons, n_het_pos);
#ifdef __DEBUG__
    fprintf(stderr, "After re-assign: n_clu: %d\n", n_clu);
    for (i = 0; i < n_clu; ++i) {
        fprintf(stderr, "%d:\tsize: %d\n", i, clu_size[i]);
    }
#endif
    for (i = 0; i < m_clu; ++i) free(clu_haps[i]); free(clu_haps); free(clu_size);
    *_m_clu = m_clu;
    return n_clu;
}

// cluster reads into _n_clu_ groups based on heterogeneous bases
int abpoa_multip_read_clu(abpoa_graph_t *abg, int src_id, int sink_id, int n_seq, int m, int max_n_cons, double min_freq, uint64_t ***clu_read_ids, int *_m_clu) {
    abpoa_set_msa_rank(abg, src_id, sink_id);
    int i, j, n_clu, m_clu, read_ids_n = (n_seq-1)/64+1;
    int msa_l = abg->node_id_to_msa_rank[sink_id]-1, min_w = MAX_OF_TWO(1, n_seq * min_freq);
    
    // read_ids: support reads for each base (A/C/G/T) at each position
    uint64_t ***read_ids = (uint64_t***)_err_malloc(sizeof(uint64_t**) * msa_l);
    for (i = 0; i < msa_l; ++i) {
        read_ids[i] = (uint64_t**)_err_malloc(sizeof(uint64_t*) * m);
        for (j = 0; j < m; ++j) read_ids[i][j] = (uint64_t*)_err_calloc(read_ids_n, sizeof(uint64_t));
    }

    // is rc_weight necessary?
    int **rc_weight = (int**)_err_malloc(msa_l * sizeof(int*));
    for (i = 0; i < msa_l; ++i) {
        rc_weight[i] = (int*)_err_calloc(m, sizeof(int)); // ACGT
        rc_weight[i][m-1] = n_seq;
    } 
    // find min set of het nodes
    int *het_poss = (int*)_err_calloc(msa_l, sizeof(int));
    int n_het_pos = abpoa_set_het_row_column_ids_weight(abg, read_ids, het_poss, rc_weight, msa_l, n_seq, m, min_w, read_ids_n);
    
    if (n_het_pos < 1) n_clu = 1;
    // collect at most _max_n_cons_ haplotypes and corresponding read ids
    else n_clu = abpoa_collect_clu_hap_read_ids(het_poss, n_het_pos, read_ids, read_ids_n, n_seq, m, min_w, max_n_cons, clu_read_ids, &m_clu);

    for (i = 0; i < msa_l; ++i) {
        for (j = 0; j < m; ++j) free(read_ids[i][j]);
        free(read_ids[i]); free(rc_weight[i]);
    } free(read_ids); free(rc_weight); free(het_poss);

    *_m_clu = m_clu;
    return n_clu;
}

// should always do topological sort first, then generate consensus
void abpoa_generate_consensus(abpoa_t *ab, abpoa_para_t *abpt) {
    if (ab->abg->is_called_cons == 1) return;
    abpoa_graph_t *abg = ab->abg;
    if (abg->node_n <= 2) return;
    int i, *out_degree = (int*)_err_malloc(abg->node_n * sizeof(int));
    for (i = 0; i < abg->node_n; ++i) {
        out_degree[i] = abg->node[i].out_edge_n;
    }

    int n_clu, m_clu, n_seq = ab->abs->n_seq; uint64_t **clu_read_ids;
    int read_ids_n = (n_seq-1)/64+1;

    if (abpt->max_n_cons > 1) n_clu = abpoa_multip_read_clu(abg, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, n_seq, abpt->m, abpt->max_n_cons, abpt->min_freq, &clu_read_ids, &m_clu);
    else n_clu = 1;

    abpoa_cons_t *abc = ab->abc;
    abpoa_allocate_cons(abc, abg->node_n, ab->abs->n_seq, n_clu);
    if (n_clu > 1) {
         abpoa_multip_heaviest_bundling(abg, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, out_degree, n_clu, read_ids_n, clu_read_ids, abc);
        for (i = 0; i < m_clu; ++i) free(clu_read_ids[i]); free(clu_read_ids);
    } else {
        abpoa_heaviest_bundling(abg, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, out_degree, abc);
    }
    abg->is_called_cons = 1; free(out_degree);
}