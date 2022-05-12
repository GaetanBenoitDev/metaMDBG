#ifndef ABPOA_H
#define ABPOA_H

#include <stdint.h>
#include "simd_instruction.h"

#define ABPOA_GLOBAL_MODE 0
#define ABPOA_LOCAL_MODE  1
#define ABPOA_EXTEND_MODE 2
//#define ABPOA_SEMI_MODE 3

// gap mode
#define ABPOA_LINEAR_GAP 0
#define ABPOA_AFFINE_GAP 1
#define ABPOA_CONVEX_GAP 2

#define ABPOA_EXTRA_B 10
#define ABPOA_EXTRA_F 0.01

#define ABPOA_CIGAR_STR "MIDXSH"
#define ABPOA_CMATCH     0
#define ABPOA_CINS       1
#define ABPOA_CDEL       2
#define ABPOA_CDIFF      3
#define ABPOA_CSOFT_CLIP 4
#define ABPOA_CHARD_CLIP 5

#define ABPOA_SRC_NODE_ID  0
#define ABPOA_SINK_NODE_ID 1

#define ABPOA_OUT_CONS     0
#define ABPOA_OUT_MSA      1
#define ABPOA_OUT_CONS_MSA 2
#define ABPOA_OUT_GFA      3
#define ABPOA_OUT_CONS_GFA 4
#define ABPOA_OUT_CONS_FQ  5

#define ABPOA_HB 0
#define ABPOA_HC 1

// NOTE: upper boundary of in_edge_n is pow(2,30)
// for MATCH/MISMATCH: node_id << 34  | query_id << 4 | op
// for INSERTION:      query_id << 34 | op_len << 4   | op
// for DELETION:       node_id << 34  | op_len << 4   | op // op_len is always equal to 1
// for CLIP            query_id << 34 | op_len << 4   | op 
#define abpoa_cigar_t uint64_t 

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
    int n_cigar, m_cigar; abpoa_cigar_t *graph_cigar;
    int node_s, node_e, query_s, query_e; // for local and  extension mode
    int n_aln_bases, n_matched_bases;
    int32_t best_score; 
    // uint8_t is_rc:1; // is_rc: best_score is from the reverse complement
                        // now is_rc is determined based on minimizer-based seeding and chaining
} abpoa_res_t;

typedef struct {
    int m; int *mat; char *mat_fn; // score matrix
    int use_score_matrix; // set _mat_ based on score matrix file, then _match_/_mismatch_ is not used.
    int match, max_mat, mismatch, min_mis, gap_open1, gap_open2, gap_ext1, gap_ext2; int inf_min;
    // minimizer seeding parameter
    int k, w, min_w;
    int wb; float wf; // extra band width
    int zdrop, end_bonus; // from minimap2
    // int simd_flag; // available SIMD instruction
    // alignment mode
    uint8_t ret_cigar:1, rev_cigar:1, out_msa:1, out_cons:1, out_gfa:1, out_fq:1, use_read_ids:1, amb_strand:1;
    uint8_t disable_seeding:1, progressive_poa:1;
    char *incr_fn, *out_pog;
    int align_mode, gap_mode, max_n_cons;
    double min_freq; // for multiploid data
    int verbose; // to control output msg

    // char LogTable65536[65536];
    // char bit_table16[65536];
} abpoa_para_t;

typedef struct {
    int node_id;
    int in_edge_n, in_edge_m, *in_id;
    int out_edge_n, out_edge_m, *out_id; int *out_weight;
    uint64_t **read_ids; int read_ids_n; // for each edge

    int aligned_node_n, aligned_node_m, *aligned_node_id; // mismatch; aligned node will have same rank
    // int heaviest_weight, heaviest_out_id; // for consensus
    uint8_t base; // 0~m
    // ID, pos ???
} abpoa_node_t;

typedef struct {
    abpoa_node_t *node; int node_n, node_m, index_rank_m; 
    int *index_to_node_id;
    int *node_id_to_index, *node_id_to_max_pos_left, *node_id_to_max_pos_right, *node_id_to_max_remain, *node_id_to_msa_rank;
    uint8_t is_topological_sorted:1, is_called_cons:1, is_set_msa_rank:1;
} abpoa_graph_t;

typedef struct {
    int n_cons, n_seq, msa_len; // # cons, # of total seq, length of row-column MSA (including gaps)
    int *clu_n_seq;      // # of reads in each read cluster/group, size: n_cons
    int **clu_read_ids; // read ids for each cluster/group, size: n_cons * clu_n_seq[i]
    int *cons_len;       // length of each consensus sequence, size: n_cons
    int **cons_node_ids; // node id of each consensus, size: n_cons * cons_len[i]
    uint8_t **cons_base; // sequence base of each consensus, size: n_cons * cons_len[i]
    uint8_t **msa_base;  // sequence base of RC-MSA, size: (n_seq + n_cons) * msa_len
    int **cons_cov;      // coverage of each consensus base, size: n_cons * cons_len[i]
    int **cons_phred_score; // phred score for each consensus base, size: n_cons * cons_len[i]
} abpoa_cons_t;

typedef struct {
    int l, m; char *s;
} abpoa_str_t;

typedef struct {
    int n_seq, m_seq;
    abpoa_str_t *seq, *name, *comment, *qual;
    uint8_t *is_rc;
} abpoa_seq_t;

typedef struct {
    SIMDi *s_mem; uint64_t s_msize; // qp, DP_HE, dp_f OR qp, DP_H, dp_f : based on (qlen, num_of_value, m, node_n)
    int *dp_beg, *dp_end, *dp_beg_sn, *dp_end_sn, rang_m; // if band : based on (node_m)
} abpoa_simd_matrix_t;

typedef struct {
    abpoa_graph_t *abg;
    abpoa_seq_t *abs;
    abpoa_simd_matrix_t *abm;
    abpoa_cons_t *abc;
} abpoa_t;

// init for abpoa parameters
abpoa_para_t *abpoa_init_para(void);
void abpoa_set_mat_from_file(abpoa_para_t *abpt, char *mat_fn);
void abpoa_post_set_para(abpoa_para_t *abpt);
void abpoa_free_para(abpoa_para_t *abpt);

// init for alignment
abpoa_t *abpoa_init(void);
void abpoa_free(abpoa_t *ab);

// perform msa
int abpoa_msa(abpoa_t *ab, abpoa_para_t *abpt, int n_seqs, char **seq_names, int *seq_lens, uint8_t **seqs, FILE *out_fp);

int abpoa_msa1(abpoa_t *ab, abpoa_para_t *abpt, char *read_fn, FILE *out_fp);

// clean alignment graph
void abpoa_reset(abpoa_t *ab, abpoa_para_t *abpt, int qlen);

// restore graph from GFA/FASTA file
abpoa_t *abpoa_restore_graph(abpoa_t *ab, abpoa_para_t *abpt);

// for development:
// align a sequence to a graph
int abpoa_align_sequence_to_graph(abpoa_t *ab, abpoa_para_t *abpt, uint8_t *query, int qlen, abpoa_res_t *res);
// align a sequence to a graph between beg_node_id and end_node_id (both are excluded)
void abpoa_subgraph_nodes(abpoa_t *ab, abpoa_para_t *abpt, int inc_beg, int inc_end, int *exc_beg, int *exc_end);
int abpoa_align_sequence_to_subgraph(abpoa_t *ab, abpoa_para_t *abpt, int beg_node_id, int end_node_id, uint8_t *query, int qlen, abpoa_res_t *res);

// add a node to a graph
// para:
//   base: 0123 for ACGT
int abpoa_add_graph_node(abpoa_graph_t *abg, uint8_t base);

// add an edge to a graph
// para:
//   from_id/to_id: ids of from and to nodes
//   check_edge: set as 1 if this edge maybe alread exist and only need to update weight, set as 0 if the edge is new
//   add_read_id: set as 1 if read_id is used (to use row-column algorithm/generate MSA result/multiple consensus)
//   read_id: is of sequence
//   read_ids_n: size of read_id array, each one is 64-bit (1+(tot_read_n-1)/64)
int abpoa_add_graph_edge(abpoa_graph_t *abg, int from_id, int to_id, int check_edge, int w, uint8_t add_read_id, int read_id, int read_ids_n);

// add an alignment to a graph
// para:
//   query: 0123 for ACGT
//   qlen: query length
//   n_cigar/abpoa_cigar: from alignment result (abpoa_res_t)
//   read_id: id of sequence
//   tot_read_n: total number of sequence
int abpoa_add_graph_alignment(abpoa_t *ab, abpoa_para_t *abpt, uint8_t *query, int qlen, int *qpos_to_node_id, abpoa_res_t res, int read_id, int tot_read_n, int inc_both_ends);
int abpoa_add_subgraph_alignment(abpoa_t *ab, abpoa_para_t *abpt, int beg_node_id, int end_node_id, uint8_t *query, int qlen, int *qpos_to_node_id, abpoa_res_t res, int read_id, int tot_read_n, int inc_both_ends);

void abpoa_BFS_set_node_index(abpoa_graph_t *abg, int src_id, int sink_id);
void abpoa_BFS_set_node_remain(abpoa_graph_t *abg, int src_id, int sink_id);

// topological sortting of graph
void abpoa_topological_sort(abpoa_graph_t *abg, abpoa_para_t *abpt);

// generate consensus sequence from graph
// para:
//   out_fp: consensus sequence output in FASTA format, set as NULL to disable
//   cons_seq, cons_l, cons_n: store consensus sequences in variables, set cons_n as NULL to disable. 
//     cons_seq: store consensus sequences
//     cons_l: store consensus sequences length
//     cons_n: store number of consensus sequences
//     Note: cons_seq and cons_l need to be freed by user.
void abpoa_generate_consensus(abpoa_t *ab, abpoa_para_t *abpt);
void abpoa_output_fx_consensus(abpoa_t *ab, abpoa_para_t *abpt, FILE *out_fp);

// generate column multiple sequence alignment from graph
void abpoa_generate_rc_msa(abpoa_t *ab, abpoa_para_t *abpt);
void abpoa_output_rc_msa(abpoa_t *ab, abpoa_para_t *abpt, FILE *out_fp);

// generate graph in GFA format to _out_fp_
void abpoa_generate_gfa(abpoa_t *ab, abpoa_para_t *abpt, FILE *out_fp);

// output cons/msa
void abpoa_output(abpoa_t *ab, abpoa_para_t *abpt, FILE *out_fp);

// generate DOT graph plot and dump graph into PDF/PNG format file
void abpoa_dump_pog(abpoa_t *ab, abpoa_para_t *abpt);

#ifdef __cplusplus
}
#endif

#endif
