#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "abpoa.h"
#include "simd_abpoa_align.h"
#include "abpoa_align.h"
#include "abpoa_seq.h"
#include "abpoa_output.h"
#include "utils.h"
#include "abpoa_seed.h"

void gen_simple_mat(abpoa_para_t *abpt) {
    int m = abpt->m, i, j;
    int match = abpt->match < 0 ? -abpt->match : abpt->match;
    int mismatch = abpt->mismatch > 0? -abpt->mismatch : abpt->mismatch;
    for (i = 0; i < m - 1; ++i) {
        for (j = 0; j < m - 1; ++j)
            abpt->mat[i * m + j] = i == j ? match : mismatch;
        abpt->mat[i * m + m - 1] = 0;
    }
    for (j = 0; j < m; ++j)
        abpt->mat[(m - 1) * m + j] = 0;
    abpt->max_mat = match;
    abpt->min_mis = -mismatch;
}

extern char ab_nt4_table[256];
extern char ab_nt256_table[256];
extern char ab_aa26_table[256];
extern char ab_aa256_table[256];
extern char ab_char26_table[256];
extern char ab_char256_table[256];

void parse_mat_first_line(char *l, int *order) {
    int i, n;
    for (i = n = 0; l[i]; ++i) {
        if (isspace(l[i])) continue;
        order[n++] = ab_char26_table[(int)l[i]];
    }
}

void parse_mat_score_line(char *l, int *order, int m, int *mat) {
    int n, is_base=1, _i=-1; long s; char *str = l, *pEnd=NULL;
    for (n = 0; *str; ++str) {
        if (!isalpha(*str) && !isdigit(*str) && *str != '+' && *str != '-') continue;
        if (is_base) { // get base
            _i = ab_char26_table[(int)*str];
            if (_i >= m) err_fatal(__func__, "Unknown base: \"%c\" (%d).\n", *str, _i);
            is_base = 0;
        } else { // get score
            if (n == m) 
                err_fatal_simple("Too many scores in matrix.\n");
            s = strtol(str, &pEnd, 10);
            str = pEnd;
            mat[_i *m + order[n]] = s;
            n++;
        }
    }
}

void abpoa_set_mat_from_file(abpoa_para_t *abpt, char *mat_fn) {
    char *l = (char*)_err_malloc(1024 * sizeof(char)); FILE *fp;
    if ((fp = fopen(mat_fn, "r")) == NULL) err_fatal(__func__, "Unable to open scoring matrix file: \"%s\"\n", mat_fn);
    int first_line = 1;
    int *order = (int*)_err_malloc(abpt->m * sizeof(int));
    while (fgets(l, 1024, fp) != NULL) {
        if (l[0] == '#') continue;
        if (first_line) {
            first_line = 0;
            // get A/C/G/T/N bases
            parse_mat_first_line(l, order);
        } else {
            // get match/mismatch scores
            parse_mat_score_line(l, order, abpt->m, abpt->mat);
        }
    }
    int i; abpt->min_mis = 0, abpt->max_mat = 0;
    for (i = 0; i < abpt->m * abpt->m; ++i) {
        if (abpt->mat[i] > abpt->max_mat)
            abpt->max_mat = abpt->mat[i];
        if (-abpt->mat[i] > abpt->min_mis) 
            abpt->min_mis = -abpt->mat[i];
    }
    free(l); free(order); fclose(fp);
}

void abpoa_set_gap_mode(abpoa_para_t *abpt) {
    if (abpt->gap_open1 == 0) abpt->gap_mode = ABPOA_LINEAR_GAP;
    else if (abpt->gap_open1 > 0 && abpt->gap_open2 == 0) abpt->gap_mode = ABPOA_AFFINE_GAP;
    else abpt->gap_mode = ABPOA_CONVEX_GAP;
}

abpoa_para_t *abpoa_init_para(void) {
    abpoa_para_t *abpt = (abpoa_para_t*)_err_malloc(sizeof(abpoa_para_t));
    abpt->align_mode = ABPOA_GLOBAL_MODE;
    abpt->gap_mode = ABPOA_CONVEX_GAP;
    abpt->zdrop = -1;     // disable zdrop
    abpt->end_bonus = -1; // disable end bouns
    abpt->wb = ABPOA_EXTRA_B; // extra bandwidth
    abpt->wf = ABPOA_EXTRA_F; // extra bandwidth

    abpt->amb_strand = 0; // ambiguous strand
    abpt->ret_cigar = 1;  // return cigar
    abpt->rev_cigar = 0;  // reverse cigar
    abpt->out_cons = 1;   // output consensus sequence in fasta
    abpt->out_fq = 0;     // output consensus sequence in fastq
    abpt->out_gfa = 0;    // out graph in GFA format
    abpt->out_msa = 0;    // output msa
    abpt->max_n_cons = 1; // number of max. generated consensus sequence
    abpt->min_freq = MULTIP_MIN_FREQ; 
    abpt->use_read_ids = 0;
    abpt->incr_fn = NULL; // incrementally align seq to an existing graph
    abpt->out_pog = NULL; // dump partial order graph to file

    // number of residue types
    abpt->m = 5; // nucleotide
    abpt->mat = (int*)_err_malloc(abpt->m * abpt->m * sizeof(int));

    // score matrix
    abpt->use_score_matrix = 0;
    abpt->mat_fn = NULL;
    abpt->match = ABPOA_MATCH;
    abpt->mismatch = ABPOA_MISMATCH;
    abpt->gap_open1 = ABPOA_GAP_OPEN1;
    abpt->gap_open2 = ABPOA_GAP_OPEN2;
    abpt->gap_ext1 = ABPOA_GAP_EXT1;
    abpt->gap_ext2 = ABPOA_GAP_EXT2;

    abpt->disable_seeding = 1; // no seeding by default
    abpt->k = ABPOA_MMK;
    abpt->w = ABPOA_MMW;
    abpt->min_w = ABPOA_MIN_POA_WIN;
    abpt->progressive_poa = 0; // progressive partial order alignment

    abpt->verbose = 0;

    // abpt->simd_flag = simd_check();

    return abpt;
}

void abpoa_post_set_para(abpoa_para_t *abpt) {
    abpoa_set_gap_mode(abpt);
    if (abpt->out_msa || abpt->out_gfa || abpt->max_n_cons > 1) {
        abpt->use_read_ids = 1;
        set_65536_table();
        if (abpt->max_n_cons > 1) set_bit_table16();
    }
    if (abpt->align_mode == ABPOA_LOCAL_MODE) abpt->wb = -1;
    int i;
    if (abpt->m > 5) { // for aa sequence
        for (i = 0; i < 256; ++i) {
            ab_char26_table[i] = ab_aa26_table[i];
            ab_char256_table[i] = ab_aa256_table[i];
        }
        if (abpt->k > 11) {
            abpt->k = 7, abpt->w = 4;
        }
    } else {
        for (i = 0; i < 256; ++i) {
            ab_char26_table[i] = ab_nt4_table[i];
            ab_char256_table[i] = ab_nt256_table[i];
        }
    }
    if (abpt->use_score_matrix == 0) gen_simple_mat(abpt);
    else abpoa_set_mat_from_file(abpt, abpt->mat_fn);
}

void abpoa_free_para(abpoa_para_t *abpt) {
    if (abpt->mat != NULL) free(abpt->mat);
    if (abpt->mat_fn != NULL) free(abpt->mat_fn);
    if (abpt->out_pog != NULL) free(abpt->out_pog);
    if (abpt->incr_fn != NULL) free(abpt->incr_fn);
    free(abpt);
}

int abpoa_align_sequence_to_subgraph(abpoa_t *ab, abpoa_para_t *abpt, int exc_beg_node_id, int exc_end_node_id, uint8_t *query, int qlen, abpoa_res_t *res) {
    if (ab->abg->node_n <= 2) return -1;
    if (ab->abg->is_topological_sorted == 0) abpoa_topological_sort(ab->abg, abpt);
    simd_abpoa_align_sequence_to_subgraph(ab, abpt, exc_beg_node_id, exc_end_node_id, query, qlen, res);
    return 0;
}

int abpoa_align_sequence_to_graph(abpoa_t *ab, abpoa_para_t *abpt, uint8_t *query, int qlen, abpoa_res_t *res) {
    if (ab->abg->node_n <= 2) return -1;
    if (ab->abg->is_topological_sorted == 0) abpoa_topological_sort(ab->abg, abpt);
    simd_abpoa_align_sequence_to_graph(ab, abpt, query, qlen, res);
    return 0;
}

int abpoa_anchor_poa(abpoa_t *ab, abpoa_para_t *abpt, uint8_t **seqs, int *seq_lens, ab_u64_v par_anchors, int *par_c, int *tpos_to_node_id, int *qpos_to_node_id, int *read_id_map, int exist_n_seq, int n_seq) {
    // err_func_format_printf(__func__, "Performing POA between anchors ...");
    abpoa_res_t res; int read_id, last_read_id = -1, m_c = 0, k = abpt->k, qlen;
    abpoa_seq_t *abs = ab->abs;
    int *tmp;
    int i, _i, ai, j, tot_n_seq = exist_n_seq + n_seq;
    uint8_t *qseq; abpoa_res_t whole_res;
    // uint8_t *seq1;
    for (_i = 0; _i < n_seq; ++_i) {
        i = read_id_map[_i]; read_id = exist_n_seq + i; qlen = seq_lens[i]; whole_res.n_cigar = 0, whole_res.m_cigar = 0, whole_res.graph_cigar = 0;
#ifdef __DEBUG__
        fprintf(stderr, "seq: # %d\n", i);
#endif
        // seq-to-graph alignment and add alignment within each split window
        if (_i == 0) ai = 0; else ai = par_c[_i-1];

        int beg_id = ABPOA_SRC_NODE_ID, beg_qpos = 0, end_id=-1, end_tpos=-1, end_qpos=-1;
        if (ai < par_c[_i]) {
            abs->is_rc[read_id] = (abs->is_rc[last_read_id] ^ (par_anchors.a[ai] >> 63));
            // construct rc qseq
            if (abs->is_rc[read_id]) {
                qseq = (uint8_t*)_err_malloc(qlen * sizeof(uint8_t));
                for (j = 0; j < qlen; ++j) {
                    if (seqs[i][qlen-j-1] < 4) qseq[j] = 3 - seqs[i][qlen-j-1];
                    else qseq[j] = 4;
                }
                if (abs->is_rc[last_read_id]) { // reset tpos/qpos in par_anchors
                    int last_qlen = seq_lens[read_id_map[_i-1]];
                    for (j = ai; j < par_c[_i]; ++j) {
                        end_tpos = ((par_anchors.a[j] >> 32) & 0x7fffffff); end_qpos = (uint32_t)par_anchors.a[j];
                        par_anchors.a[j] = (par_anchors.a[j] >> 63) << 63 | (uint64_t)(last_qlen-end_tpos+k) << 32 | (qlen-end_qpos+k);
                    }
                    for (j = 0; j < (par_c[_i]-ai)/2; ++j) {
                        uint64_t tmp = par_anchors.a[ai+j];
                        par_anchors.a[ai+j] = par_anchors.a[par_c[_i]-1-j];
                        par_anchors.a[par_c[_i]-1-j] = tmp;
                    }
                }
            } else { 
                qseq = seqs[i];
                if (abs->is_rc[last_read_id]) { // reset tpos/qpos in par_anchors 
                    int last_qlen = seq_lens[read_id_map[_i-1]];
                    for (j = ai; j < par_c[_i]; ++j) {
                        end_tpos = ((par_anchors.a[j] >> 32) & 0x7fffffff); end_qpos = (uint32_t)par_anchors.a[j];
                        par_anchors.a[j] = (par_anchors.a[j] >> 63) << 63 | (uint64_t)(last_qlen-end_tpos+k) << 32 | (qlen-end_qpos+k);
                    }
                    for (j = 0; j < (par_c[_i]-ai)/2; ++j) {
                        uint64_t tmp = par_anchors.a[ai+j];
                        par_anchors.a[ai+j] = par_anchors.a[par_c[_i]-1-j];
                        par_anchors.a[par_c[_i]-1-j] = tmp;
                    }
                }
            }
        } else abs->is_rc[read_id] = 0, qseq = seqs[i];

        for (; ai < par_c[_i]; ++ai) {
            end_tpos = ((par_anchors.a[ai] >> 32) & 0x7fffffff) - k + 1; end_id = tpos_to_node_id[end_tpos];
            end_qpos = (uint32_t)par_anchors.a[ai] - k + 1;

#ifdef __DEBUG__
            fprintf(stderr, "\tanchor: t: %d (id: %d), q: %d\n", end_tpos, end_id, end_qpos);
#endif

            res.graph_cigar = 0; res.n_cigar = 0;
            abpoa_align_sequence_to_subgraph(ab, abpt, beg_id, end_id, qseq+beg_qpos, end_qpos-beg_qpos, &res);
            abpoa_push_whole_cigar(&whole_res.n_cigar, &whole_res.m_cigar, &whole_res.graph_cigar, res.n_cigar, res.graph_cigar);
            if (res.n_cigar) free(res.graph_cigar);
            // abpoa_add_subgraph_alignment(ab, abpt, beg_id, end_id, qseq+beg_qpos, end_qpos-beg_qpos, qpos_to_node_id+beg_qpos, res, read_id, tot_n_seq, 1);

            // add alignment for anchors
            res.graph_cigar = (abpoa_cigar_t*)_err_malloc((k) * sizeof(abpoa_cigar_t)); res.n_cigar = 0; m_c = k;
            for (j = 0; j < k; ++j)
                res.graph_cigar = abpoa_push_cigar(&(res.n_cigar), &m_c, res.graph_cigar, ABPOA_CMATCH, 1, tpos_to_node_id[end_tpos+j], j);
            // for (j = 0; j < k; ++j) qpos_to_node_id[end_qpos+j] = tpos_to_node_id[end_tpos+j];
            // abpoa_add_subgraph_alignment(ab, abpt, end_id, tpos_to_node_id[end_tpos+k-1], qseq+end_qpos, k, NULL, res, read_id, tot_n_seq, 1);
            abpoa_push_whole_cigar(&whole_res.n_cigar, &whole_res.m_cigar, &whole_res.graph_cigar, res.n_cigar, res.graph_cigar);
            if (res.n_cigar) free(res.graph_cigar);

            // for next anchor
            beg_id = tpos_to_node_id[end_tpos+k-1]; beg_qpos = end_qpos+k;
        }
        end_id = ABPOA_SINK_NODE_ID; end_qpos = seq_lens[i];

#ifdef __DEBUG__
            fprintf(stderr, "\tanchor: t: %d (id: %d), q: %d\n", end_tpos, end_id, end_qpos);
#endif
        res.graph_cigar = 0; res.n_cigar = 0;
        abpoa_align_sequence_to_subgraph(ab, abpt, beg_id, end_id, qseq+beg_qpos, end_qpos-beg_qpos, &res);
        abpoa_push_whole_cigar(&whole_res.n_cigar, &whole_res.m_cigar, &whole_res.graph_cigar, res.n_cigar, res.graph_cigar);
        if (res.n_cigar) free(res.graph_cigar);

        abpoa_add_subgraph_alignment(ab, abpt, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, qseq, qlen, qpos_to_node_id, whole_res, read_id, tot_n_seq, 1);
        if (abs->is_rc[read_id]) free(qseq);

        if (whole_res.n_cigar) free(whole_res.graph_cigar);

        tmp = qpos_to_node_id; qpos_to_node_id = tpos_to_node_id; tpos_to_node_id = tmp;
        last_read_id = read_id;
    }
    // err_func_format_printf(__func__, "Performing POA between anchors done.");
    return 0;
}

// simply partial order alignment, no seeding-based anchor or progressive tree
int abpoa_poa(abpoa_t *ab, abpoa_para_t *abpt, uint8_t **seqs, int *seq_lens, int exist_n_seq, int n_seq) {
    // err_func_format_printf(__func__, "Performing POA ...");
    abpoa_seq_t *abs = ab->abs;
    abpoa_res_t res; int i, j, read_id, qlen, tot_n_seq = exist_n_seq + n_seq;
    uint8_t *qseq, *rc_qseq;
    // uint8_t *seq1;
    for (i = 0; i < n_seq; ++i) {
        qlen = seq_lens[i]; qseq = seqs[i]; read_id = exist_n_seq + i;
#ifdef __DEBUG__
        fprintf(stderr, "seq: # %d\n", i);
#endif
        res.graph_cigar = 0; res.n_cigar = 0;
        if (abpoa_align_sequence_to_graph(ab, abpt, qseq, qlen, &res) >= 0) {
            if (abpt->amb_strand && (res.best_score < MIN_OF_TWO(qlen, ab->abg->node_n-2) * abpt->max_mat * .3333)) { // TODO .3333
                rc_qseq = (uint8_t*)_err_malloc(sizeof(uint8_t) * qlen);
                for (j = 0; j < qlen; ++j) {
                    if (qseq[qlen-j-1] < 4) rc_qseq[j] = 3 - qseq[qlen-j-1];
                    else rc_qseq[j] = 4;
                }
                abpoa_res_t rc_res; rc_res.n_cigar = 0, rc_res.graph_cigar = 0;
                simd_abpoa_align_sequence_to_graph(ab, abpt, rc_qseq, qlen, &rc_res);
                if (rc_res.best_score > res.best_score) {
                    abpoa_res_copy(&res, &rc_res);
                    qseq = rc_qseq;
                    abs->is_rc[read_id] = 1;
                } else free(rc_qseq);
                if (rc_res.n_cigar) free(rc_res.graph_cigar);
            } 
        }
        abpoa_add_graph_alignment(ab, abpt, qseq, qlen, NULL, res, read_id, tot_n_seq, 1);
        if (abs->is_rc[read_id]) free(qseq);
        if (res.n_cigar) free(res.graph_cigar);
    }
    // err_func_format_printf(__func__, "Performing POA ... done.");
    return 0;
}

void abpoa_output(abpoa_t *ab, abpoa_para_t *abpt, FILE *out_fp) {
    // generate & output GFA
    if (abpt->out_gfa) abpoa_generate_gfa(ab, abpt, out_fp);
    else {
        // generate rc-msa/cons
        if (abpt->out_msa) abpoa_generate_rc_msa(ab, abpt);
        if (abpt->out_cons) {
            abpoa_generate_consensus(ab, abpt);
            if (ab->abg->is_called_cons == 0) err_printf("Warning: no consensus sequence generated.\n");
        }
        // output cons/rc-msa
        if (abpt->out_msa) abpoa_output_rc_msa(ab, abpt, out_fp);
        else if (abpt->out_cons) abpoa_output_fx_consensus(ab, abpt, out_fp);
    }
    // plot partial-order graph using dot
    if (abpt->out_pog) abpoa_dump_pog(ab, abpt);
}

// do msa for a set of input sequences
// @function: 
//    generate consensus sequence
//    generate rc-msa (row column multiple sequence alignment)
// @para:
//    ab/abpt: abpoa related variable and parameter 
//    n_seq: number of input sequences
//    seq_len: array of input sequence length, size: seq_n
//    seqs: array of input sequences, 0123 for ACGT, size: seq_n * seq_len[]
int abpoa_msa(abpoa_t *ab, abpoa_para_t *abpt, int n_seq, char **seq_names, int *seq_lens, uint8_t **seqs, FILE *out_fp) {
    if ((!abpt->out_msa && !abpt->out_cons && !abpt->out_gfa) || n_seq <= 0) return 0;
    abpoa_reset(ab, abpt, 1024);
    if (abpt->incr_fn) abpoa_restore_graph(ab, abpt); // restore existing graph
    abpoa_seq_t *abs = ab->abs; int i, exist_n_seq = abs->n_seq;

    // set ab->abs, name
    abs->n_seq += n_seq; abpoa_realloc_seq(abs);

    if (seq_names) {
        for (i = 0; i < n_seq; ++i) {
            abpoa_cpy_str(abs->name+exist_n_seq+i, seq_names[i], strlen(seq_names[i]));
        }
    } else {
        for (i = 0; i < n_seq; ++i) {
            abs->name[exist_n_seq+i].l = 0; abs->name[exist_n_seq+i].m = 0;
        }
    }

    // always reset graph before perform POA
    int max_len = 0;
    for (i = 0; i < n_seq; ++i) {
        if (seq_lens[i] > max_len) max_len = seq_lens[i];
    }


    if ((abpt->disable_seeding && abpt->progressive_poa==0) || abpt->align_mode != ABPOA_GLOBAL_MODE) {
        abpoa_poa(ab, abpt, seqs, seq_lens, exist_n_seq, n_seq);
    } else {
        // sequence pos to node id
        int *tpos_to_node_id = (int*)_err_calloc(max_len, sizeof(int)), *qpos_to_node_id = (int*)_err_calloc(max_len, sizeof(int));
        // seeding, build guide tree, and partition into small windows
        int *read_id_map = (int*)_err_malloc(sizeof(int) * n_seq); // guide tree order -> input order
        ab_u64_v par_anchors = {0, 0, 0}; int *par_c = (int*)_err_calloc(n_seq, sizeof(int));

        abpoa_build_guide_tree_partition(seqs, seq_lens, n_seq, abpt, read_id_map, &par_anchors, par_c);
        if (abpt->incr_fn) { // collect anchors between last one path and first seq
            // anchors
            // new_par_anchors
            // push anchors 
            // free(par_anchors.a);
            // par_anchors = new_par_anchors;
            // collect tpos_to_node_id for last one path
        }

        // perform partial order alignment
        abpoa_anchor_poa(ab, abpt, seqs, seq_lens, par_anchors, par_c, tpos_to_node_id, qpos_to_node_id, read_id_map, exist_n_seq, n_seq);
        free(read_id_map); free(tpos_to_node_id); free(qpos_to_node_id); free(par_c);
        if (par_anchors.m > 0) free(par_anchors.a);
    }

    // output
    abpoa_output(ab, abpt, out_fp);
    return 0;
}

int abpoa_msa1(abpoa_t *ab, abpoa_para_t *abpt, char *read_fn, FILE *out_fp) {
    if (!abpt->out_msa && !abpt->out_cons && !abpt->out_gfa) return 0;
    abpoa_reset(ab, abpt, 1024);
    if (abpt->incr_fn) abpoa_restore_graph(ab, abpt); // restore existing graph
    abpoa_seq_t *abs = ab->abs; int exist_n_seq = abs->n_seq;

    // read seq from read_fn
    gzFile readfp = xzopen(read_fn, "r"); kseq_t *ks = kseq_init(readfp);
    int i, j, n_seq = abpoa_read_seq(abs, ks);

    // always reset graph before perform POA
    int max_len = 0;
    for (i = 0; i < abs->n_seq; ++i) {
        if (abs->seq[i].l > max_len) max_len = abs->seq[i].l;
    }

    // set seqs, seq_lens
    extern char ab_char26_table[256];
    uint8_t **seqs = (uint8_t**)_err_malloc(n_seq * sizeof(uint8_t*)); int *seq_lens = (int*)_err_malloc(n_seq * sizeof(int));
    for (i = 0; i < n_seq; ++i) {
        seq_lens[i] = abs->seq[exist_n_seq+i].l;
        seqs[i] = (uint8_t*)_err_malloc(sizeof(uint8_t) * seq_lens[i]);
        for (j = 0; j < seq_lens[i]; ++j) seqs[i][j] = ab_char26_table[(int)abs->seq[exist_n_seq+i].s[j]];
    }
    if ((abpt->disable_seeding && abpt->progressive_poa==0) || abpt->align_mode != ABPOA_GLOBAL_MODE) {
        abpoa_poa(ab, abpt, seqs, seq_lens, exist_n_seq, n_seq);
    } else {
        // sequence pos to node id
        int *tpos_to_node_id = (int*)_err_calloc(max_len, sizeof(int)), *qpos_to_node_id = (int*)_err_calloc(max_len, sizeof(int));
        // seeding, build guide tree, and partition into small windows
        int *read_id_map = (int*)_err_malloc(sizeof(int) * n_seq); // guide tree order -> input order
        ab_u64_v par_anchors = {0, 0, 0}; int *par_c = (int*)_err_calloc(n_seq, sizeof(int));

        abpoa_build_guide_tree_partition(seqs, seq_lens, n_seq, abpt, read_id_map, &par_anchors, par_c);
        if (abpt->incr_fn) { // TODO collect anchors between last one path and first seq
            // anchors
            // new_par_anchors
            // push anchors 
            // free(par_anchors.a);
            // par_anchors = new_par_anchors;
            // collect tpos_to_node_id for last one path
            // set tpos_to_node_id
            //
        }
        abpoa_anchor_poa(ab, abpt, seqs, seq_lens, par_anchors, par_c, tpos_to_node_id, qpos_to_node_id, read_id_map, exist_n_seq, n_seq);
        free(read_id_map); free(tpos_to_node_id); free(qpos_to_node_id); free(par_c);
        if (par_anchors.m > 0) free(par_anchors.a);
    }

    // output
    abpoa_output(ab, abpt, out_fp);

    kseq_destroy(ks); gzclose(readfp);
    for (i = 0; i < n_seq; ++i) free(seqs[i]); free(seqs); free(seq_lens);
    return 0;
}
