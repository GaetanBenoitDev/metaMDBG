#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "abpoa_align.h"
#include "simd_instruction.h"
#include "utils.h"

typedef struct {
    const int reg_n, bits_n, log_num, num_of_value, size;
    int inf_min; // based on penalty of mismatch and GAP_OE1
} SIMD_para_t;

#define SIMDShiftOneNi8  1
#define SIMDShiftOneNi16 2
#define SIMDShiftOneNi32 4
#define SIMDShiftOneNi64 8

#ifdef __AVX512F__
SIMD_para_t _simd_p8  = {512,  8, 6, 64, 64, -1};
SIMD_para_t _simd_p16 = {512, 16, 5, 32, 64, -1};
SIMD_para_t _simd_p32 = {512, 32, 4, 16, 64, -1};
SIMD_para_t _simd_p64 = {512, 64, 3,  8, 64, -1};
#define SIMDTotalBytes 64
#elif defined(__AVX2__)
SIMD_para_t _simd_p8  = {256,  8, 5, 32, 32, -1};
SIMD_para_t _simd_p16 = {256, 16, 4, 16, 32, -1};
SIMD_para_t _simd_p32 = {256, 32, 3,  8, 32, -1};
SIMD_para_t _simd_p64 = {256, 64, 2,  4, 32, -1};
#define SIMDTotalBytes 32
#else
SIMD_para_t _simd_p8  = {128,  8, 4, 16, 16, -1};
SIMD_para_t _simd_p16 = {128, 16, 3,  8, 16, -1};
SIMD_para_t _simd_p32 = {128, 32, 2,  4, 16, -1};
SIMD_para_t _simd_p64 = {128, 64, 1,  2, 16, -1};
#define SIMDTotalBytes 16
#endif

#define print_simd(s, str, score_t) {                                \
    int _i; score_t *_a = (score_t*)(s);                             \
    fprintf(stderr, "%s\t", str);                                    \
    for (_i = 0; _i < SIMDTotalBytes / (int)sizeof(score_t); ++_i) { \
        fprintf(stderr, "%d\t", _a[_i]);                             \
    } fprintf(stderr, "\n");                                         \
}


#define simd_abpoa_print_lg_matrix(score_t, beg_index, end_index) { \
    for (j = 0; j < end_index-beg_index; ++j) {                     \
        fprintf(stderr, "index: %d\t", j);                          \
        dp_h = DP_H + j * dp_sn;                                    \
        _dp_h = (score_t*)dp_h;                                     \
        for (i = dp_beg[j]; i <= dp_end[j]; ++i) {                  \
            fprintf(stderr, "%d:(%d)\t", i, _dp_h[i]);              \
        } fprintf(stderr, "\n");                                    \
    }                                                               \
}

#define simd_abpoa_print_ag_matrix(score_t, beg_index, end_index) {  \
    for (j = beg_index; j < end_index; ++j) {                        \
        fprintf(stderr, "index: %d\t", j);                           \
        dp_h = DP_HEF + j * 3 * dp_sn; dp_e1 = dp_h + dp_sn;         \
        _dp_h = (score_t*)dp_h, _dp_e1 = (score_t*)dp_e1;            \
        for (i = dp_beg[j]; i <= dp_end[j]; ++i) {                   \
            fprintf(stderr, "%d:(%d,%d)\t", i, _dp_h[i], _dp_e1[i]); \
        } fprintf(stderr, "\n");                                     \
    }                                                                \
}

#define debug_simd_abpoa_print_cg_matrix_row(str, score_t, index_i) {                                    \
    score_t *_dp_h = (score_t*)dp_h, *_dp_e1 = (score_t*)dp_e1;                                          \
    score_t *_dp_e2 = (score_t*)dp_e2, *_dp_f1 = (score_t*)dp_f1, *_dp_f2 = (score_t*)dp_f2;             \
    fprintf(stderr, "%s\tindex: %d\t", str, index_i);                                                    \
    for (i = dp_beg[index_i]; i <= (dp_end[index_i]/16+1)*16-1; ++i) {                                   \
        fprintf(stderr, "%d:(%d,%d,%d,%d,%d)\t", i, _dp_h[i], _dp_e1[i],_dp_e2[i], _dp_f1[i],_dp_f2[i]); \
    } fprintf(stderr, "\n");                                                                             \
}

#define simd_abpoa_print_cg_matrix(score_t, beg_index, end_index) {                                          \
    for (j = 0; j < end_index-beg_index; ++j) {                                                              \
        fprintf(stderr, "index: %d\t", j);                                                                   \
        dp_h=DP_H2E2F+j*5*dp_sn; dp_e1=dp_h+dp_sn; dp_e2=dp_e1+dp_sn; dp_f1=dp_e2+dp_sn; dp_f2=dp_f1+dp_sn;  \
        score_t *_dp_h=(score_t*)dp_h, *_dp_e1=(score_t*)dp_e1, *_dp_e2=(score_t*)dp_e2;                     \
        score_t *_dp_f1=(score_t*)dp_f1, *_dp_f2=(score_t*)dp_f2;                                            \
        for (i = dp_beg[j]; i <= dp_end[j]; ++i) {                                                           \
            fprintf(stderr, "%d:(%d,%d,%d,%d,%d)\t", i, _dp_h[i], _dp_e1[i],_dp_e2[i], _dp_f1[i],_dp_f2[i]); \
        } fprintf(stderr, "\n");                                                                             \
    }                                                                                                        \
}

/* max_pos_left/right: left/right boundary of max column index for each row, based on the pre_nodes' DP score */
/* === workflow of alignment === */
/* a. global:
 * (1) alloc mem
 * (2) init for first row
 * (3) DP for each row
 * (3.2) if use_ada, update max_pos_left/right
 * (4) find best_i/j, backtrack
 * b. extend:
 * (1) alloc mem
 * (2) init for first row
 * (3) DP for each row
 * (3.2) find max of current row
 * (3.3) z-drop, set_max_score
 * (3.4) if use_ada, update max_pos_left/right
 */

// backtrack order:
// Match/Mismatch, Deletion, Insertion
#define simd_abpoa_lg_backtrack(score_t) {                                                                  \
    int i, j, k, pre_i, n_c = 0, s, is_match, m_c = 0, hit, id, _start_i, _start_j;                         \
    SIMDi *dp_h; score_t *_dp_h=NULL, *_pre_dp_h; abpoa_cigar_t *cigar = 0;                                 \
    i = best_i, j = best_j, _start_i = best_i, _start_j = best_j;                                           \
    id = abpoa_graph_index_to_node_id(graph, i+beg_index);                                                  \
    if (best_j < qlen) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, qlen-j, -1, qlen-1);         \
    dp_h = DP_H + i * dp_sn; _dp_h = (score_t*)dp_h;                                                        \
    int indel_first = 1; /* prefer to keep gaps at the end */                                               \
    while (i > 0 && j > 0) {                                                                                \
        if (abpt->align_mode == ABPOA_LOCAL_MODE && _dp_h[j] == 0) break;                                   \
        _start_i = i, _start_j = j;                                                                         \
        int *pre_index_i = pre_index[i];                                                                    \
        s = mat[m * graph->node[id].base + query[j-1]]; hit = 0;                                            \
        is_match = graph->node[id].base == query[j-1];                                                      \
        if (indel_first == 0) { /* match/mismatch */                                                        \
            for (k = 0; k < pre_n[i]; ++k) { /* match/mismatch */                                           \
                pre_i = pre_index_i[k];                                                                     \
                if (j-1 < dp_beg[pre_i] || j-1 > dp_end[pre_i]) continue;                                   \
                _pre_dp_h = (score_t*)(DP_H + pre_i * dp_sn);                                               \
                if (_pre_dp_h[j-1] + s == _dp_h[j]) {                                                       \
                    cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, j-1);                  \
                    i = pre_i; --j; hit = 1; id = abpoa_graph_index_to_node_id(graph, i+beg_index);         \
                    dp_h = DP_H + i * dp_sn; _dp_h = (score_t*)dp_h;                                        \
                    ++res->n_aln_bases; res->n_matched_bases += is_match ? 1 : 0;                           \
                    break;                                                                                  \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0) { /* deletion */                                                                      \
            for (k = 0; k < pre_n[i]; ++k) {                                                                \
                pre_i = pre_index_i[k];                                                                     \
                if (j < dp_beg[pre_i] || j > dp_end[pre_i]) continue;                                       \
                _pre_dp_h = (score_t*)( DP_H + pre_i * dp_sn);                                              \
                if (_pre_dp_h[j] - gap_ext1 == _dp_h[j]) {                                                  \
                    cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);                    \
                    i = pre_i; hit = 1; id = abpoa_graph_index_to_node_id(graph, i+beg_index);              \
                    dp_h = DP_H + i * dp_sn; _dp_h = (score_t*)dp_h;                                        \
                    break;                                                                                  \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0) { /* insertion */                                                                     \
            if (_dp_h[j-1] - gap_ext1 == _dp_h[j]) {                                                        \
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, 1, id, j-1); j--;                   \
                hit = 1; ++res->n_aln_bases;                                                                \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0 && indel_first == 1) { /* match/mismatch */                                            \
            for (k = 0; k < pre_n[i]; ++k) {                                                                \
                pre_i = pre_index_i[k];                                                                     \
                if (j-1 < dp_beg[pre_i] || j-1 > dp_end[pre_i]) continue;                                   \
                _pre_dp_h = (score_t*)(DP_H + pre_i * dp_sn);                                               \
                if (_pre_dp_h[j-1] + s == _dp_h[j]) { /* match/mismatch */                                  \
                    cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, j-1);                  \
                    i = pre_i; --j; hit = 1; id = abpoa_graph_index_to_node_id(graph, i+beg_index);         \
                    dp_h = DP_H + i * dp_sn; _dp_h = (score_t*)dp_h;                                        \
                    ++res->n_aln_bases; res->n_matched_bases += is_match ? 1 : 0;                           \
                    indel_first = 0;                                                                        \
                    break;                                                                                  \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0) err_fatal_simple("Error in lg_backtrack.");                                           \
     /* fprintf(stderr, "%d, %d, (%d)\n", i, j, indel_first); */                                            \
    }                                                                                                       \
    if (j > 0) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, j, -1, j-1);                         \
    /* reverse cigar */                                                                                     \
    res->graph_cigar = abpt->rev_cigar ? cigar : abpoa_reverse_cigar(n_c, cigar);                           \
    res->n_cigar = n_c; res->m_cigar = m_c;                                                                 \
    res->node_e = abpoa_graph_index_to_node_id(graph, best_i+beg_index), res->query_e=best_j-1; /*0-based*/ \
    res->node_s = abpoa_graph_index_to_node_id(graph, _start_i+beg_index), res->query_s=_start_j-1;         \
    /*abpoa_print_cigar(n_c, *graph_cigar, graph);*/                                                        \
}

#define simd_abpoa_ag_backtrack(score_t) {                                                                  \
    int i, j, k, pre_i, n_c = 0, s, is_match, m_c = 0, id, hit, cur_op = ABPOA_ALL_OP, _start_i, _start_j;  \
    score_t *_dp_h, *_dp_e1, *_dp_f1, *_pre_dp_h, *_pre_dp_e1; abpoa_cigar_t *cigar = 0;                    \
    i = best_i, j = best_j; _start_i = best_i, _start_j = best_j;                                           \
    id = abpoa_graph_index_to_node_id(graph, i+beg_index);                                                  \
    if (best_j < qlen) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, qlen-j, -1, qlen-1);         \
    SIMDi *dp_h = DP_HEF + dp_sn * i * 3; _dp_h = (score_t*)dp_h;                                           \
    int indel_first = 1; /* prefer to keep gaps at the end */                                               \
    while (i > 0 && j > 0) {                                                                                \
        if (abpt->align_mode == ABPOA_LOCAL_MODE && _dp_h[j] == 0) break;                                   \
        _start_i = i, _start_j = j;                                                                         \
        int *pre_index_i = pre_index[i];                                                                    \
        s = mat[m * graph->node[id].base + query[j-1]]; hit = 0;                                            \
        is_match = graph->node[id].base == query[j-1];                                                      \
        if (cur_op & ABPOA_M_OP && indel_first == 0) { /* match/mismatch */                                 \
            for (k = 0; k < pre_n[i]; ++k) {                                                                \
                pre_i = pre_index_i[k];                                                                     \
                if (j-1 < dp_beg[pre_i] || j-1 > dp_end[pre_i]) continue;                                   \
                _pre_dp_h = (score_t*)(DP_HEF + dp_sn * pre_i * 3);                                         \
                if (_pre_dp_h[j-1] + s == _dp_h[j]) {                                                       \
                    cur_op = ABPOA_ALL_OP; hit = 1;                                                         \
                    cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, j-1);                  \
                    i = pre_i; --j; id = abpoa_graph_index_to_node_id(graph, i+beg_index);                  \
                    dp_h = DP_HEF + dp_sn * i * 3; _dp_h = (score_t*)dp_h;                                  \
                    ++res->n_aln_bases; res->n_matched_bases += is_match ? 1 : 0;                           \
                    break;                                                                                  \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0 && cur_op & ABPOA_E1_OP) { /* deletion */                                              \
            for (k = 0; k < pre_n[i]; ++k) {                                                                \
                pre_i = pre_index_i[k];                                                                     \
                if (j < dp_beg[pre_i] || j > dp_end[pre_i]) continue;                                       \
                _pre_dp_e1 = (score_t*)(DP_HEF + dp_sn * (pre_i * 3 + 1));                                  \
                if (cur_op & ABPOA_M_OP) {                                                                  \
                    if (_dp_h[j] == _pre_dp_e1[j]) {                                                        \
                        _pre_dp_h = (score_t*)(DP_HEF + dp_sn * pre_i * 3);                                 \
                        if (_pre_dp_h[j] - gap_oe1 == _pre_dp_e1[j]) cur_op = ABPOA_M_OP | ABPOA_F_OP;      \
                        else cur_op = ABPOA_E1_OP;                                                          \
                        hit = 1;                                                                            \
                        cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);                \
                        i = pre_i; id = abpoa_graph_index_to_node_id(graph, i+beg_index);                   \
                        dp_h = DP_HEF + dp_sn * i * 3; _dp_h = (score_t*)dp_h;                              \
                        break;                                                                              \
                    }                                                                                       \
                } else {                                                                                    \
                    _dp_e1 = (score_t*)(dp_h + dp_sn);                                                      \
                    if (_dp_e1[j] == _pre_dp_e1[j] - gap_ext1) {                                            \
                        _pre_dp_h = (score_t*)(DP_HEF + dp_sn * pre_i * 3);                                 \
                        if (_pre_dp_h[j] - gap_oe1 == _pre_dp_e1[j]) cur_op = ABPOA_M_OP | ABPOA_F_OP;      \
                        else cur_op = ABPOA_E1_OP;                                                          \
                        hit = 1;                                                                            \
                        cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);                \
                        i = pre_i; id = abpoa_graph_index_to_node_id(graph, i+beg_index);                   \
                        dp_h = DP_HEF + dp_sn * i * 3; _dp_h = (score_t*)dp_h;                              \
                        break;                                                                              \
                    }                                                                                       \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0 && cur_op & ABPOA_F_OP) { /* insertion */                                              \
            _dp_f1 = (score_t*)(dp_h + dp_sn * 2);                                                          \
            if (cur_op & ABPOA_M_OP) {                                                                      \
                if (_dp_h[j] == _dp_f1[j]) {                                                                \
                    if (_dp_h[j-1] - gap_oe1 == _dp_f1[j]) cur_op = ABPOA_M_OP | ABPOA_E_OP, hit = 1;       \
                    else if (_dp_f1[j-1] - gap_ext1 == _dp_f1[j]) cur_op = ABPOA_F1_OP, hit = 1;            \
                }                                                                                           \
            } else {                                                                                        \
                if (_dp_h[j-1] - gap_oe1 == _dp_f1[j]) cur_op = ABPOA_M_OP | ABPOA_E_OP, hit = 1;           \
                else if (_dp_f1[j-1] - gap_ext1 == _dp_f1[j]) cur_op = ABPOA_F1_OP, hit = 1;                \
            }                                                                                               \
            if (hit == 1) {                                                                                 \
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, 1, id, j-1); --j;                   \
                ++res->n_aln_bases;                                                                         \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0 && cur_op & ABPOA_M_OP && indel_first == 1) {                                          \
            for (k = 0; k < pre_n[i]; ++k) {                                                                \
                pre_i = pre_index_i[k];                                                                     \
                if (j-1 < dp_beg[pre_i] || j-1 > dp_end[pre_i]) continue;                                   \
                _pre_dp_h = (score_t*)(DP_HEF + dp_sn * pre_i * 3);                                         \
                if (_pre_dp_h[j-1] + s == _dp_h[j]) {                                                       \
                    cur_op = ABPOA_ALL_OP; hit = 1;                                                         \
                    cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, j-1);                  \
                    i = pre_i; --j; id = abpoa_graph_index_to_node_id(graph, i+beg_index);                  \
                    dp_h = DP_HEF + dp_sn * i * 3; _dp_h = (score_t*)dp_h;                                  \
                    ++res->n_aln_bases; res->n_matched_bases += is_match ? 1 : 0;                           \
                    indel_first = 0;                                                                        \
                    break;                                                                                  \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0) err_fatal_simple("Error in ag_backtrack.");                                           \
     /* fprintf(stderr, "%d, %d, %d (indel_first: %d)\n", i, j, cur_op, indel_first); */                    \
    }                                                                                                       \
    if (j > 0) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, j, -1, j-1);                         \
    /* reverse cigar */                                                                                     \
    res->graph_cigar = abpt->rev_cigar ? cigar : abpoa_reverse_cigar(n_c, cigar);                           \
    res->n_cigar = n_c; res->m_cigar = m_c;                                                                 \
    res->node_e = abpoa_graph_index_to_node_id(graph, best_i+beg_index), res->query_e=best_j-1; /*0-based*/ \
    res->node_s = abpoa_graph_index_to_node_id(graph, _start_i+beg_index), res->query_s=_start_j-1;         \
    /*abpoa_print_cigar(n_c, *graph_cigar, graph);*/                                                        \
}

#define simd_abpoa_cg_backtrack(score_t) {                                                                  \
    int i, j, k, pre_i, n_c = 0, s, is_match, m_c = 0, id, hit, cur_op = ABPOA_ALL_OP, _start_i, _start_j;  \
    score_t *_dp_h, *_dp_e1, *_dp_e2, *_dp_f1, *_dp_f2, *_pre_dp_h, *_pre_dp_e1, *_pre_dp_e2;               \
    abpoa_cigar_t *cigar = 0;                                                                               \
    i = best_i, j = best_j, _start_i = best_i, _start_j = best_j;                                           \
    id = abpoa_graph_index_to_node_id(graph, i+beg_index);                                                  \
    if (best_j < qlen) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, qlen-j, -1, qlen-1);         \
    SIMDi *dp_h = DP_H2E2F + dp_sn * i * 5; _dp_h = (score_t*)dp_h;                                         \
    int indel_first = 1; /* prefer to keep gaps at the end */                                               \
    while (i > 0 && j > 0) {                                                                                \
        if (abpt->align_mode == ABPOA_LOCAL_MODE && _dp_h[j] == 0) break;                                   \
        _start_i = i, _start_j = j;                                                                         \
        int *pre_index_i = pre_index[i];                                                                    \
        s = mat[m * graph->node[id].base + query[j-1]]; hit = 0;                                            \
        is_match = graph->node[id].base == query[j-1];                                                      \
        if (cur_op & ABPOA_M_OP && indel_first == 0) { /* match/mismatch */                                 \
            for (k = 0; k < pre_n[i]; ++k) {                                                                \
                pre_i = pre_index_i[k];                                                                     \
                if (j-1 < dp_beg[pre_i] || j-1 > dp_end[pre_i]) continue;                                   \
                _pre_dp_h = (score_t*)(DP_H2E2F + dp_sn * pre_i * 5);                                       \
                if (_pre_dp_h[j-1] + s == _dp_h[j]) {                                                       \
                    cur_op = ABPOA_ALL_OP; hit = 1;                                                         \
                    cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, j-1);                  \
                    i = pre_i; --j; id = abpoa_graph_index_to_node_id(graph, i+beg_index); hit = 1;         \
                    dp_h = DP_H2E2F + dp_sn * i * 5; _dp_h = (score_t*)dp_h;                                \
                    ++res->n_aln_bases; res->n_matched_bases += is_match ? 1 : 0;                           \
                    break;                                                                                  \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0 && cur_op & ABPOA_E_OP) { /* deletion */                                               \
            for (k = 0; k < pre_n[i]; ++k) {                                                                \
                pre_i = pre_index_i[k];                                                                     \
                if (j < dp_beg[pre_i] || j > dp_end[pre_i]) continue;                                       \
                if (cur_op & ABPOA_E1_OP) {                                                                 \
                    _pre_dp_e1 = (score_t*)(DP_H2E2F + dp_sn * (pre_i * 5 + 1));                            \
                    if (cur_op & ABPOA_M_OP) {                                                              \
                        if (_dp_h[j] == _pre_dp_e1[j]) {                                                    \
                            _pre_dp_h = (score_t*)(DP_H2E2F + dp_sn * pre_i * 5);                           \
                            if (_pre_dp_h[j] - gap_oe1 == _pre_dp_e1[j]) cur_op = ABPOA_M_OP | ABPOA_F_OP;  \
                            else cur_op = ABPOA_E1_OP;                                                      \
                            hit = 1; cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);   \
                            i = pre_i; id = abpoa_graph_index_to_node_id(graph, i+beg_index);               \
                            dp_h = DP_H2E2F + dp_sn * i * 5; _dp_h = (score_t*)dp_h;                        \
                            break;                                                                          \
                        }                                                                                   \
                    } else {                                                                                \
                        _dp_e1 = (score_t*)(dp_h + dp_sn);                                                  \
                        if (_dp_e1[j] == _pre_dp_e1[j] - gap_ext1) {                                        \
                            _pre_dp_h = (score_t*)(DP_H2E2F + dp_sn * pre_i * 5);                           \
                            if (_pre_dp_h[j] - gap_oe1 == _pre_dp_e1[j]) cur_op = ABPOA_M_OP | ABPOA_F_OP;  \
                            else cur_op = ABPOA_E1_OP;                                                      \
                            hit = 1; cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);   \
                            i = pre_i; id = abpoa_graph_index_to_node_id(graph, i+beg_index);               \
                            dp_h = DP_H2E2F + dp_sn * i * 5; _dp_h = (score_t*)dp_h;                        \
                            break;                                                                          \
                        }                                                                                   \
                    }                                                                                       \
                }                                                                                           \
                if (cur_op & ABPOA_E2_OP) {                                                                 \
                    _pre_dp_e2 = (score_t*)(DP_H2E2F + dp_sn * (pre_i * 5 + 2));                            \
                    if (cur_op & ABPOA_M_OP) {                                                              \
                        if (_dp_h[j] == _pre_dp_e2[j]) {                                                    \
                            _pre_dp_h = (score_t*)(DP_H2E2F + dp_sn * pre_i * 5);                           \
                            if (_pre_dp_h[j] - gap_oe2 == _pre_dp_e2[j]) cur_op = ABPOA_M_OP | ABPOA_F_OP;  \
                            else cur_op = ABPOA_E2_OP;                                                      \
                            hit = 1; cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);   \
                            i = pre_i; id = abpoa_graph_index_to_node_id(graph, i+beg_index);               \
                            dp_h = DP_H2E2F + dp_sn * i * 5; _dp_h = (score_t*)dp_h;                        \
                            break;                                                                          \
                        }                                                                                   \
                    } else {                                                                                \
                        _dp_e2 = (score_t*)(dp_h + dp_sn * 2);                                              \
                        if (_dp_e2[j] == _pre_dp_e2[j] - gap_ext2) {                                        \
                            _pre_dp_h = (score_t*)(DP_H2E2F + dp_sn * pre_i * 5);                           \
                            if (_pre_dp_h[j] - gap_oe2 == _pre_dp_e2[j]) cur_op = ABPOA_M_OP | ABPOA_F_OP;  \
                            else cur_op = ABPOA_E2_OP;                                                      \
                            hit = 1; cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, j-1);   \
                            i = pre_i; id = abpoa_graph_index_to_node_id(graph, i+beg_index);               \
                            dp_h = DP_H2E2F + dp_sn * i * 5; _dp_h = (score_t*)dp_h;                        \
                            break;                                                                          \
                        }                                                                                   \
                    }                                                                                       \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0 && cur_op & ABPOA_F_OP) { /* insertion */                                              \
            if (cur_op & ABPOA_F1_OP) {                                                                     \
                _dp_f1 = (score_t*)(dp_h + dp_sn * 3);                                                      \
                if (cur_op & ABPOA_M_OP) {                                                                  \
                    if (_dp_h[j] == _dp_f1[j]) {                                                            \
                        if (_dp_h[j-1] - gap_oe1 == _dp_f1[j]) cur_op = ABPOA_M_OP | ABPOA_E_OP, hit = 1;   \
                        else if (_dp_f1[j-1] - gap_ext1 == _dp_f1[j]) cur_op = ABPOA_F1_OP, hit = 1;        \
                    }                                                                                       \
                } else {                                                                                    \
                    if (_dp_h[j-1] - gap_oe1 == _dp_f1[j]) cur_op = ABPOA_M_OP | ABPOA_E_OP, hit = 1;       \
                    else if (_dp_f1[j-1] - gap_ext1 == _dp_f1[j]) cur_op = ABPOA_F1_OP, hit = 1;            \
                }                                                                                           \
            }                                                                                               \
            if (hit == 0 && cur_op & ABPOA_F2_OP) {                                                         \
                _dp_f2 = (score_t*)(dp_h + dp_sn * 4);                                                      \
                if (cur_op & ABPOA_M_OP) {                                                                  \
                    if (_dp_h[j] == _dp_f2[j]) {                                                            \
                        if (_dp_h[j-1] - gap_oe2 == _dp_f2[j]) cur_op = ABPOA_M_OP | ABPOA_E_OP, hit = 1;   \
                        else if (_dp_f2[j-1] - gap_ext2 == _dp_f2[j]) cur_op = ABPOA_F2_OP, hit = 1;        \
                    }                                                                                       \
                } else {                                                                                    \
                    if (_dp_h[j-1] - gap_oe2 == _dp_f2[j]) cur_op = ABPOA_M_OP | ABPOA_E_OP, hit = 1;       \
                    else if (_dp_f2[j-1] - gap_ext2 == _dp_f2[j]) cur_op = ABPOA_F2_OP, hit = 1;            \
                }                                                                                           \
            }                                                                                               \
            if (hit == 1) {                                                                                 \
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, 1, id, j-1); --j;                   \
                ++res->n_aln_bases;                                                                         \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0 && cur_op & ABPOA_M_OP && indel_first == 1) { /* match/mismatch */                     \
            for (k = 0; k < pre_n[i]; ++k) {                                                                \
                pre_i = pre_index_i[k];                                                                     \
                if (j-1 < dp_beg[pre_i] || j-1 > dp_end[pre_i]) continue;                                   \
                _pre_dp_h = (score_t*)(DP_H2E2F + dp_sn * pre_i * 5);                                       \
                if (_pre_dp_h[j-1] + s == _dp_h[j]) {                                                       \
                    cur_op = ABPOA_ALL_OP; hit = 1;                                                         \
                    cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, j-1);                  \
                    i = pre_i; --j; id = abpoa_graph_index_to_node_id(graph, i+beg_index); hit = 1;         \
                    dp_h = DP_H2E2F + dp_sn * i * 5; _dp_h = (score_t*)dp_h;                                \
                    ++res->n_aln_bases; res->n_matched_bases += is_match ? 1 : 0;                           \
                    indel_first = 0;                                                                        \
                    break;                                                                                  \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
        if (hit == 0) err_fatal_simple("Error in cg_backtrack.");                                           \
     /* fprintf(stderr, "%d, %d, %d\n", i, j, cur_op); */                                                   \
    }                                                                                                       \
    if (j > 0) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, j, -1, j-1);                         \
    /* reverse cigar */                                                                                     \
    res->graph_cigar = abpt->rev_cigar ? cigar : abpoa_reverse_cigar(n_c, cigar);                           \
    res->n_cigar = n_c; res->m_cigar = m_c;                                                                 \
    res->node_e = abpoa_graph_index_to_node_id(graph, best_i+beg_index), res->query_e=best_j-1; /*0-based*/ \
    res->node_s = abpoa_graph_index_to_node_id(graph, _start_i+beg_index), res->query_s=_start_j-1;         \
    /*abpoa_print_cigar(n_c, *graph_cigar, graph);*/                                                        \
}

// simd_abpoa_va
// simd_abpoa_ag_only_var
// sim_abpoa_init_var
#define simd_abpoa_var(score_t, sp, SIMDSetOne, SIMDShiftOneN)                                      \
    /* int tot_dp_sn = 0; */                                                                        \
    abpoa_graph_t *graph = ab->abg; abpoa_simd_matrix_t *abm = ab->abm;                             \
    int matrix_row_n = end_index-beg_index+1, matrix_col_n = qlen + 1;                              \
    int **pre_index, *pre_n, _pre_index, _pre_n, pre_i;                                             \
    int i, j, k, *dp_beg, *dp_beg_sn, *dp_end, *dp_end_sn, node_id, index_i, dp_i;                  \
    int beg, end, beg_sn, end_sn, _beg_sn, _end_sn, pre_beg_sn, pre_end, sn_i;                      \
    int pn, log_n, size, qp_sn, dp_sn; /* pn: value per SIMDi, qp_sn/dp_sn/d_sn: segmented length*/ \
    SIMDi *dp_h, *pre_dp_h, *qp, *qi=NULL;                                                          \
    score_t *_dp_h=NULL, *_qi, best_score = sp.inf_min, inf_min = sp.inf_min;                       \
    int *mat = abpt->mat, m = abpt->m; score_t gap_ext1 = abpt->gap_ext1;                           \
    int w = abpt->wb < 0 ? qlen : abpt->wb+(int)(abpt->wf*qlen); /* when w < 0, do whole global */  \
    int best_i = 0, best_j = 0, best_id = 0, max, max_i=-1;                                         \
    SIMDi zero = SIMDSetZeroi(), SIMD_INF_MIN = SIMDSetOne(inf_min);                                \
    pn = sp.num_of_value; qp_sn = dp_sn = (matrix_col_n + pn - 1) / pn;                             \
    log_n = sp.log_num, size = sp.size; qp = abm->s_mem;                                            \
    int set_num; SIMDi *PRE_MASK, *SUF_MIN, *PRE_MIN;                                               \
    PRE_MASK = (SIMDi*)SIMDMalloc((pn+1) * size, size);                                             \
    SUF_MIN = (SIMDi*)SIMDMalloc((pn+1) * size, size);                                              \
    PRE_MIN = (SIMDi*)SIMDMalloc(pn * size, size);                                                  \
    for (i = 0; i < pn; ++i) {                                                                      \
        score_t *pre_mask = (score_t*)(PRE_MASK+i);                                                 \
        for (j = 0; j <= i; ++j) pre_mask[j] = -1;                                                  \
        for (j = i+1; j < pn; ++j) pre_mask[j] = 0;                                                 \
    } PRE_MASK[pn] = PRE_MASK[pn-1];                                                                \
    SUF_MIN[0] = SIMDShiftLeft(SIMD_INF_MIN, SIMDShiftOneN);                                        \
    for (i = 1; i < pn; ++i)                                                                        \
        SUF_MIN[i] = SIMDShiftLeft(SUF_MIN[i-1], SIMDShiftOneN); SUF_MIN[pn] = SUF_MIN[pn-1];       \
    for (i = 1; i < pn; ++i) {                                                                      \
        score_t *pre_min = (score_t*)(PRE_MIN + i);                                                 \
        for (j = 0; j < i; ++j) pre_min[j] = inf_min;                                               \
        for (j = i; j < pn; ++j) pre_min[j] = 0;                                                    \
    }

#define simd_abpoa_lg_only_var(score_t, SIMDSetOne, SIMDAdd)              \
    SIMDi *DP_H = qp + qp_sn * abpt->m; qi = DP_H + dp_sn * matrix_row_n; \
    SIMDi GAP_E1 = SIMDSetOne(gap_ext1);                                  \
    SIMDi *GAP_E1S =  (SIMDi*)SIMDMalloc(log_n * size, size);             \
    GAP_E1S[0] = GAP_E1;                                                  \
    for (i = 1; i < log_n; ++i) {                                         \
        GAP_E1S[i] = SIMDAdd(GAP_E1S[i-1], GAP_E1S[i-1]);                 \
    }

#define simd_abpoa_ag_only_var(score_t, SIMDSetOne, SIMDAdd)                                            \
    score_t *_dp_e1, *_dp_f1, gap_open1 = abpt->gap_open1, gap_oe1 = abpt->gap_open1 + abpt->gap_ext1;  \
    SIMDi *DP_HEF, *dp_e1, *pre_dp_e1, *dp_f1; int pre_end_sn;                                          \
    DP_HEF = qp + qp_sn * abpt->m; qi = DP_HEF + dp_sn * matrix_row_n * 3;                              \
    SIMDi GAP_O1 = SIMDSetOne(gap_open1), GAP_E1 = SIMDSetOne(gap_ext1), GAP_OE1 = SIMDSetOne(gap_oe1); \
    SIMDi *GAP_E1S =  (SIMDi*)SIMDMalloc(log_n * size, size);  GAP_E1S[0] = GAP_E1;                     \
    for (i = 1; i < log_n; ++i) {                                                                       \
        GAP_E1S[i] = SIMDAdd(GAP_E1S[i-1], GAP_E1S[i-1]);                                               \
    }

#define simd_abpoa_cg_only_var(score_t, SIMDSetOne, SIMDAdd)                                                      \
    score_t *_dp_e1, *_dp_e2, *_dp_f1, *_dp_f2, gap_open1 = abpt->gap_open1, gap_oe1 = gap_open1 + gap_ext1;      \
    score_t gap_open2 = abpt->gap_open2, gap_ext2 = abpt->gap_ext2, gap_oe2 = gap_open2 + gap_ext2;               \
    SIMDi *DP_H2E2F, *dp_e1, *dp_e2, *dp_f1, *dp_f2, *pre_dp_e1, *pre_dp_e2; int pre_end_sn;                      \
    SIMDi GAP_O1 = SIMDSetOne(gap_open1), GAP_O2 = SIMDSetOne(gap_open2);                                         \
    SIMDi GAP_E1 = SIMDSetOne(gap_ext1), GAP_E2 = SIMDSetOne(gap_ext2);                                           \
    SIMDi GAP_OE1 = SIMDSetOne(gap_oe1), GAP_OE2 = SIMDSetOne(gap_oe2);                                           \
    DP_H2E2F = qp + qp_sn * abpt->m; qi = DP_H2E2F + dp_sn * matrix_row_n * 5;                                    \
    SIMDi *GAP_E1S =  (SIMDi*)SIMDMalloc(log_n * size, size), *GAP_E2S =  (SIMDi*)SIMDMalloc(log_n * size, size); \
    GAP_E1S[0] = GAP_E1; GAP_E2S[0] = GAP_E2;                                                                     \
    for (i = 1; i < log_n; ++i) {                                                                                 \
        GAP_E1S[i] = SIMDAdd(GAP_E1S[i-1], GAP_E1S[i-1]);                                                         \
        GAP_E2S[i] = SIMDAdd(GAP_E2S[i-1], GAP_E2S[i-1]);                                                         \
    }

#define simd_abpoa_init_var(score_t) {                                                             \
    /* generate the query profile */                                                               \
    for (i = 0; i < qp_sn * abpt->m; ++i) qp[i] = SIMD_INF_MIN;                                    \
    for (k = 0; k < abpt->m; ++k) { /* SIMD parallelization */                                     \
        int *p = &mat[k * abpt->m];                                                                \
        score_t *_qp = (score_t*)(qp + k * qp_sn); _qp[0] = 0;                                     \
        for (j = 0; j < qlen; ++j) _qp[j+1] = (score_t)p[query[j]];                                \
        for (j = qlen+1; j < qp_sn * pn; ++j) _qp[j] = 0;                                          \
    }                                                                                              \
    if (abpt->wb>=0 || abpt->align_mode==ABPOA_LOCAL_MODE || abpt->align_mode==ABPOA_EXTEND_MODE){ \
        _qi = (score_t*)qi; /* query index */                                                      \
        for (i = 0; i <= qlen; ++i) _qi[i] = i;                                                    \
        for (i = qlen+1; i < (qlen/pn+1) * pn; ++i) _qi[i] = -1;                                   \
    }                                                                                              \
    /* for backtrack */                                                                            \
    dp_beg=abm->dp_beg, dp_end=abm->dp_end, dp_beg_sn=abm->dp_beg_sn, dp_end_sn=abm->dp_end_sn;    \
    /* index of pre-node */                                                                        \
    pre_index = (int**)_err_calloc(matrix_row_n, sizeof(int*));                                    \
    pre_n = (int*)_err_calloc(matrix_row_n, sizeof(int));                                          \
    for (index_i=beg_index+1, dp_i=1; index_i<=end_index; ++index_i, ++dp_i) {                     \
        node_id = abpoa_graph_index_to_node_id(graph, index_i);                                    \
        pre_n[dp_i] = graph->node[node_id].in_edge_n;                                              \
        pre_index[dp_i] = (int*)_err_malloc(pre_n[dp_i] * sizeof(int));                            \
        for (j = _pre_n = 0; j < pre_n[dp_i]; ++j) {                                               \
            _pre_index = abpoa_graph_node_id_to_index(graph, graph->node[node_id].in_id[j]);       \
            if (index_map[_pre_index]) pre_index[dp_i][_pre_n++] = _pre_index-beg_index;           \
        }                                                                                          \
        pre_n[dp_i] = _pre_n;                                                                      \
    }                                                                                              \
}

#define simd_abpoa_free_var {                                                            \
    for (i = 0; i < matrix_row_n; ++i) free(pre_index[i]); free(pre_index); free(pre_n); \
    SIMDFree(PRE_MASK); SIMDFree(SUF_MIN); SIMDFree(PRE_MIN);                            \
}                                                                                        \

#define simd_abpoa_lg_var(score_t, sp, SIMDSetOne, SIMDShiftOneN, SIMDAdd) \
    simd_abpoa_var(score_t, sp, SIMDSetOne, SIMDShiftOneN);                \
    simd_abpoa_lg_only_var(score_t, SIMDSetOne, SIMDAdd);                  \
    simd_abpoa_init_var(score_t);

#define simd_abpoa_ag_var(score_t, sp, SIMDSetOne, SIMDShiftOneN, SIMDAdd) \
    simd_abpoa_var(score_t, sp, SIMDSetOne, SIMDShiftOneN);                \
    simd_abpoa_ag_only_var(score_t, SIMDSetOne, SIMDAdd);                  \
    simd_abpoa_init_var(score_t);

#define simd_abpoa_cg_var(score_t, sp, SIMDSetOne, SIMDShiftOneN, SIMDAdd) \
    simd_abpoa_var(score_t, sp, SIMDSetOne, SIMDShiftOneN);                \
    simd_abpoa_cg_only_var(score_t, SIMDSetOne, SIMDAdd);                  \
    simd_abpoa_init_var(score_t);

#define simd_abpoa_lg_first_row {                                                                       \
    /* fill the first row */                                                                            \
    if (abpt->wb >= 0) {                                                                                \
        graph->node_id_to_max_pos_left[beg_node_id] = graph->node_id_to_max_pos_right[beg_node_id] = 0; \
        for (i = 0; i < graph->node[beg_node_id].out_edge_n; ++i) { /* set max pos for out_id */        \
            int out_id = graph->node[beg_node_id].out_id[i];                                            \
            if (index_map[abpoa_graph_node_id_to_index(graph, out_id)])                                 \
                graph->node_id_to_max_pos_left[out_id] = graph->node_id_to_max_pos_right[out_id] = 1;   \
        }                                                                                               \
        dp_beg[0] = 0, dp_end[0] = GET_AD_DP_END(graph, w, beg_node_id, end_node_id, qlen);             \
    } else {                                                                                            \
        dp_beg[0] = 0, dp_end[0] = qlen;                                                                \
    }                                                                                                   \
    dp_beg_sn[0] = (dp_beg[0])/pn; dp_end_sn[0] = (dp_end[0])/pn;                                       \
    dp_beg[0] = dp_beg_sn[0] * pn; dp_end[0] = (dp_end_sn[0]+1)*pn-1;                                   \
    dp_h = DP_H; _end_sn = MIN_OF_TWO(dp_end_sn[0]+1, dp_sn-1);                                         \
}

#define simd_abpoa_ag_first_row {                                                                       \
    /* fill the first row */                                                                            \
    if (abpt->wb >= 0) {                                                                                \
        graph->node_id_to_max_pos_left[beg_node_id] = graph->node_id_to_max_pos_right[beg_node_id] = 0; \
        for (i = 0; i < graph->node[beg_node_id].out_edge_n; ++i) { /* set max pos for out_id */        \
            int out_id = graph->node[beg_node_id].out_id[i];                                            \
            if (index_map[abpoa_graph_node_id_to_index(graph, out_id)])                                 \
                graph->node_id_to_max_pos_left[out_id] = graph->node_id_to_max_pos_right[out_id] = 1;   \
        }                                                                                               \
        dp_beg[0] = 0, dp_end[0] = GET_AD_DP_END(graph, w, beg_node_id, end_node_id, qlen);             \
    } else {                                                                                            \
        dp_beg[0] = 0, dp_end[0] = qlen;                                                                \
    }                                                                                                   \
    dp_beg_sn[0] = (dp_beg[0])/pn; dp_end_sn[0] = (dp_end[0])/pn;                                       \
    dp_beg[0] = dp_beg_sn[0] * pn; dp_end[0] = (dp_end_sn[0]+1)*pn-1;                                   \
    dp_h = DP_HEF; dp_e1 = dp_h + dp_sn; dp_f1 = dp_e1 + dp_sn;                                         \
    _end_sn = MIN_OF_TWO(dp_end_sn[0]+1, dp_sn-1);                                                      \
}

#define simd_abpoa_cg_first_row {                                                                       \
    /* fill the first row */                                                                            \
    if (abpt->wb >= 0) {                                                                                \
        graph->node_id_to_max_pos_left[beg_node_id] = graph->node_id_to_max_pos_right[beg_node_id] = 0; \
        for (i = 0; i < graph->node[beg_node_id].out_edge_n; ++i) { /* set max pos for out_id */        \
            int out_id = graph->node[beg_node_id].out_id[i];                                            \
            if (index_map[abpoa_graph_node_id_to_index(graph, out_id)])                                 \
                graph->node_id_to_max_pos_left[out_id] = graph->node_id_to_max_pos_right[out_id] = 1;   \
        }                                                                                               \
        dp_beg[0] = 0, dp_end[0] = GET_AD_DP_END(graph, w, beg_node_id, end_node_id, qlen);             \
    } else {                                                                                            \
        dp_beg[0] = 0, dp_end[0] = qlen;                                                                \
    }                                                                                                   \
    dp_beg_sn[0] = (dp_beg[0])/pn; dp_end_sn[0] = (dp_end[0])/pn;                                       \
    dp_beg[0] = dp_beg_sn[0] * pn; dp_end[0] = (dp_end_sn[0]+1)*pn-1;                                   \
    dp_h = DP_H2E2F; dp_e1 = dp_h+dp_sn; dp_e2 = dp_e1+dp_sn; dp_f1 = dp_e2+dp_sn; dp_f2 = dp_f1+dp_sn; \
    _end_sn = MIN_OF_TWO(dp_end_sn[0]+1, dp_sn-1);                                                      \
}

#define simd_abpoa_lg_first_dp(score_t) {                                   \
    simd_abpoa_lg_first_row;                                                \
    if (abpt->align_mode == ABPOA_LOCAL_MODE) {                             \
        for (i = 0; i < _end_sn; ++i)                                       \
            dp_h[i] = zero;                                                 \
    } else {                                                                \
        for (i = 0; i <= _end_sn; ++i) {                                    \
            dp_h[i] = SIMD_INF_MIN;                                         \
        }                                                                   \
        _dp_h = (score_t*)dp_h;                                             \
        for (i = 0; i <= dp_end[0]; ++i) { /* no SIMD parallelization */    \
            _dp_h[i] = -gap_ext1 * i;                                       \
        }                                                                   \
    }                                                                       \
}

#define simd_abpoa_ag_first_dp(score_t) {                                   \
    simd_abpoa_ag_first_row;                                                \
    if (abpt->align_mode == ABPOA_LOCAL_MODE) {                             \
        for (i = 0; i < _end_sn; ++i)                                       \
            dp_h[i] = dp_e1[i] = dp_f1[i] = zero;                           \
    } else {                                                                \
        for (i = 0; i <= _end_sn; ++i) {                                    \
            dp_h[i] = SIMD_INF_MIN; dp_e1[i] = SIMD_INF_MIN;                \
        }                                                                   \
        _dp_h=(score_t*)dp_h,_dp_e1=(score_t*)dp_e1,_dp_f1=(score_t*)dp_f1; \
        _dp_h[0] = 0; _dp_e1[0] = -(gap_oe1), _dp_f1[0] = inf_min;          \
        for (i = 1; i <= dp_end[0]; ++i) { /* no SIMD parallelization */    \
            _dp_f1[i] = -gap_open1 - gap_ext1 * i;                          \
            _dp_h[i] = -gap_open1 - gap_ext1 * i;                           \
        }                                                                   \
    }                                                                       \
}

#define simd_abpoa_cg_first_dp(score_t) {                                             \
    simd_abpoa_cg_first_row;                                                          \
    if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                       \
        for (i = 0; i < _end_sn; ++i)                                                 \
            dp_h[i] = dp_e1[i] = dp_e2[i] = dp_f1[i] = dp_f2[i] = zero;               \
    } else {                                                                          \
        for (i = 0; i <= _end_sn; ++i) {                                              \
            dp_h[i] = SIMD_INF_MIN; dp_e1[i] = SIMD_INF_MIN; dp_e2[i] = SIMD_INF_MIN; \
        }                                                                             \
        _dp_h = (score_t*)dp_h, _dp_e1 = (score_t*)dp_e1, _dp_e2 = (score_t*)dp_e2;   \
        _dp_f1 = (score_t*)dp_f1, _dp_f2 = (score_t*)dp_f2;                           \
        _dp_h[0] = 0; _dp_e1[0] = -(gap_oe1); _dp_e2[0] = -(gap_oe2);                 \
        _dp_f1[0] = _dp_f2[0] = inf_min;                                              \
        for (i = 1; i <= dp_end[0]; ++i) { /* no SIMD parallelization */              \
            _dp_f1[i] = -gap_open1 - gap_ext1 * i;                                    \
            _dp_f2[i] = -gap_open2 - gap_ext2 * i;                                    \
            _dp_h[i] = MAX_OF_TWO(_dp_f1[i], _dp_f2[i]);                              \
        }                                                                             \
    }                                                                                 \
}

// mask[pn], suf_min[pn], pre_min[logN]
#define SIMD_SET_F(F, log_n, set_num, PRE_MIN, PRE_MASK, SUF_MIN, GAP_E1S, SIMDMax, SIMDAdd, SIMDSub, SIMDShiftOneN) { \
    if (set_num == pn) {                                                                                               \
        F = SIMDMax(F, SIMDOri(SIMDShiftLeft(SIMDSub(F, GAP_E1S[0]), SIMDShiftOneN), PRE_MIN[1]));                     \
        if (log_n > 1) {                                                                                               \
            F = SIMDMax(F, SIMDOri(SIMDShiftLeft(SIMDSub(F, GAP_E1S[1]), SIMDShiftOneN<<1), PRE_MIN[2]));              \
        } if (log_n > 2) {                                                                                             \
            F = SIMDMax(F, SIMDOri(SIMDShiftLeft(SIMDSub(F, GAP_E1S[2]), SIMDShiftOneN<<2), PRE_MIN[4]));              \
        } if (log_n > 3) {                                                                                             \
            F = SIMDMax(F, SIMDOri(SIMDShiftLeft(SIMDSub(F, GAP_E1S[3]), SIMDShiftOneN<<3), PRE_MIN[8]));              \
        } if (log_n > 4) {                                                                                             \
            F = SIMDMax(F, SIMDOri(SIMDShiftLeft(SIMDSub(F, GAP_E1S[4]), SIMDShiftOneN<<4), PRE_MIN[16]));             \
        } if (log_n > 5) {                                                                                             \
            F = SIMDMax(F, SIMDOri(SIMDShiftLeft(SIMDSub(F, GAP_E1S[5]), SIMDShiftOneN<<5), PRE_MIN[32]));             \
        }                                                                                                              \
    } else { /*suffix MIN_INF*/                                                                                                                                    \
        int cov_bit = set_num;                                                                                                                                     \
        F = SIMDMax(F, SIMDOri(SIMDAndi(SIMDShiftLeft(SIMDSub(F, GAP_E1S[0]), SIMDShiftOneN), PRE_MASK[cov_bit]), SIMDOri(SUF_MIN[cov_bit], PRE_MIN[1])));         \
        if (log_n > 1) {                                                                                                                                           \
            cov_bit += 2;                                                                                                                                          \
            F = SIMDMax(F, SIMDOri(SIMDAndi(SIMDShiftLeft(SIMDSub(F, GAP_E1S[1]), SIMDShiftOneN<<1), PRE_MASK[cov_bit]), SIMDOri(SUF_MIN[cov_bit], PRE_MIN[2])));  \
        } if (log_n > 2) {                                                                                                                                         \
            cov_bit += 4;                                                                                                                                          \
            F = SIMDMax(F, SIMDOri(SIMDAndi(SIMDShiftLeft(SIMDSub(F, GAP_E1S[2]), SIMDShiftOneN<<2), PRE_MASK[cov_bit]), SIMDOri(SUF_MIN[cov_bit], PRE_MIN[4])));  \
        } if (log_n > 3) {                                                                                                                                         \
            cov_bit += 8;                                                                                                                                          \
            F = SIMDMax(F, SIMDOri(SIMDAndi(SIMDShiftLeft(SIMDSub(F, GAP_E1S[3]), SIMDShiftOneN<<3), PRE_MASK[cov_bit]), SIMDOri(SUF_MIN[cov_bit], PRE_MIN[8])));  \
        } if (log_n > 4) {                                                                                                                                         \
            cov_bit += 16;                                                                                                                                         \
            F = SIMDMax(F, SIMDOri(SIMDAndi(SIMDShiftLeft(SIMDSub(F, GAP_E1S[4]), SIMDShiftOneN<<4), PRE_MASK[cov_bit]), SIMDOri(SUF_MIN[cov_bit], PRE_MIN[16]))); \
        } if (log_n > 5) {                                                                                                                                         \
            cov_bit += 32;                                                                                                                                         \
            F = SIMDMax(F, SIMDOri(SIMDAndi(SIMDShiftLeft(SIMDSub(F, GAP_E1S[5]), SIMDShiftOneN<<5), PRE_MASK[cov_bit]), SIMDOri(SUF_MIN[cov_bit], PRE_MIN[32]))); \
        }                                                                                                                                                          \
    }                                                                                                                                                              \
}

#define simd_abpoa_lg_dp(score_t, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub) {                       \
    node_id = abpoa_graph_index_to_node_id(graph, index_i);                                         \
    SIMDi *q = qp + graph->node[node_id].base * qp_sn, first, remain;                               \
    dp_h = &DP_H[dp_i * dp_sn]; _dp_h = (score_t*)dp_h;                                             \
    int min_pre_beg_sn, max_pre_end_sn;                                                             \
    if (abpt->wb < 0) {                                                                             \
        beg = dp_beg[dp_i] = 0, end = dp_end[dp_i] = qlen;                                          \
        beg_sn = dp_beg_sn[dp_i] = (dp_beg[dp_i])/pn; end_sn = dp_end_sn[dp_i] = (dp_end[dp_i])/pn; \
        min_pre_beg_sn = 0, max_pre_end_sn = end_sn;                                                \
    } else {                                                                                        \
        beg = GET_AD_DP_BEGIN(graph, w, node_id, end_node_id, qlen), end = GET_AD_DP_END(graph, w, node_id, end_node_id, qlen); \
        beg_sn = beg / pn; min_pre_beg_sn = INT32_MAX, max_pre_end_sn = -1;                         \
        for (i = 0; i < pre_n[dp_i]; ++i) {                                                         \
            pre_i = pre_index[dp_i][i];                                                             \
            if (min_pre_beg_sn > dp_beg_sn[pre_i]) min_pre_beg_sn = dp_beg_sn[pre_i];               \
            if (max_pre_end_sn < dp_end_sn[pre_i]) max_pre_end_sn = dp_end_sn[pre_i];               \
        } if (beg_sn < min_pre_beg_sn) beg_sn = min_pre_beg_sn;                                     \
        dp_beg_sn[dp_i] = beg_sn; beg = dp_beg[dp_i] = dp_beg_sn[dp_i] * pn;                        \
        end_sn = dp_end_sn[dp_i] = end/pn; end = dp_end[dp_i] = (dp_end_sn[dp_i]+1)*pn-1;           \
    }                                                                                               \
    /* loop query */                                                                                                         \
    /* first pre_node */                                                                                                     \
    pre_i = pre_index[dp_i][0];                                                                                              \
    pre_dp_h = DP_H + pre_i * dp_sn;                                                                                         \
    pre_end = dp_end[pre_i];                                                                                                 \
    pre_beg_sn = dp_beg_sn[pre_i];                                                                                           \
    /* set M from (pre_i, q_i-1), E from (pre_i, q_i) */                                                                     \
    if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                              \
        _beg_sn = 0, _end_sn = end_sn; first = SIMDShiftRight(zero, SIMDTotalBytes-SIMDShiftOneN);                           \
    } else {                                                                                                                 \
        if (pre_beg_sn < beg_sn) _beg_sn = beg_sn, first = SIMDShiftRight(pre_dp_h[beg_sn-1], SIMDTotalBytes-SIMDShiftOneN); \
        else _beg_sn = pre_beg_sn, first = SIMDShiftRight(SIMD_INF_MIN, SIMDTotalBytes-SIMDShiftOneN);                       \
        _end_sn = MIN_OF_THREE((pre_end+1)/pn, end_sn, dp_sn-1);                                                             \
        for (i = beg_sn; i < _beg_sn; ++i) dp_h[i] = SIMD_INF_MIN;                                                           \
        for (i = _end_sn+1; i <= MIN_OF_TWO(end_sn+1, dp_sn-1); ++i) dp_h[i] = SIMD_INF_MIN;                                 \
    }                                                                                                                        \
    for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */                                               \
        remain = SIMDShiftLeft(pre_dp_h[sn_i], SIMDShiftOneN);                                                               \
        dp_h[sn_i] = SIMDMax(SIMDAdd(SIMDOri(first, remain), q[sn_i]), SIMDSub(pre_dp_h[sn_i], GAP_E1));                     \
        first = SIMDShiftRight(pre_dp_h[sn_i], SIMDTotalBytes-SIMDShiftOneN);                                                \
    }                                                                                                                        \
    /* get max m and e */                                                                                                         \
    for (i = 1; i < pre_n[dp_i]; ++i) {                                                                                           \
        pre_i = pre_index[dp_i][i];                                                                                               \
        pre_dp_h = DP_H + pre_i * dp_sn;                                                                                          \
        pre_end = dp_end[pre_i];                                                                                                  \
        pre_beg_sn = dp_beg_sn[pre_i];                                                                                            \
        /* set M from (pre_i, q_i-1), E from (pre_i, q_i) */                                                                      \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                               \
            first = SIMDShiftRight(zero, SIMDTotalBytes-SIMDShiftOneN);                                                           \
        } else {                                                                                                                  \
            if (pre_beg_sn < beg_sn) _beg_sn = beg_sn, first = SIMDShiftRight(pre_dp_h[beg_sn-1], SIMDTotalBytes-SIMDShiftOneN);  \
            else _beg_sn = pre_beg_sn, first = SIMDShiftRight(SIMD_INF_MIN, SIMDTotalBytes-SIMDShiftOneN);                        \
            _end_sn = MIN_OF_THREE((pre_end+1)/pn, end_sn, dp_sn-1);                                                              \
        }                                                                                                                         \
        for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */                                                \
            remain = SIMDShiftLeft(pre_dp_h[sn_i], SIMDShiftOneN);                                                                \
            dp_h[sn_i] = SIMDMax(SIMDAdd(SIMDOri(first, remain), q[sn_i]), SIMDMax(SIMDSub(pre_dp_h[sn_i], GAP_E1), dp_h[sn_i])); \
            first = SIMDShiftRight(pre_dp_h[sn_i], SIMDTotalBytes-SIMDShiftOneN);                                                 \
        } /* now we have max(h,e) stored at dp_h */                                                                               \
    }                                                                                                                             \
    /* new F start */                                                                                                                  \
    first = SIMDOri(SIMDAndi(dp_h[beg_sn], PRE_MASK[0]), SUF_MIN[0]);                                                                  \
    for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) {                                                                                      \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                                    \
            set_num = pn;                                                                                                              \
        } else {                                                                                                                       \
            if (sn_i < min_pre_beg_sn) {                                                                                               \
                _err_fatal_simple(__func__, "sn_i < min_pre_beg_sn\n");                                                                \
            } else if (sn_i > max_pre_end_sn) {                                                                                        \
                set_num = sn_i == max_pre_end_sn+1 ? 1 : 0;                                                                            \
            } else set_num = pn;                                                                                                       \
        }                                                                                                                              \
        dp_h[sn_i] = SIMDMax(dp_h[sn_i], first);                                                                                       \
        SIMD_SET_F(dp_h[sn_i], log_n, set_num, PRE_MIN, PRE_MASK, SUF_MIN, GAP_E1S, SIMDMax, SIMDAdd, SIMDSub, SIMDShiftOneN);         \
        first = SIMDOri(SIMDAndi(SIMDShiftRight(SIMDSub(dp_h[sn_i], GAP_E1), SIMDTotalBytes-SIMDShiftOneN), PRE_MASK[0]), SUF_MIN[0]); \
    }                                                                                                                                  \
    if (abpt->align_mode == ABPOA_LOCAL_MODE) for (sn_i = 0; sn_i <= end_sn; ++sn_i) dp_h[sn_i] = SIMDMax(zero, dp_h[sn_i]);           \
}

#define simd_abpoa_ag_dp(score_t, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub, SIMDGetIfGreater, SIMDSetIfGreater, SIMDSetIfEqual) { \
    node_id = abpoa_graph_index_to_node_id(graph, index_i);                                                                       \
    SIMDi *q = qp + graph->node[node_id].base * qp_sn, first, remain;                                                             \
    dp_h = DP_HEF + dp_i * 3 * dp_sn; dp_e1 = dp_h + dp_sn; dp_f1 = dp_e1 + dp_sn;                                                \
    _dp_h = (score_t*)dp_h, _dp_e1 = (score_t*)dp_e1, _dp_f1 = (score_t*)dp_f1;                                                   \
    int min_pre_beg_sn, max_pre_end_sn;                                                                                           \
    if (abpt->wb < 0) {                                                                                                           \
        beg = dp_beg[dp_i] = 0, end = dp_end[dp_i] = qlen;                                                                        \
        beg_sn = dp_beg_sn[dp_i] = (dp_beg[dp_i])/pn; end_sn = dp_end_sn[dp_i] = (dp_end[dp_i])/pn;                               \
        min_pre_beg_sn = 0, max_pre_end_sn = end_sn;                                                                              \
    } else {                                                                                                                      \
        beg = GET_AD_DP_BEGIN(graph, w, node_id, end_node_id, qlen), end = GET_AD_DP_END(graph, w, node_id, end_node_id, qlen);   \
        beg_sn = beg / pn; min_pre_beg_sn = INT32_MAX, max_pre_end_sn = -1;                                                       \
        for (i = 0; i < pre_n[dp_i]; ++i) {                                                                                       \
            pre_i = pre_index[dp_i][i];                                                                                           \
            if (min_pre_beg_sn > dp_beg_sn[pre_i]) min_pre_beg_sn = dp_beg_sn[pre_i];                                             \
            if (max_pre_end_sn < dp_end_sn[pre_i]) max_pre_end_sn = dp_end_sn[pre_i];                                             \
        } if (beg_sn < min_pre_beg_sn) beg_sn = min_pre_beg_sn;                                                                   \
        dp_beg_sn[dp_i] = beg_sn; beg = dp_beg[dp_i] = dp_beg_sn[dp_i] * pn;                                                      \
        end_sn = dp_end_sn[dp_i] = end/pn; end = dp_end[dp_i] = (dp_end_sn[dp_i]+1)*pn-1;                                         \
    }                                                                                                                             \
    /* loop query */                                                                                                              \
    /* first pre_node */                                                                                                          \
    pre_i = pre_index[dp_i][0];                                                                                                   \
    pre_dp_h = DP_HEF + pre_i * 3 * dp_sn; pre_dp_e1 = pre_dp_h + dp_sn;                                                          \
    pre_end = dp_end[pre_i]; pre_beg_sn = dp_beg_sn[pre_i]; pre_end_sn = dp_end_sn[pre_i];                                        \
    /* set M from (pre_i, q_i-1) */                                                                                               \
    if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                                   \
        _beg_sn = 0, _end_sn = end_sn; first = SIMDShiftRight(zero, SIMDTotalBytes-SIMDShiftOneN);                                \
    } else {                                                                                                                      \
        if (pre_beg_sn < beg_sn) _beg_sn = beg_sn, first = SIMDShiftRight(pre_dp_h[beg_sn-1], SIMDTotalBytes-SIMDShiftOneN);      \
        else _beg_sn = pre_beg_sn, first = SIMDShiftRight(SIMD_INF_MIN, SIMDTotalBytes-SIMDShiftOneN);                            \
        _end_sn = MIN_OF_THREE((pre_end+1)/pn, end_sn, dp_sn-1);                                                                  \
        for (i = beg_sn; i < _beg_sn; ++i) dp_h[i] = SIMD_INF_MIN;                                                                \
        for (i = _end_sn+1; i <= MIN_OF_TWO(end_sn+1, dp_sn-1); ++i) dp_h[i] = SIMD_INF_MIN;                                      \
    }                                                                                                                             \
    for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */                                                    \
        remain = SIMDShiftLeft(pre_dp_h[sn_i], SIMDShiftOneN);                                                                    \
        dp_h[sn_i] = SIMDOri(first, remain);                                                                                      \
        first = SIMDShiftRight(pre_dp_h[sn_i], SIMDTotalBytes-SIMDShiftOneN);                                                     \
    }                                                                                                                             \
    /* set E from (pre_i, q_i) */                                                                                                 \
    if (abpt->align_mode != ABPOA_LOCAL_MODE) {                                                                                   \
        _end_sn = MIN_OF_TWO(pre_end_sn, end_sn);                                                                                 \
        for (i = beg_sn; i < _beg_sn; ++i) dp_e1[i] = SIMD_INF_MIN;                                                               \
        for (i = _end_sn+1; i <= end_sn; ++i) dp_e1[i] = SIMD_INF_MIN;                                                            \
    }                                                                                                                             \
    for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i)   /* SIMD parallelization */                                                    \
        dp_e1[sn_i] = pre_dp_e1[sn_i];                                                                                            \
    /* get max m and e */                                                                                                         \
    for (i = 1; i < pre_n[dp_i]; ++i) {                                                                                           \
        pre_i = pre_index[dp_i][i];                                                                                               \
        pre_dp_h = DP_HEF + pre_i * 3 * dp_sn; pre_dp_e1 = pre_dp_h + dp_sn;                                                      \
        pre_end = dp_end[pre_i]; pre_beg_sn = dp_beg_sn[pre_i]; pre_end_sn = dp_end_sn[pre_i];                                    \
        /* set M from (pre_i, q_i-1) */                                                                                           \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                               \
            first = SIMDShiftRight(zero, SIMDTotalBytes-SIMDShiftOneN);                                                           \
        } else {                                                                                                                  \
            if (pre_beg_sn < beg_sn) _beg_sn = beg_sn, first = SIMDShiftRight(pre_dp_h[beg_sn-1], SIMDTotalBytes-SIMDShiftOneN);  \
            else _beg_sn = pre_beg_sn, first = SIMDShiftRight(SIMD_INF_MIN, SIMDTotalBytes-SIMDShiftOneN);                        \
            _end_sn = MIN_OF_THREE((pre_end+1)/pn, end_sn, dp_sn-1);                                                              \
        }                                                                                                                         \
        for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */                                                \
            remain = SIMDShiftLeft(pre_dp_h[sn_i], SIMDShiftOneN);                                                                \
            dp_h[sn_i] = SIMDMax(SIMDOri(first, remain), dp_h[sn_i]);                                                             \
            first = SIMDShiftRight(pre_dp_h[sn_i], SIMDTotalBytes-SIMDShiftOneN);                                                 \
        }                                                                                                                         \
        /* set E from (pre_i, q_i) */                                                                                             \
        _end_sn = MIN_OF_TWO(pre_end_sn, end_sn);                                                                                 \
        for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i)   /* SIMD parallelization */                                                \
            dp_e1[sn_i] = SIMDMax(pre_dp_e1[sn_i], dp_e1[sn_i]);                                                                  \
    }                                                                                                                             \
    /* compare M, E, and F */                                                                                                     \
    for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { /* SIMD parallelization */                                                      \
        dp_h[sn_i] = SIMDAdd(dp_h[sn_i], q[sn_i]);                                                                                \
    }                                                                                                                             \
    /* new F start */                                                                                                             \
    first = SIMDShiftRight(SIMDShiftLeft(dp_h[beg_sn], SIMDTotalBytes-SIMDShiftOneN), SIMDTotalBytes-SIMDShiftOneN);              \
    for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) {                                                                                 \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                               \
            set_num  = pn;                                                                                                        \
        } else {                                                                                                                  \
            if (sn_i < min_pre_beg_sn) {                                                                                          \
                _err_fatal_simple(__func__, "sn_i < min_pre_beg_sn\n");                                                           \
            } else if (sn_i > max_pre_end_sn) {                                                                                   \
                set_num = sn_i == max_pre_end_sn+1 ? 2 : 1;                                                                       \
            } else set_num = pn;                                                                                                  \
        }                                                                                                                         \
        /* F = (H << 1 | x) - OE */                                                                                               \
        dp_f1[sn_i] = SIMDSub(SIMDOri(SIMDShiftLeft(dp_h[sn_i], SIMDShiftOneN), first), GAP_OE1);                                 \
        /* F = max{F, (F-e)<<1}, F = max{F, (F-2e)<<2} ... */                                                                     \
        SIMD_SET_F(dp_f1[sn_i], log_n, set_num, PRE_MIN, PRE_MASK, SUF_MIN, GAP_E1S, SIMDMax, SIMDAdd, SIMDSub, SIMDShiftOneN);   \
        /* x = max{H, F+o} */                                                                                                     \
        first = SIMDShiftRight(SIMDMax(dp_h[sn_i], SIMDAdd(dp_f1[sn_i], GAP_O1)), SIMDTotalBytes-SIMDShiftOneN);                  \
        /* H = max{H, F} */                                                                                                       \
        dp_h[sn_i] = SIMDMax(dp_h[sn_i], dp_e1[sn_i]); SIMDi tmp = dp_h[sn_i];                                                    \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                               \
            dp_h[sn_i] = SIMDMax(zero, SIMDMax(dp_h[sn_i], dp_f1[sn_i]));                                                         \
            SIMDSetIfEqual(dp_e1[sn_i], dp_h[sn_i],tmp, SIMDMax(SIMDSub(dp_e1[sn_i],GAP_E1), SIMDSub(dp_h[sn_i],GAP_OE1)),zero);  \
        } else {                                                                                                                  \
            dp_h[sn_i] = SIMDMax(dp_h[sn_i], dp_f1[sn_i]);                                                                        \
            SIMDSetIfEqual(dp_e1[sn_i], dp_h[sn_i],tmp, SIMDMax(SIMDSub(dp_e1[sn_i],GAP_E1), SIMDSub(dp_h[sn_i],GAP_OE1)),SIMD_INF_MIN); \
        }                                                                                                                         \
    }                                                                                                                             \
}

#define simd_abpoa_cg_dp(score_t, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub, SIMDGetIfGreater, SIMDSetIfGreater, SIMDSetIfEqual) { \
    node_id = abpoa_graph_index_to_node_id(graph, index_i);                                                                       \
    SIMDi *q = qp + graph->node[node_id].base * qp_sn, first, remain;                                                             \
    dp_h = DP_H2E2F+dp_i*5*dp_sn; dp_e1 = dp_h+dp_sn; dp_e2 = dp_e1+dp_sn; dp_f1 = dp_e2+dp_sn; dp_f2 = dp_f1+dp_sn;              \
    _dp_h=(score_t*)dp_h, _dp_e1=(score_t*)dp_e1, _dp_e2=(score_t*)dp_e2, _dp_f1=(score_t*)dp_f1, _dp_f2=(score_t*)dp_f2;         \
    int min_pre_beg_sn, max_pre_end_sn;                                                                                           \
    if (abpt->wb < 0) {                                                                                                           \
        beg = dp_beg[dp_i] = 0, end = dp_end[dp_i] = qlen;                                                                        \
        beg_sn = dp_beg_sn[dp_i] = beg/pn; end_sn = dp_end_sn[dp_i] = end/pn;                                                     \
        min_pre_beg_sn = 0, max_pre_end_sn = end_sn;                                                                              \
    } else {                                                                                                                      \
        beg = GET_AD_DP_BEGIN(graph, w, node_id, end_node_id, qlen), end = GET_AD_DP_END(graph, w, node_id, end_node_id, qlen);   \
        beg_sn = beg / pn; min_pre_beg_sn = INT32_MAX, max_pre_end_sn = -1;                                                       \
        for (i = 0; i < pre_n[dp_i]; ++i) {                                                                                       \
            pre_i = pre_index[dp_i][i];                                                                                           \
            if (min_pre_beg_sn > dp_beg_sn[pre_i]) min_pre_beg_sn = dp_beg_sn[pre_i];                                             \
            if (max_pre_end_sn < dp_end_sn[pre_i]) max_pre_end_sn = dp_end_sn[pre_i];                                             \
        } if (beg_sn < min_pre_beg_sn) beg_sn = min_pre_beg_sn;                                                                   \
        dp_beg_sn[dp_i] = beg_sn; beg = dp_beg[dp_i] = dp_beg_sn[dp_i] * pn;                                                      \
        end_sn = dp_end_sn[dp_i] = end/pn; end = dp_end[dp_i] = (dp_end_sn[dp_i]+1)*pn-1;                                         \
     /* fprintf(stderr, "index: %d, beg: %d, end: %d, beg_sn: %d, end_sn: %d\n", index_i, beg, end, beg_sn, end_sn); */           \
    }                                                                                                                             \
 /* fprintf(stderr, "%d: beg, end: %d, %d\n", index_i, beg, end); */                                                              \
    /* tot_dp_sn += (end_sn - beg_sn + 1); */                                                                                     \
    /* loop query */                                                                                                              \
    /* first pre_node */                                                                                                          \
    pre_i = pre_index[dp_i][0];                                                                                                   \
    pre_dp_h = DP_H2E2F + pre_i * 5 * dp_sn; pre_dp_e1 = pre_dp_h + dp_sn; pre_dp_e2 = pre_dp_e1 + dp_sn;                         \
    pre_end = dp_end[pre_i]; pre_beg_sn = dp_beg_sn[pre_i]; pre_end_sn = dp_end_sn[pre_i];                                        \
    /* set M from (pre_i, q_i-1) */                                                                                               \
    if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                                   \
        _beg_sn = 0, _end_sn = end_sn; first = SIMDShiftRight(zero, SIMDTotalBytes-SIMDShiftOneN);                                \
    } else {                                                                                                                      \
        if (pre_beg_sn < beg_sn) _beg_sn = beg_sn, first = SIMDShiftRight(pre_dp_h[beg_sn-1], SIMDTotalBytes-SIMDShiftOneN);      \
        else _beg_sn = pre_beg_sn, first = SIMDShiftRight(SIMD_INF_MIN, SIMDTotalBytes-SIMDShiftOneN);                            \
        _end_sn = MIN_OF_THREE((pre_end+1)/pn, end_sn, dp_sn-1);                                                                  \
        for (i = beg_sn; i < _beg_sn; ++i) dp_h[i] = SIMD_INF_MIN;                                                                \
        for (i = _end_sn+1; i <= MIN_OF_TWO(end_sn+1, dp_sn-1); ++i) dp_h[i] = SIMD_INF_MIN;                                      \
    }                                                                                                                             \
 /* fprintf(stderr, "1 index_i: %d, beg_sn: %d, end_sn: %d\n", index_i, _beg_sn, _end_sn); */                                     \
    for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */                                                    \
        remain = SIMDShiftLeft(pre_dp_h[sn_i], SIMDShiftOneN);                                                                    \
        dp_h[sn_i] = SIMDOri(first, remain);                                                                                      \
        first = SIMDShiftRight(pre_dp_h[sn_i], SIMDTotalBytes-SIMDShiftOneN);                                                     \
    }                                                                                                                             \
    /* set E from (pre_i, q_i) */                                                                                                 \
    if (abpt->align_mode != ABPOA_LOCAL_MODE) {                                                                                   \
        _end_sn = MIN_OF_TWO(pre_end_sn, end_sn);                                                                                 \
        for (i = beg_sn; i < _beg_sn; ++i) dp_e1[i] = SIMD_INF_MIN, dp_e2[i] = SIMD_INF_MIN;                                      \
        for (i = _end_sn+1; i <= end_sn; ++i) dp_e1[i] = SIMD_INF_MIN, dp_e2[i] = SIMD_INF_MIN;                                   \
    }                                                                                                                             \
 /* fprintf(stderr, "2 index_i: %d, beg_sn: %d, end_sn: %d\n", index_i, _beg_sn, _end_sn); */                                     \
    for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */                                                    \
        dp_e1[sn_i] = pre_dp_e1[sn_i];                                                                                            \
        dp_e2[sn_i] = pre_dp_e2[sn_i];                                                                                            \
    }                                                                                                                             \
    /* get max m and e */                                                                                                         \
    for (i = 1; i < pre_n[dp_i]; ++i) {                                                                                           \
        pre_i = pre_index[dp_i][i];                                                                                               \
        pre_dp_h = DP_H2E2F + (pre_i * 5) * dp_sn; pre_dp_e1 = pre_dp_h + dp_sn; pre_dp_e2 = pre_dp_e1 + dp_sn;                   \
        pre_end = dp_end[pre_i]; pre_beg_sn = dp_beg_sn[pre_i]; pre_end_sn = dp_end_sn[pre_i];                                    \
        /* set M from (pre_i, q_i-1) */                                                                                           \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                               \
            first = SIMDShiftRight(zero, SIMDTotalBytes-SIMDShiftOneN);                                                           \
        } else {                                                                                                                  \
            if (pre_beg_sn < beg_sn) _beg_sn = beg_sn, first = SIMDShiftRight(pre_dp_h[beg_sn-1], SIMDTotalBytes-SIMDShiftOneN);  \
            else _beg_sn = pre_beg_sn, first = SIMDShiftRight(SIMD_INF_MIN, SIMDTotalBytes-SIMDShiftOneN);                        \
            _end_sn = MIN_OF_THREE((pre_end+1)/pn, end_sn, dp_sn-1);                                                              \
        }                                                                                                                         \
     /* fprintf(stderr, "3 index_i: %d, beg_sn: %d, end_sn: %d\n", index_i, _beg_sn, _end_sn); */                                 \
        for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */                                                \
            remain = SIMDShiftLeft(pre_dp_h[sn_i], SIMDShiftOneN);                                                                \
            dp_h[sn_i] = SIMDMax(SIMDOri(first, remain), dp_h[sn_i]);                                                             \
            first = SIMDShiftRight(pre_dp_h[sn_i], SIMDTotalBytes-SIMDShiftOneN);                                                 \
        }                                                                                                                         \
        /* set E from (pre_i, q_i) */                                                                                             \
        _end_sn = MIN_OF_TWO(pre_end_sn, end_sn);                                                                                 \
     /* fprintf(stderr, "4 index_i: %d, beg_sn: %d, end_sn: %d\n", index_i, _beg_sn, _end_sn); */                                 \
        for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */                                                \
            dp_e1[sn_i] = SIMDMax(pre_dp_e1[sn_i], dp_e1[sn_i]);                                                                  \
            dp_e2[sn_i] = SIMDMax(pre_dp_e2[sn_i], dp_e2[sn_i]);                                                                  \
        }                                                                                                                         \
    }                                                                                                                             \
    /* compare M, E, and F */                                                                                                     \
 /* fprintf(stderr, "5 index_i: %d, beg_sn: %d, end_sn: %d\n", index_i, _beg_sn, _end_sn); */                                     \
    for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { /* SIMD parallelization */                                                      \
        dp_h[sn_i] = SIMDAdd(dp_h[sn_i], q[sn_i]);                                                                                \
    }                                                                                                                             \
    /* new F start */                                                                                                             \
    first = SIMDShiftRight(SIMDShiftLeft(dp_h[beg_sn], SIMDTotalBytes-SIMDShiftOneN), SIMDTotalBytes-SIMDShiftOneN);              \
    SIMDi first2 = first;                                                                                                         \
    for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) {                                                                                 \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) set_num = pn;                                                                   \
        else {                                                                                                                    \
            if (sn_i < min_pre_beg_sn) {                                                                                          \
                _err_fatal_simple(__func__, "sn_i < min_pre_beg_sn\n");                                                           \
            } else if (sn_i > max_pre_end_sn) {                                                                                   \
                set_num = sn_i == max_pre_end_sn+1 ? 2 : 1;                                                                       \
            } else set_num = pn;                                                                                                  \
        }                                                                                                                         \
        /* H = max{H, E} */                                                                                                       \
        dp_h[sn_i] =  SIMDMax(SIMDMax(dp_h[sn_i], dp_e1[sn_i]), dp_e2[sn_i]);                                                     \
        /* F = (H << 1 | x) - OE */                                                                                               \
        dp_f1[sn_i] = SIMDSub(SIMDOri(SIMDShiftLeft(dp_h[sn_i], SIMDShiftOneN), first), GAP_OE1);                                 \
        dp_f2[sn_i] = SIMDSub(SIMDOri(SIMDShiftLeft(dp_h[sn_i], SIMDShiftOneN), first2), GAP_OE2);                                \
        /* F = max{F, (F-e)<<1}, F = max{F, (F-2e)<<2} ... */                                                                     \
        SIMD_SET_F(dp_f1[sn_i], log_n, set_num, PRE_MIN, PRE_MASK, SUF_MIN, GAP_E1S, SIMDMax, SIMDAdd, SIMDSub, SIMDShiftOneN);   \
        SIMD_SET_F(dp_f2[sn_i], log_n, set_num, PRE_MIN, PRE_MASK, SUF_MIN, GAP_E2S, SIMDMax, SIMDAdd, SIMDSub, SIMDShiftOneN);   \
        /* x = max{H, F+o} */                                                                                                     \
        first = SIMDShiftRight(SIMDMax(dp_h[sn_i], SIMDAdd(dp_f1[sn_i], GAP_O1)), SIMDTotalBytes-SIMDShiftOneN);                  \
        first2 = SIMDShiftRight(SIMDMax(dp_h[sn_i], SIMDAdd(dp_f2[sn_i], GAP_O2)), SIMDTotalBytes-SIMDShiftOneN);                 \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                                               \
            dp_h[sn_i] = SIMDMax(zero, SIMDMax(dp_h[sn_i], SIMDMax(dp_f1[sn_i], dp_f2[sn_i])));                                   \
            dp_e1[sn_i] = SIMDMax(zero,SIMDMax(SIMDSub(dp_e1[sn_i],GAP_E1),SIMDSub(dp_h[sn_i],GAP_OE1)));                         \
            dp_e2[sn_i] = SIMDMax(zero,SIMDMax(SIMDSub(dp_e2[sn_i],GAP_E2),SIMDSub(dp_h[sn_i],GAP_OE2)));                         \
        } else {                                                                                                                  \
            /* H = max{H, F}    */                                                                                                \
            dp_h[sn_i] = SIMDMax(dp_h[sn_i], SIMDMax(dp_f1[sn_i], dp_f2[sn_i]));                                                  \
            /* e for next cell */                                                                                                 \
            dp_e1[sn_i] = SIMDMax(SIMDSub(dp_e1[sn_i],GAP_E1),SIMDSub(dp_h[sn_i],GAP_OE1));                                       \
            dp_e2[sn_i] = SIMDMax(SIMDSub(dp_e2[sn_i],GAP_E2),SIMDSub(dp_h[sn_i],GAP_OE2));                                       \
        }                                                                                                                         \
    }                                                                                                                             \
}

#define set_global_max_score(score, i, j) {         \
    if (score > best_score) {                       \
        best_score = score; best_i = i; best_j = j; \
    }                                               \
}

#define set_extend_max_score(score, i, j) {                                                              \
    if (score > best_score) {                                                                            \
        best_score = score; best_i = i; best_j = j; best_id = node_id;                                   \
    } else if (abpt->zdrop > 0) {                                                                        \
        int delta_index = graph->node_id_to_max_remain[best_id] - graph->node_id_to_max_remain[node_id]; \
        if (best_score - score > abpt->zdrop + gap_ext1 * abs(delta_index-(j-best_j)))                   \
            break;                                                                                       \
    }                                                                                                    \
}

#define simd_abpoa_global_get_max(score_t, DP_M, dp_sn) {      \
    int end, in_id, in_index, in_dp_i;                         \
    for (i = 0; i < graph->node[end_node_id].in_edge_n; ++i) { \
        in_id = graph->node[end_node_id].in_id[i];             \
        in_index = abpoa_graph_node_id_to_index(graph, in_id); \
        if (index_map[in_index] == 0) continue;                \
        in_dp_i = in_index - beg_index;                        \
        dp_h = DP_M + in_dp_i * dp_sn;                         \
        _dp_h = (score_t*)dp_h;                                \
        if (qlen > dp_end[in_dp_i]) end = dp_end[in_dp_i];     \
        else end = qlen;                                       \
        set_global_max_score(_dp_h[end], in_dp_i, end);        \
    }                                                          \
}

#define simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater) { \
    /* select max dp_h */                                                    \
    max = inf_min, max_i = -1;                                               \
    SIMDi a = dp_h[end_sn], b = qi[end_sn];                                  \
    if (end_sn == qlen / pn) SIMDSetIfGreater(a, zero, b, SIMD_INF_MIN, a);  \
    for (i = beg_sn; i < end_sn; ++i) {                                      \
        SIMDGetIfGreater(b, a, dp_h[i], a, qi[i], b);                        \
    }                                                                        \
    _dp_h = (score_t*)&a, _qi = (score_t*)&b;                                \
    for (i = 0; i < pn; ++i) {                                               \
        if (_dp_h[i] > max) {                                                \
            max = _dp_h[i]; max_i = _qi[i];                                  \
        }                                                                    \
    }                                                                        \
}

#define simd_abpoa_ada_max_i   {                                                                                        \
    /* set max_pos_left/right for next nodes */                                                                         \
    int out_i = max_i + 1;                                                                                              \
    for (i = 0; i < graph->node[node_id].out_edge_n; ++i) {                                                             \
        int out_node_id = graph->node[node_id].out_id[i];                                                               \
        if (out_i > graph->node_id_to_max_pos_right[out_node_id]) graph->node_id_to_max_pos_right[out_node_id] = out_i; \
        if (out_i < graph->node_id_to_max_pos_left[out_node_id]) graph->node_id_to_max_pos_left[out_node_id] = out_i;   \
    }                                                                                                                   \
}

// TODO end_bonus for extension
// linear gap penalty: gap_open1 == 0
#define simd_abpoa_lg_align_sequence_to_graph_core(score_t, sp, SIMDSetOne, SIMDMax, SIMDAdd,   \
        SIMDSub, SIMDShiftOneN, SIMDSetIfGreater, SIMDGetIfGreater) {                           \
    simd_abpoa_lg_var(score_t, sp, SIMDSetOne, SIMDShiftOneN, SIMDAdd);                         \
    simd_abpoa_lg_first_dp(score_t);                                                            \
    for (index_i=beg_index+1, dp_i=1; index_i<end_index; ++index_i, ++dp_i) {                   \
        if (index_map[index_i] == 0) continue;                                                  \
        simd_abpoa_lg_dp(score_t, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub);                    \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                             \
            simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater);                 \
            set_global_max_score(max, dp_i, max_i);                                             \
        }                                                                                       \
        if (abpt->align_mode == ABPOA_EXTEND_MODE) {                                            \
            simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater);                 \
            set_extend_max_score(max, dp_i, max_i);                                             \
        }                                                                                       \
        if (abpt->wb >= 0) {                                                                    \
            if (abpt->align_mode == ABPOA_GLOBAL_MODE) {                                        \
                simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater);             \
            }                                                                                   \
            simd_abpoa_ada_max_i;                                                               \
        }                                                                                       \
    }                                                                                           \
    if (abpt->align_mode == ABPOA_GLOBAL_MODE) simd_abpoa_global_get_max(score_t, DP_H, dp_sn); \
    res->best_score = best_score;                                                               \
 /* simd_abpoa_print_lg_matrix(score_t, beg_index, end_index); printf("best_score: (%d, %d) -> %d\n", best_i, best_j, best_score); */                          \
    if (abpt->ret_cigar) simd_abpoa_lg_backtrack(score_t);                                      \
    simd_abpoa_free_var; SIMDFree(GAP_E1S);                                                     \
}

// affine gap penalty: gap_open1 > 0
#define simd_abpoa_ag_align_sequence_to_graph_core(score_t, sp, SIMDSetOne, SIMDMax, SIMDAdd,       \
        SIMDSub, SIMDShiftOneN, SIMDSetIfGreater, SIMDGetIfGreater, SIMDSetIfEqual) {               \
    simd_abpoa_ag_var(score_t, sp, SIMDSetOne, SIMDShiftOneN, SIMDAdd);                             \
    simd_abpoa_ag_first_dp(score_t);                                                                \
    for (index_i=beg_index+1, dp_i=1; index_i<end_index; ++index_i, ++dp_i) {                       \
        if (index_map[index_i] == 0) continue;                                                      \
        simd_abpoa_ag_dp(score_t, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub, SIMDGetIfGreater, SIMDSetIfGreater, SIMDSetIfEqual); \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                 \
            simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater);                     \
            set_global_max_score(max, dp_i, max_i);                                                 \
        } else if (abpt->align_mode == ABPOA_EXTEND_MODE) {                                         \
            simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater);                     \
            set_extend_max_score(max, dp_i, max_i);                                                 \
        }                                                                                           \
        if (abpt->wb >= 0) {                                                                        \
            if (abpt->align_mode == ABPOA_GLOBAL_MODE) {                                            \
                simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater);                 \
            }                                                                                       \
            simd_abpoa_ada_max_i;                                                                   \
        }                                                                                           \
    }                                                                                               \
    if (abpt->align_mode == ABPOA_GLOBAL_MODE) simd_abpoa_global_get_max(score_t, DP_HEF, 3*dp_sn); \
    res->best_score = best_score;                                                                   \
 /* simd_abpoa_print_ag_matrix(score_t, beg_index, end_index); fprintf(stderr, "best_score: (%d, %d) -> %d\n", best_i, best_j, best_score); */ \
    if (abpt->ret_cigar) simd_abpoa_ag_backtrack(score_t);                                          \
    simd_abpoa_free_var; SIMDFree(GAP_E1S);                                                         \
}

// convex gap penalty: gap_open1 > 0 && gap_open2 > 0
#define simd_abpoa_cg_align_sequence_to_graph_core(score_t, sp, SIMDSetOne, SIMDMax, SIMDAdd,         \
        SIMDSub, SIMDShiftOneN, SIMDSetIfGreater, SIMDGetIfGreater, SIMDSetIfEqual) {                 \
    simd_abpoa_cg_var(score_t, sp, SIMDSetOne, SIMDShiftOneN, SIMDAdd);                               \
    simd_abpoa_cg_first_dp(score_t);                                                                  \
    for (index_i=beg_index+1, dp_i=1; index_i<end_index; ++index_i, ++dp_i) {                         \
        if (index_map[index_i] == 0) continue;                                                        \
        simd_abpoa_cg_dp(score_t, SIMDShiftOneN, SIMDMax, SIMDAdd, SIMDSub, SIMDGetIfGreater, SIMDSetIfGreater, SIMDSetIfEqual); \
        if (abpt->align_mode == ABPOA_LOCAL_MODE) {                                                   \
            simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater);                       \
            set_global_max_score(max, dp_i, max_i);                                                   \
        } else if (abpt->align_mode == ABPOA_EXTEND_MODE) {                                           \
            simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater);                       \
            set_extend_max_score(max, dp_i, max_i);                                                   \
        }                                                                                             \
        if (abpt->wb >= 0) {                                                                          \
            if (abpt->align_mode == ABPOA_GLOBAL_MODE) {                                              \
                simd_abpoa_max_in_row(score_t, SIMDSetIfGreater, SIMDGetIfGreater);                   \
            }                                                                                         \
            simd_abpoa_ada_max_i;                                                                     \
        }                                                                                             \
    }                                                                                                 \
 /* printf("dp_sn: %d\n", tot_dp_sn); */                                                              \
    if (abpt->align_mode == ABPOA_GLOBAL_MODE) simd_abpoa_global_get_max(score_t, DP_H2E2F, 5*dp_sn); \
    res->best_score = best_score;                                                                     \
/* simd_abpoa_print_cg_matrix(score_t, beg_index, end_index); fprintf(stderr,"best_score: (%d, %d) -> %d\n",best_i,best_j,best_score); */ \
    if (abpt->ret_cigar) simd_abpoa_cg_backtrack(score_t);                                            \
    simd_abpoa_free_var; SIMDFree(GAP_E1S); SIMDFree(GAP_E2S);                                        \
}

abpoa_simd_matrix_t *abpoa_init_simd_matrix(void) {
    abpoa_simd_matrix_t *abm = (abpoa_simd_matrix_t*)_err_malloc(sizeof(abpoa_simd_matrix_t));
    abm->s_msize = 0; abm->s_mem = NULL; abm->rang_m = 0;
    abm->dp_beg = NULL; abm->dp_end = NULL; abm->dp_beg_sn = NULL; abm->dp_end_sn = NULL;
    return abm;
}

void abpoa_free_simd_matrix(abpoa_simd_matrix_t *abm) {
    if (abm->s_mem) SIMDFree(abm->s_mem);
    if (abm->dp_beg) {
        free(abm->dp_beg); free(abm->dp_end); free(abm->dp_beg_sn); free(abm->dp_end_sn);
    } free(abm);
}

// realloc memory everytime the graph is updated (nodes are updated already)
// * index_to_node_id/node_id_to_index/node_id_to_max_remain, max_pos_left/right
// * qp, DP_HE/H (if ag/lg), dp_f, qi (if ada/extend)
// * dp_beg/end, dp_beg/end_sn if band
// * pre_n, pre_index
int simd_abpoa_realloc(abpoa_t *ab, int gn, int qlen, abpoa_para_t *abpt, SIMD_para_t sp) {
    uint64_t pn = sp.num_of_value, size = sp.size, sn = (qlen + sp.num_of_value) / pn;
    uint64_t s_msize = sn * abpt->m * size; // qp

    if (abpt->gap_mode == ABPOA_LINEAR_GAP) s_msize += (sn * gn * size); // DP_H, linear
    else if (abpt->gap_mode == ABPOA_AFFINE_GAP) s_msize += (sn * gn * 3 * size); // DP_HEF, affine
    else s_msize += (sn * gn * 5 * size); // DP_H2E2F, convex

    if (abpt->wb >= 0 || abpt->align_mode == ABPOA_EXTEND_MODE) // qi
        s_msize += sn * size;

    // if (s_msize > UINT32_MAX) {
        // err_func_format_printf(__func__, "Warning: Graph is too large or query is too long.\n");
        // return 1;
    // }
    // fprintf(stderr, "%lld, %lld, %lld\n", (long long)node_n, (long long)ab->abm->s_msize, (long long)s_msize);
    if (s_msize > ab->abm->s_msize) {
        if (ab->abm->s_mem) SIMDFree(ab->abm->s_mem);
        kroundup64(s_msize); ab->abm->s_msize = s_msize;
        ab->abm->s_mem = (SIMDi*)SIMDMalloc(ab->abm->s_msize, size);
    }

    if (gn > ab->abm->rang_m) {
        ab->abm->rang_m = gn; kroundup32(ab->abm->rang_m);
        ab->abm->dp_beg = (int*)_err_realloc(ab->abm->dp_beg, ab->abm->rang_m * sizeof(int));
        ab->abm->dp_end = (int*)_err_realloc(ab->abm->dp_end, ab->abm->rang_m * sizeof(int));
        ab->abm->dp_beg_sn = (int*)_err_realloc(ab->abm->dp_beg_sn, ab->abm->rang_m * sizeof(int));
        ab->abm->dp_end_sn = (int*)_err_realloc(ab->abm->dp_end_sn, ab->abm->rang_m * sizeof(int));
    }
    return 0;
}

void abpoa_init_var(abpoa_para_t *abpt, uint8_t *query, int qlen, SIMDi *qp, SIMDi *qi, int *mat, int qp_sn, int pn, SIMDi SIMD_INF_MIN) {
    int i, j, k; int32_t *_qi;
    /* generate the query profile */
    for (i = 0; i < qp_sn * abpt->m; ++i) qp[i] = SIMD_INF_MIN;
    for (k = 0; k < abpt->m; ++k) { /* SIMD parallelization */
        int *p = &mat[k * abpt->m];
        int32_t *_qp = (int32_t*)(qp + k * qp_sn); _qp[0] = 0;
        for (j = 0; j < qlen; ++j) _qp[j+1] = (int32_t)p[query[j]];
        for (j = qlen+1; j < qp_sn * pn; ++j) _qp[j] = 0;
    }
    if (abpt->wb >= 0 || abpt->align_mode == ABPOA_EXTEND_MODE) { /* query index */
        _qi = (int32_t*)qi;
        for (i = 0; i <= qlen; ++i) _qi[i] = i;
        for (i = qlen+1; i < (qlen/pn+1) * pn; ++i) _qi[i] = -1;
    }
}

void abpoa_cg_first_dp(abpoa_para_t *abpt, abpoa_graph_t *graph, uint8_t *index_map, int beg_node_id, int end_node_id, int *dp_beg, int *dp_end, int *dp_beg_sn, int *dp_end_sn, int pn, int qlen, int w, int dp_sn, SIMDi *DP_H2E2F, SIMDi SIMD_INF_MIN, int32_t inf_min, int gap_open1, int gap_ext1, int gap_open2, int gap_ext2, int gap_oe1, int gap_oe2) {
    int i, _end_sn;
    if (abpt->wb >= 0) {
        graph->node_id_to_max_pos_left[beg_node_id] = graph->node_id_to_max_pos_right[beg_node_id] = 0;
        for (i = 0; i < graph->node[beg_node_id].out_edge_n; ++i) { /* set min/max rank for next_id */
            int out_id = graph->node[beg_node_id].out_id[i];
            if (index_map[abpoa_graph_node_id_to_index(graph, out_id)])
                graph->node_id_to_max_pos_left[out_id] = graph->node_id_to_max_pos_right[out_id] = 1;
        }
        dp_beg[0] = GET_AD_DP_BEGIN(graph, w, beg_node_id, end_node_id, qlen), dp_end[0] = GET_AD_DP_END(graph, w, beg_node_id, end_node_id, qlen);
    } else {
        dp_beg[0] = 0, dp_end[0] = qlen;
    }
    dp_beg_sn[0] = (dp_beg[0])/pn; dp_end_sn[0] = (dp_end[0])/pn;
    dp_beg[0] = dp_beg_sn[0] * pn; dp_end[0] = (dp_end_sn[0]+1)*pn-1;
    SIMDi *dp_h = DP_H2E2F; SIMDi *dp_e1 = dp_h + dp_sn; SIMDi *dp_e2 = dp_e1 + dp_sn, *dp_f1 = dp_e2 + dp_sn, *dp_f2 = dp_f1 + dp_sn;
    _end_sn = MIN_OF_TWO(dp_end_sn[0]+1, dp_sn-1);

    for (i = 0; i <= _end_sn; ++i) {
        dp_h[i] = SIMD_INF_MIN; dp_e1[i] = SIMD_INF_MIN; dp_e2[i] = SIMD_INF_MIN;
    }
    int32_t *_dp_h = (int32_t*)dp_h, *_dp_e1 = (int32_t*)dp_e1, *_dp_e2 = (int32_t*)dp_e2, *_dp_f1 = (int32_t*)dp_f1, *_dp_f2 = (int32_t*)dp_f2;
    _dp_h[0] = 0; _dp_e1[0] = -(gap_oe1); _dp_e2[0] = -(gap_oe2); _dp_f1[0] = _dp_f2[0] = inf_min;
    for (i = 1; i <= dp_end[0]; ++i) { /* no SIMD parallelization */
        _dp_f1[i] = -(gap_open1 + gap_ext1 * i);
        _dp_f2[i] = -(gap_open2 + gap_ext2 * i);
        _dp_h[i] = MAX_OF_TWO(_dp_f1[i], _dp_f2[i]); // -MIN_OF_TWO(gap_open1+gap_ext1*i, gap_open2+gap_ext2*i);
    }
}

int abpoa_max(SIMDi SIMD_INF_MIN, SIMDi zero, int inf_min, SIMDi *dp_h, SIMDi *qi, int qlen, int pn, int beg_sn, int end_sn) {
    /* select max dp_h */
    int max = inf_min, max_i = -1, i;
    SIMDi a = dp_h[end_sn], b = qi[end_sn];
    if (end_sn == qlen / pn) SIMDSetIfGreateri32(a, zero, b, SIMD_INF_MIN, a);
    for (i = beg_sn; i < end_sn; ++i) {
        SIMDGetIfGreateri32(b, a, dp_h[i], a, qi[i], b);
    }
    int32_t *_dp_h = (int32_t*)&a, *_qi = (int32_t*)&b;
    for (i = 0; i < pn; ++i) {
        if (_dp_h[i] > max) {
            max = _dp_h[i]; max_i = _qi[i];
        }
    }
    return max_i;
}

void abpoa_ada_max_i(int max_i, abpoa_graph_t *graph, int node_id) {
    /* set max_pos_left/right for next nodes */
    int out_i = max_i + 1; int i;
    for (i = 0; i < graph->node[node_id].out_edge_n; ++i) {
        int out_node_id = graph->node[node_id].out_id[i];
        if (out_i > graph->node_id_to_max_pos_right[out_node_id]) graph->node_id_to_max_pos_right[out_node_id] = out_i;
        if (out_i < graph->node_id_to_max_pos_left[out_node_id]) graph->node_id_to_max_pos_left[out_node_id] = out_i;
    }
}

void abpoa_global_get_max(abpoa_graph_t *graph, int beg_index, int end_node_id, uint8_t *index_map, SIMDi *DP_H_HE, int dp_sn, int qlen, int *dp_end, int32_t *best_score, int *best_i, int *best_j) {
    int in_id, in_index, dp_i, i;
    for (i = 0; i < graph->node[end_node_id].in_edge_n; ++i) {
        in_id = graph->node[end_node_id].in_id[i];
        in_index = abpoa_graph_node_id_to_index(graph, in_id);
        if (index_map[in_index] == 0) continue;
        dp_i = in_index - beg_index;
        SIMDi *dp_h = DP_H_HE + dp_i * dp_sn;
        int32_t *_dp_h = (int32_t*)dp_h;
        int end;
        if (qlen > dp_end[dp_i]) end = dp_end[dp_i];
        else end = qlen;
        if (_dp_h[end] > *best_score) {
            *best_score = _dp_h[end]; *best_i = dp_i; *best_j = end;
        }
    }
}

int abpoa_cg_dp(SIMDi *q, SIMDi *dp_h, SIMDi *dp_e1, SIMDi *dp_e2, SIMDi *dp_f1, SIMDi *dp_f2, int **pre_index, int *pre_n, int index_i, int dp_i, abpoa_graph_t *graph, abpoa_para_t *abpt, int dp_sn, int pn, int qlen, int w, SIMDi *DP_H2E2F, SIMDi SIMD_INF_MIN, SIMDi GAP_O1, SIMDi GAP_O2, SIMDi GAP_E1, SIMDi GAP_E2, SIMDi GAP_OE1, SIMDi GAP_OE2, SIMDi* GAP_E1S, SIMDi* GAP_E2S, SIMDi *PRE_MIN, SIMDi *PRE_MASK, SIMDi *SUF_MIN, int log_n, int *dp_beg, int *dp_end, int *dp_beg_sn, int *dp_end_sn, int end_node_id) {
    int tot_dp_sn = 0, i, pre_i, node_id = abpoa_graph_index_to_node_id(graph, index_i);
    int min_pre_beg_sn, max_pre_end_sn, beg, end, beg_sn, end_sn, pre_end, pre_end_sn, pre_beg_sn, sn_i;
    if (abpt->wb < 0) {
        beg = dp_beg[dp_i] = 0, end = dp_end[dp_i] = qlen;
        beg_sn = dp_beg_sn[dp_i] = beg/pn; end_sn = dp_end_sn[dp_i] = end/pn;
        min_pre_beg_sn = 0, max_pre_end_sn = end_sn;
    } else {
        beg = GET_AD_DP_BEGIN(graph, w, node_id, end_node_id, qlen), end = GET_AD_DP_END(graph, w, node_id, end_node_id, qlen);
        beg_sn = beg / pn; min_pre_beg_sn = INT32_MAX, max_pre_end_sn = -1;
        for (i = 0; i < pre_n[dp_i]; ++i) {
            pre_i = pre_index[dp_i][i];
            if (min_pre_beg_sn > dp_beg_sn[pre_i]) min_pre_beg_sn = dp_beg_sn[pre_i];
            if (max_pre_end_sn < dp_end_sn[pre_i]) max_pre_end_sn = dp_end_sn[pre_i];
        } if (beg_sn < min_pre_beg_sn) beg_sn = min_pre_beg_sn;
        dp_beg_sn[dp_i] = beg_sn; beg = dp_beg[dp_i] = dp_beg_sn[dp_i] * pn;
        end_sn = dp_end_sn[dp_i] = end/pn; end = dp_end[dp_i] = (dp_end_sn[dp_i]+1)*pn-1;
#ifdef __DEBUG__
    fprintf(stderr, "index: %d (node: %d): beg: %d, end: %d, beg_sn: %d, end_sn: %d\n", index_i, node_id, beg, end, beg_sn, end_sn);
#endif
    }
    tot_dp_sn += (end_sn - beg_sn + 1);
    /* loop query */
    // new init start
    int _beg_sn, _end_sn;
    // first pre_node
    pre_i = pre_index[dp_i][0];
    SIMDi *pre_dp_h = DP_H2E2F + (pre_i * 5) * dp_sn; SIMDi *pre_dp_e1 = pre_dp_h + dp_sn; SIMDi *pre_dp_e2 = pre_dp_e1 + dp_sn;
    pre_end = dp_end[pre_i]; pre_beg_sn = dp_beg_sn[pre_i]; pre_end_sn = dp_end_sn[pre_i];
    SIMDi first, remain;
    /* set M from (pre_i, q_i-1) */
    if (pre_beg_sn < beg_sn) _beg_sn = beg_sn, first = SIMDShiftRight(pre_dp_h[beg_sn-1], SIMDTotalBytes-SIMDShiftOneNi32);
    else _beg_sn = pre_beg_sn, first = SIMDShiftRight(SIMD_INF_MIN, SIMDTotalBytes-SIMDShiftOneNi32);
    _end_sn = MIN_OF_THREE((pre_end+1)/pn, end_sn, dp_sn-1);
    for (i = beg_sn; i < _beg_sn; ++i) dp_h[i] = SIMD_INF_MIN;
    for (i = _end_sn+1; i <= MIN_OF_TWO(end_sn+1, dp_sn-1); ++i) dp_h[i] = SIMD_INF_MIN;
    for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */
        remain = SIMDShiftLeft(pre_dp_h[sn_i], SIMDShiftOneNi32);
        dp_h[sn_i] = SIMDOri(first, remain);
        first = SIMDShiftRight(pre_dp_h[sn_i], SIMDTotalBytes-SIMDShiftOneNi32);
    }
    /* set E from (pre_i, q_i) */
    _end_sn = MIN_OF_TWO(pre_end_sn, end_sn);
    for (i = beg_sn; i < _beg_sn; ++i) dp_e1[i] = SIMD_INF_MIN, dp_e2[i] = SIMD_INF_MIN;
    for (i = _end_sn+1; i <= end_sn; ++i) dp_e1[i] = SIMD_INF_MIN, dp_e2[i] = SIMD_INF_MIN;
    for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */
        dp_e1[sn_i] = pre_dp_e1[sn_i];
        dp_e2[sn_i] = pre_dp_e2[sn_i];
    }
    // if (index_i == 13095) debug_simd_abpoa_print_cg_matrix_row("1", int32_t, index_i);
    // new init end
    /* get max m and e */
    for (i = 1; i < pre_n[dp_i]; ++i) {
        pre_i = pre_index[dp_i][i];
        pre_dp_h = DP_H2E2F + (pre_i * 5) * dp_sn; pre_dp_e1 = pre_dp_h + dp_sn; pre_dp_e2 = pre_dp_e1 + dp_sn;
        pre_end = dp_end[pre_i]; pre_beg_sn = dp_beg_sn[pre_i]; pre_end_sn = dp_end_sn[pre_i];
        /* set M from (pre_i, q_i-1) */
        if (pre_beg_sn < beg_sn) _beg_sn = beg_sn, first = SIMDShiftRight(pre_dp_h[beg_sn-1], SIMDTotalBytes-SIMDShiftOneNi32);
        else _beg_sn = pre_beg_sn, first = SIMDShiftRight(SIMD_INF_MIN, SIMDTotalBytes-SIMDShiftOneNi32);
        _end_sn = MIN_OF_THREE((pre_end+1)/pn, end_sn, dp_sn-1);
        for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */
            remain = SIMDShiftLeft(pre_dp_h[sn_i], SIMDShiftOneNi32);
            dp_h[sn_i] = SIMDMaxi32(SIMDOri(first, remain), dp_h[sn_i]);
            first = SIMDShiftRight(pre_dp_h[sn_i], SIMDTotalBytes-SIMDShiftOneNi32);
        }
        /* set E from (pre_i, q_i) */
        _end_sn = MIN_OF_TWO(pre_end_sn, end_sn);
        for (sn_i = _beg_sn; sn_i <= _end_sn; ++sn_i) { /* SIMD parallelization */
            dp_e1[sn_i] = SIMDMaxi32(pre_dp_e1[sn_i], dp_e1[sn_i]);
            dp_e2[sn_i] = SIMDMaxi32(pre_dp_e2[sn_i], dp_e2[sn_i]);
        }
    }
    // debug_simd_abpoa_print_cg_matrix_row("2", int32_t, index_i);
    /* compare M, E, and F */
    for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) { /* SIMD parallelization */
        dp_h[sn_i] =  SIMDAddi32(dp_h[sn_i], q[sn_i]);
    }
    // debug_simd_abpoa_print_cg_matrix_row("3", int32_t, index_i);
    /* new F start */
    first = SIMDShiftRight(SIMDShiftLeft(dp_h[beg_sn], SIMDTotalBytes-SIMDShiftOneNi32), SIMDTotalBytes-SIMDShiftOneNi32);
    int set_num; SIMDi first2 = first;//, tmp;
    for (sn_i = beg_sn; sn_i <= end_sn; ++sn_i) {
        if (sn_i < min_pre_beg_sn) {
            _err_fatal_simple(__func__, "sn_i < min_pre_beg_sn\n");
        } else if (sn_i > max_pre_end_sn) {
            set_num = sn_i == max_pre_end_sn+1 ? 2 : 1;
        } else set_num = pn;
        /* H = max{H, E} */
        dp_h[sn_i] = SIMDMaxi32(SIMDMaxi32(dp_h[sn_i], dp_e1[sn_i]), dp_e2[sn_i]); // tmp = dp_h[sn_i];
        /* F = (H << 1 | x) - OE */
        // if (sn_i==beg_sn) debug_simd_abpoa_print_cg_matrix_row("4.1", int32_t, index_i);
        dp_f1[sn_i] = SIMDSubi32(SIMDOri(SIMDShiftLeft(dp_h[sn_i], SIMDShiftOneNi32), first), GAP_OE1);
        dp_f2[sn_i] = SIMDSubi32(SIMDOri(SIMDShiftLeft(dp_h[sn_i], SIMDShiftOneNi32), first2), GAP_OE2);
        /* F = max{F, (F-e)<<1}, F = max{F, (F-2e)<<2} ... */
        // if (sn_i==beg_sn) debug_simd_abpoa_print_cg_matrix_row("4.2", int32_t, index_i);
        SIMD_SET_F(dp_f1[sn_i], log_n, set_num, PRE_MIN, PRE_MASK, SUF_MIN, GAP_E1S, SIMDMaxi32, SIMDAddi32, SIMDSubi32, SIMDShiftOneNi32);
        SIMD_SET_F(dp_f2[sn_i], log_n, set_num, PRE_MIN, PRE_MASK, SUF_MIN, GAP_E2S, SIMDMaxi32, SIMDAddi32, SIMDSubi32, SIMDShiftOneNi32);
        /* x = max{H, F+o} */
        // if (sn_i==beg_sn) debug_simd_abpoa_print_cg_matrix_row("4.3", int32_t, index_i);
        first = SIMDShiftRight(SIMDMaxi32(dp_h[sn_i], SIMDAddi32(dp_f1[sn_i], GAP_O1)), SIMDTotalBytes-SIMDShiftOneNi32);
        first2 = SIMDShiftRight(SIMDMaxi32(dp_h[sn_i], SIMDAddi32(dp_f2[sn_i], GAP_O2)), SIMDTotalBytes-SIMDShiftOneNi32);
        /* H = max{H, F}    */
        dp_h[sn_i] = SIMDMaxi32(SIMDMaxi32(dp_h[sn_i], dp_f1[sn_i]), dp_f2[sn_i]);
        // if (sn_i==beg_sn) debug_simd_abpoa_print_cg_matrix_row("4.4", int32_t, index_i);
        /* e for next cell */
        // SIMDSetIfEquali32(dp_e1[sn_i], dp_h[sn_i], tmp, SIMDMaxi32(SIMDSubi32(dp_e1[sn_i], GAP_E1), SIMDSubi32(dp_h[sn_i], GAP_OE1)), SIMD_INF_MIN);
        dp_e1[sn_i] = SIMDMaxi32(SIMDSubi32(dp_e1[sn_i], GAP_E1), SIMDSubi32(dp_h[sn_i], GAP_OE1));
        // SIMDSetIfEquali32(dp_e2[sn_i], dp_h[sn_i], tmp, SIMDMaxi32(SIMDSubi32(dp_e2[sn_i], GAP_E2), SIMDSubi32(dp_h[sn_i], GAP_OE2)), SIMD_INF_MIN);
        dp_e2[sn_i] = SIMDMaxi32(SIMDSubi32(dp_e2[sn_i], GAP_E2), SIMDSubi32(dp_h[sn_i], GAP_OE2));
    }
    return tot_dp_sn;
}

void abpoa_cg_backtrack(SIMDi *DP_H2E2F, int **pre_index, int *pre_n, int *dp_beg, int *dp_end, int dp_sn, int m, int *mat, int gap_ext1, int gap_ext2, int gap_oe1, int gap_oe2, int beg_index, int best_dp_i, int best_dp_j, int qlen, abpoa_graph_t *graph, abpoa_para_t *abpt, uint8_t *query, abpoa_res_t *res) {
    int dp_i, dp_j, k, pre_i, n_c = 0, m_c = 0, id, hit, cur_op = ABPOA_ALL_OP, _start_i, _start_j;
    SIMDi *dp_h;
    int32_t s, is_match, *_dp_h=NULL, *_dp_e1, *_dp_e2, *_pre_dp_h, *_pre_dp_e1, *_pre_dp_e2, *_dp_f1, *_dp_f2; abpoa_cigar_t *cigar = 0;
    dp_i = best_dp_i, dp_j = best_dp_j, _start_i = best_dp_i, _start_j = best_dp_j;
    id = abpoa_graph_index_to_node_id(graph, dp_i+beg_index);
    if (best_dp_j < qlen) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, qlen-best_dp_j, -1, qlen-1);
    dp_h = DP_H2E2F + dp_sn * (dp_i * 5); _dp_h = (int32_t*)dp_h;
    int indel_first = 1; /* prefer to keep gaps at the end */
    while (dp_i > 0 && dp_j > 0) {
        _start_i = dp_i, _start_j = dp_j;
        int *pre_index_i = pre_index[dp_i];
        s = mat[m * graph->node[id].base + query[dp_j-1]]; hit = 0;
        is_match = graph->node[id].base == query[dp_j-1];
        if ((cur_op & ABPOA_M_OP) && (indel_first == 0)) {
            for (k = 0; k < pre_n[dp_i]; ++k) {
                pre_i = pre_index_i[k];
                if (dp_j-1 < dp_beg[pre_i] || dp_j-1 > dp_end[pre_i]) continue;
                _pre_dp_h = (int32_t*)(DP_H2E2F + dp_sn * (pre_i * 5));
                if (_pre_dp_h[dp_j-1] + s == _dp_h[dp_j]) { /* match/mismatch */
                    cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, dp_j-1);
                    dp_i = pre_i; --dp_j; id = abpoa_graph_index_to_node_id(graph, dp_i+beg_index); hit = 1;
                    dp_h = DP_H2E2F + dp_sn * (dp_i * 5); _dp_h = (int32_t*)dp_h;
                    cur_op = ABPOA_ALL_OP;
                    ++res->n_aln_bases; res->n_matched_bases += is_match ? 1 : 0;
                    break;
                }
            }
        }
        if (hit == 0 && cur_op & ABPOA_E_OP) {
            _dp_e1 = (int32_t*)(dp_h+dp_sn), _dp_e2 = (int32_t*)(dp_h+dp_sn*2);
            for (k = 0; k < pre_n[dp_i]; ++k) {
                pre_i = pre_index_i[k];
                if (dp_j < dp_beg[pre_i] || dp_j > dp_end[pre_i]) continue;
                _pre_dp_h = (int32_t*)(DP_H2E2F + dp_sn * (pre_i * 5));
                if (cur_op & ABPOA_E1_OP) {
                    _pre_dp_e1 = (int32_t*)(DP_H2E2F + dp_sn * ((pre_i * 5) + 1));
                    if (cur_op & ABPOA_M_OP) {
                        if (_dp_h[dp_j] == _pre_dp_e1[dp_j]) {
                            if (_pre_dp_h[dp_j] - gap_oe1 == _pre_dp_e1[dp_j]) cur_op = ABPOA_M_OP | ABPOA_F_OP;
                            else cur_op = ABPOA_E1_OP;
                            cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, dp_j-1);
                            dp_i = pre_i; id = abpoa_graph_index_to_node_id(graph, dp_i+beg_index); hit = 1;
                            dp_h = DP_H2E2F + dp_sn * (dp_i * 5); _dp_h = (int32_t*)dp_h;
                            break;
                        }
                    } else {
                        if (_dp_e1[dp_j] == _pre_dp_e1[dp_j] - gap_ext1) {
                            if (_pre_dp_h[dp_j] - gap_oe1 == _pre_dp_e1[dp_j]) cur_op = ABPOA_M_OP | ABPOA_F_OP;
                            else cur_op = ABPOA_E1_OP;
                            cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, dp_j-1);
                            dp_i = pre_i; id = abpoa_graph_index_to_node_id(graph, dp_i+beg_index); hit = 1;
                            dp_h = DP_H2E2F + dp_sn * (dp_i * 5); _dp_h = (int32_t*)dp_h;
                            break;
                        }
                    }
                }
                if (cur_op & ABPOA_E2_OP) {
                    _pre_dp_e2 = (int32_t*)(DP_H2E2F + dp_sn * ((pre_i * 5) + 2));
                    if (cur_op & ABPOA_M_OP) {
                        if (_dp_h[dp_j] == _pre_dp_e2[dp_j]) {
                            if (_pre_dp_h[dp_j] - gap_oe2 == _pre_dp_e2[dp_j]) cur_op = ABPOA_M_OP | ABPOA_F_OP;
                            else cur_op = ABPOA_E2_OP;
                            cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, dp_j-1);
                            dp_i = pre_i; id = abpoa_graph_index_to_node_id(graph, dp_i+beg_index); hit = 1;
                            dp_h = DP_H2E2F + dp_sn * (dp_i * 5); _dp_h = (int32_t*)dp_h;
                            break;
                        }
                    } else {
                        if (_dp_e2[dp_j] == _pre_dp_e2[dp_j] - gap_ext2) {
                            if (_pre_dp_h[dp_j] - gap_oe2 == _pre_dp_e2[dp_j]) cur_op = ABPOA_M_OP | ABPOA_F_OP;
                            else cur_op = ABPOA_E2_OP;
                            cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CDEL, 1, id, dp_j-1);
                            dp_i = pre_i; id = abpoa_graph_index_to_node_id(graph, dp_i+beg_index); hit = 1;
                            dp_h = DP_H2E2F + dp_sn * (dp_i * 5); _dp_h = (int32_t*)dp_h;
                            break;
                        }
                    }
                }
            }
        }
        if (hit == 0 && cur_op & ABPOA_F_OP) {
            if (cur_op & ABPOA_F1_OP) {
                _dp_f1 = (int32_t*)(dp_h + dp_sn * 3);
                if (cur_op & ABPOA_M_OP) {
                    if (_dp_h[dp_j] == _dp_f1[dp_j]) {
                        if (_dp_h[dp_j-1] - gap_oe1 == _dp_f1[dp_j]) cur_op = ABPOA_M_OP | ABPOA_E_OP, hit = 1;
                        else if (_dp_f1[dp_j-1] - gap_ext1 == _dp_f1[dp_j]) cur_op = ABPOA_F1_OP, hit = 1;
                    }
                } else {
                    if (_dp_h[dp_j-1] - gap_oe1 == _dp_f1[dp_j]) cur_op = ABPOA_M_OP | ABPOA_E_OP, hit = 1;
                    else if (_dp_f1[dp_j-1] - gap_ext1 == _dp_f1[dp_j]) cur_op = ABPOA_F1_OP, hit = 1;
                }
            }
            if (hit == 0 && cur_op & ABPOA_F2_OP) {
                _dp_f2 = (int32_t*)(dp_h + dp_sn * 4);
                if (cur_op & ABPOA_M_OP) {
                    if (_dp_h[dp_j] == _dp_f2[dp_j]) {
                        if (_dp_h[dp_j-1] - gap_oe2 == _dp_f2[dp_j]) cur_op = ABPOA_M_OP | ABPOA_E_OP, hit = 1;
                        else if (_dp_f2[dp_j-1] - gap_ext2 == _dp_f2[dp_j]) cur_op = ABPOA_F2_OP, hit = 1;
                    }
                } else {
                    if (_dp_h[dp_j-1] - gap_oe2 == _dp_f2[dp_j]) cur_op = ABPOA_M_OP | ABPOA_E_OP, hit = 1;
                    else if (_dp_f2[dp_j-1] - gap_ext2 == _dp_f2[dp_j]) cur_op = ABPOA_F2_OP, hit = 1;
                }
            }
            if (hit == 1) {
                cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, 1, id, dp_j-1); --dp_j;
                ++res->n_aln_bases;
            }
        }
        if (hit == 0 && (cur_op & ABPOA_M_OP) && (indel_first == 1)) {
            for (k = 0; k < pre_n[dp_i]; ++k) {
                pre_i = pre_index_i[k];
                if (dp_j-1 < dp_beg[pre_i] || dp_j-1 > dp_end[pre_i]) continue;
                _pre_dp_h = (int32_t*)(DP_H2E2F + dp_sn * (pre_i * 5));
                if (_pre_dp_h[dp_j-1] + s == _dp_h[dp_j]) { /* match/mismatch */
                    cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CMATCH, 1, id, dp_j-1);
                    dp_i = pre_i; --dp_j; id = abpoa_graph_index_to_node_id(graph, dp_i+beg_index); hit = 1;
                    dp_h = DP_H2E2F + dp_sn * (dp_i * 5); _dp_h = (int32_t*)dp_h;
                    cur_op = ABPOA_ALL_OP;
                    ++res->n_aln_bases; res->n_matched_bases += is_match ? 1 : 0;
                    indel_first = 0;
                    break;
                }
            }
        }
        if (hit == 0) exit(1);
#ifdef __DEBUG__
        fprintf(stderr, "%d, %d, %d\n", dp_i, dp_j, cur_op);
#endif
    }
    if (dp_j > 0) cigar = abpoa_push_cigar(&n_c, &m_c, cigar, ABPOA_CINS, dp_j, -1, dp_j-1);
    /* reverse cigar */
    res->graph_cigar = abpt->rev_cigar ? cigar : abpoa_reverse_cigar(n_c, cigar);
    res->n_cigar = n_c; res->m_cigar = m_c;
    res->node_e = abpoa_graph_index_to_node_id(graph, best_dp_i+beg_index), res->query_e = best_dp_j-1; /* 0-based */
    res->node_s = abpoa_graph_index_to_node_id(graph, _start_i+beg_index), res->query_s = _start_j-1;
    /*abpoa_print_cigar(n_c, *graph_cigar, graph);*/
}

int abpoa_cg_global_align_sequence_to_graph_core(abpoa_t *ab, int beg_node_id, int beg_index, int end_node_id, int end_index, uint8_t *index_map, int qlen, uint8_t *query, abpoa_para_t *abpt, SIMD_para_t sp, abpoa_res_t *res) {
    int tot_dp_sn = 0;
    abpoa_graph_t *graph = ab->abg; abpoa_simd_matrix_t *abm = ab->abm;
    int matrix_row_n = end_index-beg_index+1, matrix_col_n = qlen + 1;
    int **pre_index, *pre_n, _pre_index, _pre_n;
    int i, j, *dp_beg, *dp_beg_sn, *dp_end, *dp_end_sn, node_id, index_i, dp_i;
    int beg_sn, end_sn;
    int pn, log_n, size, qp_sn, dp_sn; /* pn: value per SIMDi, qp_sn/dp_sn/d_sn: segmented length*/
    SIMDi *dp_h, *qp, *qi;
    int32_t best_score = sp.inf_min, inf_min = sp.inf_min;
    int *mat = abpt->mat, best_i = 0, best_j = 0; int32_t gap_ext1 = abpt->gap_ext1;
    int w = abpt->wb < 0 ? qlen : abpt->wb + (int)(abpt->wf * qlen); /* when w < 0, do whole global */
    SIMDi zero = SIMDSetZeroi(), SIMD_INF_MIN = SIMDSetOnei32(inf_min);
    pn = sp.num_of_value; qp_sn = dp_sn = (matrix_col_n + pn - 1) / pn, log_n = sp.log_num, size = sp.size;
    qp = abm->s_mem;

    int32_t gap_open1 = abpt->gap_open1, gap_oe1 = gap_open1 + gap_ext1;
    int32_t gap_open2 = abpt->gap_open2, gap_ext2 = abpt->gap_ext2, gap_oe2 = gap_open2 + gap_ext2;
    SIMDi *DP_H2E2F, *dp_e1, *dp_e2, *dp_f2, *dp_f1;
    SIMDi GAP_O1 = SIMDSetOnei32(gap_open1), GAP_O2 = SIMDSetOnei32(gap_open2), GAP_E1 = SIMDSetOnei32(gap_ext1), GAP_E2 = SIMDSetOnei32(gap_ext2), GAP_OE1 = SIMDSetOnei32(gap_oe1), GAP_OE2 = SIMDSetOnei32(gap_oe2);
    DP_H2E2F = qp + qp_sn * abpt->m; qi = DP_H2E2F + dp_sn * matrix_row_n * 5;

    // for SET_F mask[pn], suf_min[pn], pre_min[logN]
    SIMDi *PRE_MASK, *SUF_MIN, *PRE_MIN, *GAP_E1S, *GAP_E2S;
    PRE_MASK = (SIMDi*)SIMDMalloc((pn+1) * size, size), SUF_MIN = (SIMDi*)SIMDMalloc((pn+1) * size, size), PRE_MIN = (SIMDi*)SIMDMalloc(pn * size, size), GAP_E1S =  (SIMDi*)SIMDMalloc(log_n * size, size), GAP_E2S =  (SIMDi*)SIMDMalloc(log_n * size, size);
    for (i = 0; i < pn; ++i) {
        int32_t *pre_mask = (int32_t*)(PRE_MASK + i);
        for (j = 0; j <= i; ++j) pre_mask[j] = -1;
        for (j = i+1; j < pn; ++j) pre_mask[j] = 0;
    } PRE_MASK[pn] = PRE_MASK[pn-1];
    SUF_MIN[0] = SIMDShiftLeft(SIMD_INF_MIN, SIMDShiftOneNi32);
    for (i = 1; i < pn; ++i) SUF_MIN[i] = SIMDShiftLeft(SUF_MIN[i-1], SIMDShiftOneNi32); SUF_MIN[pn] = SUF_MIN[pn-1];
    for (i = 1; i < pn; ++i) {
        int32_t *pre_min = (int32_t*)(PRE_MIN + i);
        for (j = 0; j < i; ++j) pre_min[j] = inf_min;
        for (j = i; j < pn; ++j) pre_min[j] = 0;
    }
    GAP_E1S[0] = GAP_E1; GAP_E2S[0] = GAP_E2;
    for (i = 1; i < log_n; ++i) {
        GAP_E1S[i] = SIMDAddi32(GAP_E1S[i-1], GAP_E1S[i-1]); GAP_E2S[i] = SIMDAddi32(GAP_E2S[i-1], GAP_E2S[i-1]);
    }
    abpoa_init_var(abpt, query, qlen, qp, qi, mat, qp_sn, pn, SIMD_INF_MIN);
    dp_beg = abm->dp_beg, dp_end = abm->dp_end, dp_beg_sn = abm->dp_beg_sn, dp_end_sn = abm->dp_end_sn;
    /* index of pre-node */
    pre_index = (int**)_err_calloc(matrix_row_n, sizeof(int*));
    pre_n = (int*)_err_calloc(matrix_row_n, sizeof(int));
    for (index_i = beg_index+1, dp_i = 1; index_i <= end_index; ++index_i, ++dp_i) {
        node_id = abpoa_graph_index_to_node_id(graph, index_i);
        pre_n[dp_i] = graph->node[node_id].in_edge_n;
        pre_index[dp_i] = (int*)_err_malloc(pre_n[dp_i] * sizeof(int));
        _pre_n = 0;
        for (j = 0; j < pre_n[dp_i]; ++j) {
            _pre_index = abpoa_graph_node_id_to_index(graph, graph->node[node_id].in_id[j]);
            if (index_map[_pre_index]) 
                pre_index[dp_i][_pre_n++] = _pre_index-beg_index;
        }
        pre_n[dp_i] = _pre_n;
    }
    abpoa_cg_first_dp(abpt, graph, index_map, beg_node_id, end_node_id, dp_beg, dp_end, dp_beg_sn, dp_end_sn, pn, qlen, w, dp_sn, DP_H2E2F, SIMD_INF_MIN, inf_min, gap_open1, gap_ext1, gap_open2, gap_ext2, gap_oe1, gap_oe2);

    for (index_i=beg_index+1, dp_i=1; index_i<end_index; ++index_i, ++dp_i) {
        if (index_map[index_i] == 0) continue;
        node_id = abpoa_graph_index_to_node_id(graph, index_i);
        SIMDi *q = qp + graph->node[node_id].base * qp_sn;
        dp_h = DP_H2E2F + (dp_i*5) * dp_sn; dp_e1 = dp_h + dp_sn; dp_e2 = dp_e1 + dp_sn; dp_f1 = dp_e2 + dp_sn; dp_f2 = dp_f1 + dp_sn;
        tot_dp_sn += abpoa_cg_dp(q, dp_h, dp_e1, dp_e2, dp_f1, dp_f2, pre_index, pre_n, index_i, dp_i, graph, abpt, dp_sn, pn, qlen, w, DP_H2E2F, SIMD_INF_MIN, GAP_O1, GAP_O2, GAP_E1, GAP_E2, GAP_OE1, GAP_OE2, GAP_E1S, GAP_E2S, PRE_MIN, PRE_MASK, SUF_MIN, log_n, dp_beg, dp_end, dp_beg_sn, dp_end_sn, end_node_id);
        if (abpt->wb >= 0) {
            beg_sn = dp_beg_sn[dp_i], end_sn = dp_end_sn[dp_i];
            int max_i = abpoa_max(SIMD_INF_MIN, zero, inf_min, dp_h, qi, qlen, pn, beg_sn, end_sn);
            abpoa_ada_max_i(max_i, graph, node_id);
        }
    }
    // printf("dp_sn: %d\n", tot_dp_sn);
    // printf("dp_sn: %d, node_n: %d, seq_n: %d\n", tot_dp_sn, graph->node_n, qlen);
    abpoa_global_get_max(graph, beg_index, end_node_id, index_map, DP_H2E2F, 5*dp_sn, qlen, dp_end, &best_score, &best_i, &best_j);
#ifdef __DEBUG__
    simd_abpoa_print_cg_matrix(int32_t, beg_index, end_index); fprintf(stderr, "best_score: (%d, %d) -> %d\n", best_i, best_j, best_score);
#endif
    res->best_score = best_score;
    abpoa_cg_backtrack(DP_H2E2F, pre_index, pre_n, dp_beg, dp_end, dp_sn, abpt->m, mat, gap_ext1, gap_ext2, gap_oe1, gap_oe2, beg_index, best_i, best_j, qlen, graph, abpt, query, res);
    for (i = 0; i < matrix_row_n; ++i) free(pre_index[i]);
    free(pre_index); free(pre_n);
    SIMDFree(PRE_MASK); SIMDFree(SUF_MIN); SIMDFree(PRE_MIN);
    SIMDFree(GAP_E1S); SIMDFree(GAP_E2S);
    return best_score;
}

// align query to subgraph between beg_node_id and end_node_id (both are excluded)
// generally: beg/end are the SRC/SINK_node
int simd_abpoa_align_sequence_to_subgraph(abpoa_t *ab, abpoa_para_t *abpt, int beg_node_id, int end_node_id, uint8_t *query, int qlen, abpoa_res_t *res) {
    // if (abpt->simd_flag == 0) err_fatal_simple("No SIMD instruction available.");

    int i, j, beg_index = ab->abg->node_id_to_index[beg_node_id], end_index = ab->abg->node_id_to_index[end_node_id];
    int gn = end_index - beg_index + 1;
    uint8_t *index_map = (uint8_t*)_err_calloc(ab->abg->node_n, sizeof(uint8_t));
    index_map[beg_index] = index_map[end_index] = 1;
    for (i = beg_index; i < end_index-1; ++i) {
        if (index_map[i] == 0) continue;
        int node_id = abpoa_graph_index_to_node_id(ab->abg, i);
        int out_n = ab->abg->node[node_id].out_edge_n;
        for (j = 0; j < out_n; ++j) {
            int out_id = ab->abg->node[node_id].out_id[j];
            index_map[abpoa_graph_node_id_to_index(ab->abg, out_id)] = 1;
        }
    }

#ifdef __DEBUG__
    _simd_p32.inf_min = MAX_OF_TWO(abpt->gap_ext1, abpt->gap_ext2) * 31 +MAX_OF_THREE(INT32_MIN + abpt->min_mis, INT32_MIN + abpt->gap_open1 + abpt->gap_ext1, INT32_MIN + abpt->gap_open2 + abpt->gap_ext2);
    if (simd_abpoa_realloc(ab, gn, qlen, abpt, _simd_p32)) return 0;
    if (abpt->gap_mode == ABPOA_CONVEX_GAP) abpoa_cg_global_align_sequence_to_graph_core(ab, beg_node_id, beg_index, end_node_id, end_index, index_map, qlen, query, abpt, _simd_p32, res);
#else
    int32_t max_score, bits, mem_ret=0, gap_ext1 = abpt->gap_ext1, gap_ext2 = abpt->gap_ext2;
    int32_t gap_oe1 = abpt->gap_open1+gap_ext1, gap_oe2 = abpt->gap_open2+gap_ext2;
    // if (abpt->simd_flag & SIMD_AVX512F && !(abpt->simd_flag & SIMD_AVX512BW))
    //    max_score = INT16_MAX + 1; // AVX512F has no 8/16 bits operations
    // else {
    int len = qlen > gn ? qlen : gn;
    max_score = MAX_OF_TWO(qlen * abpt->max_mat, len * abpt->gap_ext1 + abpt->gap_open1);
    // }
    if (max_score <= INT16_MAX - abpt->min_mis - gap_oe1 - gap_oe2) {
        _simd_p16.inf_min = MAX_OF_THREE(INT16_MIN + abpt->min_mis, INT16_MIN + gap_oe1, INT16_MIN + gap_oe2) + 31 * MAX_OF_TWO(gap_ext1, gap_ext2);
        mem_ret = simd_abpoa_realloc(ab, gn, qlen, abpt, _simd_p16);
        bits = 16;
    } else {
        _simd_p32.inf_min = MAX_OF_THREE(INT32_MIN + abpt->min_mis, INT32_MIN + gap_oe1, INT32_MIN + gap_oe2) + 31 * MAX_OF_TWO(gap_ext1, gap_ext2);
        mem_ret = simd_abpoa_realloc(ab, gn, qlen, abpt, _simd_p32);
        bits = 32;
    }
    if (mem_ret) return 0;

    if (bits == 16) {
        if (abpt->gap_mode == ABPOA_LINEAR_GAP) {
            simd_abpoa_lg_align_sequence_to_graph_core(int16_t, _simd_p16, SIMDSetOnei16, SIMDMaxi16, \
                    SIMDAddi16, SIMDSubi16, SIMDShiftOneNi16, SIMDSetIfGreateri16, SIMDGetIfGreateri16);
        } else if (abpt->gap_mode == ABPOA_AFFINE_GAP) {
            simd_abpoa_ag_align_sequence_to_graph_core(int16_t, _simd_p16, SIMDSetOnei16, SIMDMaxi16, \
                    SIMDAddi16, SIMDSubi16, SIMDShiftOneNi16, SIMDSetIfGreateri16, SIMDGetIfGreateri16, SIMDSetIfEquali16);
        } else if (abpt->gap_mode == ABPOA_CONVEX_GAP) {
            simd_abpoa_cg_align_sequence_to_graph_core(int16_t, _simd_p16, SIMDSetOnei16, SIMDMaxi16, \
                    SIMDAddi16, SIMDSubi16, SIMDShiftOneNi16, SIMDSetIfGreateri16, SIMDGetIfGreateri16, SIMDSetIfEquali16);
        }
    } else { // 2147483647, DP_H/E/F: 32 bits
        if (abpt->gap_mode == ABPOA_LINEAR_GAP) {
            simd_abpoa_lg_align_sequence_to_graph_core(int32_t, _simd_p32, SIMDSetOnei32, SIMDMaxi32, \
                    SIMDAddi32, SIMDSubi32, SIMDShiftOneNi32, SIMDSetIfGreateri32, SIMDGetIfGreateri32);
        } else if (abpt->gap_mode == ABPOA_AFFINE_GAP) {
            simd_abpoa_ag_align_sequence_to_graph_core(int32_t, _simd_p32, SIMDSetOnei32, SIMDMaxi32, \
                    SIMDAddi32, SIMDSubi32, SIMDShiftOneNi32, SIMDSetIfGreateri32, SIMDGetIfGreateri32, SIMDSetIfEquali32);
        } else if (abpt->gap_mode == ABPOA_CONVEX_GAP) {
            simd_abpoa_cg_align_sequence_to_graph_core(int32_t, _simd_p32, SIMDSetOnei32, SIMDMaxi32, \
                    SIMDAddi32, SIMDSubi32, SIMDShiftOneNi32, SIMDSetIfGreateri32, SIMDGetIfGreateri32, SIMDSetIfEquali32);
        }
    }
#endif
    free(index_map);
    return 0;
}

int simd_abpoa_align_sequence_to_graph(abpoa_t *ab, abpoa_para_t *abpt, uint8_t *query, int qlen, abpoa_res_t *res) {
    return simd_abpoa_align_sequence_to_subgraph(ab, abpt, ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID, query, qlen, res);
}
