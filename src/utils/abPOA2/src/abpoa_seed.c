#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "abpoa.h"
#include "abpoa_seed.h"
#include "utils.h"
#include "kvec.h"
#include "ksort.h"

const char LogTable256[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
    -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
    LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
    LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2_32(uint32_t v) {
    uint32_t t, tt;
    if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
    return (t = v>>8) ? 8 + LogTable256[t] : LogTable256[v];
}

#define ab_sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(ab_128x, ab_u128_t, ab_sort_key_128x, 8)

#define ab_sort_key_128y(a) ((a).y)
KRADIX_SORT_INIT(ab_128y, ab_u128_t, ab_sort_key_128y, 8)

#define ab_sort_key_64(a) (a)
KRADIX_SORT_INIT(64, uint64_t, ab_sort_key_64, 8)

/* from lh3/minimap2/sketch.c */
/********** start *************/
static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

typedef struct { // a simplified version of kdq
    int front, count;
    int a[32];
} tiny_queue_t;

static inline void tq_push(tiny_queue_t *q, int x)
{
    q->a[((q->count++) + q->front) & 0x1f] = x;
}

static inline int tq_shift(tiny_queue_t *q)
{
    int x;
    if (q->count == 0) return -1;
    x = q->a[q->front++];
    q->front &= 0x1f;
    --q->count;
    return x;
}

/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      minimizers
 *               p->a[i].x = kMer<<8 | kmerSpan
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
void mm_sketch(void *km, const uint8_t *str, int len, int w, int k, uint32_t rid, int is_hpc, int both_strand, ab_u128_v *p)
{
    uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
    int i, j, l, buf_pos, min_pos, kmer_span = 0;
    ab_u128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
    tiny_queue_t tq;

    assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
    memset(buf, 0xff, w * 16);
    memset(&tq, 0, sizeof(tiny_queue_t));
    kv_resize(ab_u128_t, km, *p, p->n + len/w);

    for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
        int c = str[i];
        ab_u128_t info = { UINT64_MAX, UINT64_MAX };
        if (c < 4) { // not an ambiguous base
            uint32_t z;
            if (is_hpc) {
                int skip_len = 1;
                if (i + 1 < len && str[i + 1] == c) {
                    for (skip_len = 2; i + skip_len < len; ++skip_len)
                        if (str[i + skip_len] != c)
                            break;
                    i += skip_len - 1; // put $i at the end of the current homopolymer run
                }
                tq_push(&tq, skip_len);
                kmer_span += skip_len;
                if (tq.count > k) kmer_span -= tq_shift(&tq);
            } else kmer_span = l + 1 < k? l + 1 : k;
            if (both_strand) {
                kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
                kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
                if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
                z = kmer[0] < kmer[1]? 0 : 1; // strand
            } else {
                kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
                z = 0;
            }
            ++l;
            if (l >= k && kmer_span < 256) {
                info.x = hash64(kmer[z], mask) << 8 | kmer_span;
                info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
            }
        } else l = 0, tq.count = tq.front = 0, kmer_span = 0;
        buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
        if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
            for (j = buf_pos + 1; j < w; ++j)
                if (min.x == buf[j].x && buf[j].y != min.y) kv_push(ab_u128_t, km, *p, buf[j]);
            for (j = 0; j < buf_pos; ++j)
                if (min.x == buf[j].x && buf[j].y != min.y) kv_push(ab_u128_t, km, *p, buf[j]);
        }
        if (info.x <= min.x) { // a new minimum; then write the old min
            if (l >= w + k && min.x != UINT64_MAX) kv_push(ab_u128_t, km, *p, min);
            min = info, min_pos = buf_pos;
        } else if (buf_pos == min_pos) { // old min has moved outside the window
            if (l >= w + k - 1 && min.x != UINT64_MAX) kv_push(ab_u128_t, km, *p, min);
            for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
                if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
            for (j = 0; j <= buf_pos; ++j)
                if (min.x >= buf[j].x) min = buf[j], min_pos = j;
            if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
                for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
                    if (min.x == buf[j].x && min.y != buf[j].y) kv_push(ab_u128_t, km, *p, buf[j]);
                for (j = 0; j <= buf_pos; ++j)
                    if (min.x == buf[j].x && min.y != buf[j].y) kv_push(ab_u128_t, km, *p, buf[j]);
            }
        }
        if (++buf_pos == w) buf_pos = 0;
    }
    if (min.x != UINT64_MAX)
        kv_push(ab_u128_t, km, *p, min);
}

// For amino acid sequence
void mm_aa_sketch(void *km, const uint8_t *str, int len, int w, int k, uint32_t rid, int is_hpc, ab_u128_v *p)
{
    uint64_t mask = (1ULL<<5*k) - 1, kmer[2] = {0,0};
    int i, j, l, buf_pos, min_pos, kmer_span = 0;
    ab_u128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
    tiny_queue_t tq;

    assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 11)); // 56 / 5 == 11
    memset(buf, 0xff, w * 16);
    memset(&tq, 0, sizeof(tiny_queue_t));
    kv_resize(ab_u128_t, km, *p, p->n + len/w);

    for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
        int c = str[i];
        ab_u128_t info = { UINT64_MAX, UINT64_MAX };
        if (c < 26) { // not an ambiguous base
            uint32_t z;
            if (is_hpc) {
                int skip_len = 1;
                if (i + 1 < len && str[i + 1] == c) {
                    for (skip_len = 2; i + skip_len < len; ++skip_len)
                        if (str[i + skip_len] != c)
                            break;
                    i += skip_len - 1; // put $i at the end of the current homopolymer run
                }
                tq_push(&tq, skip_len);
                kmer_span += skip_len;
                if (tq.count > k) kmer_span -= tq_shift(&tq);
            } else kmer_span = l + 1 < k? l + 1 : k;
            
            kmer[0] = (kmer[0] << 5 | c) & mask; // only forward k-mer for aa seq
            z = 0;
            ++l;
            if (l >= k && kmer_span < 256) {
                info.x = hash64(kmer[z], mask) << 8 | kmer_span;
                info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
            }
        } else l = 0, tq.count = tq.front = 0, kmer_span = 0;
        buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
        if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
            for (j = buf_pos + 1; j < w; ++j)
                if (min.x == buf[j].x && buf[j].y != min.y) kv_push(ab_u128_t, km, *p, buf[j]);
            for (j = 0; j < buf_pos; ++j)
                if (min.x == buf[j].x && buf[j].y != min.y) kv_push(ab_u128_t, km, *p, buf[j]);
        }
        if (info.x <= min.x) { // a new minimum; then write the old min
            if (l >= w + k && min.x != UINT64_MAX) kv_push(ab_u128_t, km, *p, min);
            min = info, min_pos = buf_pos;
        } else if (buf_pos == min_pos) { // old min has moved outside the window
            if (l >= w + k - 1 && min.x != UINT64_MAX) kv_push(ab_u128_t, km, *p, min);
            for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
                if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
            for (j = 0; j <= buf_pos; ++j)
                if (min.x >= buf[j].x) min = buf[j], min_pos = j;
            if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
                for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
                    if (min.x == buf[j].x && min.y != buf[j].y) kv_push(ab_u128_t, km, *p, buf[j]);
                for (j = 0; j <= buf_pos; ++j)
                    if (min.x == buf[j].x && min.y != buf[j].y) kv_push(ab_u128_t, km, *p, buf[j]);
            }
        }
        if (++buf_pos == w) buf_pos = 0;
    }
    if (min.x != UINT64_MAX)
        kv_push(ab_u128_t, km, *p, min);
}
/************ end *************/

// tree_id_map: guide tree node id -> original input order id
/* mm: is unsorted
 * a.x = kMer<<8 | kmerSpan
 * a.y = rid<<32 | strand<<31 | lastPos
 */
int abpoa_build_guide_tree(abpoa_para_t *abpt, int n_seq, ab_u128_v *mm, int *tree_id_map) {
    if (mm->n == 0) return 0;

    if (abpt->verbose > 0) fprintf(stderr, "[%s] Building progressive guide tree ... ", __func__);
    size_t i, _i, j; int rid1, rid2;                                          // mm_hit_n: mimizer hits between each two sequences
                                                                              // 0: 0
                                                                              // 1: 0 1
                                                                              // 2: 0 1 2
    int *mm_hit_n = (int*)_err_calloc((n_seq * (n_seq+1)) >> 1, sizeof(int)); //  ...
                                                                              // n: 0 1 ... n-1 n
                                                                              // 
                                                                              // # total mimizers of i: mm_hit_n[(i*(i+1))/2+i]
                                                                              // # total hits for i and j (i>j): mm_hit_n[(i*(i+1)/2)+j]
    radix_sort_ab_128x(mm->a, mm->a + mm->n); // sort mm by k-mer hash values
    uint64_t last_x = mm->a[0].x;
    int *mm_cnt = (int*)_err_malloc(n_seq * sizeof(int));
    for (_i = 0, i = 1; i < mm->n; ++i) { // collect mm hits
        if (mm->a[i].x != last_x) {
            // now [_i, i-1] have the same minimizer k-mer
            memset(mm_cnt, 0, n_seq * sizeof(int));
            for (j = _i; j < i; ++j) {
                // count mm->a[j]
                rid1 = mm->a[j].y >> 32;
                ++mm_cnt[rid1];
                ++mm_hit_n[((rid1 * (rid1+1)) >> 1) + rid1];
            }
            for (rid1 = 0; rid1 < n_seq-1; ++rid1) {
                for (rid2 = rid1+1; rid2 < n_seq; ++rid2) {
                    mm_hit_n[((rid2 * (rid2 + 1)) >> 1) + rid1] += MIN_OF_TWO(mm_cnt[rid1], mm_cnt[rid2]);
                }
            }
            // next minimizer
            last_x = mm->a[i].x, _i = i;
        }
    }
    // now [_i, i-1] have the same minimizer k-mer
    memset(mm_cnt, 0, n_seq * sizeof(int));
    for (j = _i; j < i; ++j) {
        // count mm->a[j]
        rid1 = mm->a[j].y >> 32;
        ++mm_cnt[rid1];
        ++mm_hit_n[((rid1 * (rid1+1)) >> 1) + rid1];
    }
    for (rid1 = 0; rid1 < n_seq-1; ++rid1) {
        for (rid2 = rid1+1; rid2 < n_seq; ++rid2) {
            mm_hit_n[((rid2 * (rid2 + 1)) >> 1) + rid1] += MIN_OF_TWO(mm_cnt[rid1], mm_cnt[rid2]);
        }
    }
    free(mm_cnt);

    // calculate jaccard similarity between each two sequences
    double *jac_sim = (double*)_err_calloc((n_seq * (n_seq-1)) >> 1, sizeof(double));
    double max_jac = -1.0, jac; int max_i=-1, max_j=-1;
    for (i = 1; i < (size_t)n_seq; ++i) {
        for (j = 0; j < i; ++j) {
            int tot_n = mm_hit_n[((i*(i+1))>>1)+i] + mm_hit_n[((j*(j+1))>>1)+j] - mm_hit_n[((i*(i+1))>>1)+j];
            if (tot_n == 0) jac = 0;
            else if (tot_n < 0) err_fatal(__func__, "Bug in progressive tree building. (1)");
            else jac = (0.0+mm_hit_n[((i*(i+1))>>1)+j]) / tot_n;
            jac_sim[((i * (i-1)) >> 1) + j] = jac; // jac_sim[rid2*(rid2-1)/2 + rid1]
            if (jac > max_jac) {
                max_jac = jac; max_i = i, max_j = j;
            }
        }
    }

    // build guide tree
    // first pick two with the biggest jac (max_i, max_j)
    int n_in_map = 2; tree_id_map[0] = max_j, tree_id_map[1] = max_i;

    // then, pick one with biggest jac sum with existing sequence in tree_id_map
    while (n_in_map < n_seq) {
        max_jac = -1.0, max_i = n_seq;
        for (rid1 = 0; rid1 < n_seq; ++rid1) {
            jac = 0.0;
            for (i = 0; i < (size_t)n_in_map; ++i) {
                rid2 = tree_id_map[i];
                if (rid1 == rid2) { jac = -1.0; break; }
                else if (rid1 > rid2) jac += jac_sim[((rid1 * (rid1-1)) >> 1) + rid2];
                else jac += jac_sim[((rid2 * (rid2-1)) >> 1) + rid1];
            }
            if (jac > max_jac) {
                max_jac = jac;
                max_i = rid1;
            }
        }
        if (max_i == n_seq) err_fatal(__func__, "Bug in progressive tree building. (2)");
        tree_id_map[n_in_map++] = max_i;
    }

    free(mm_hit_n); free(jac_sim);
    if (abpt->verbose > 0) fprintf(stderr, "done!\n");
    return 0;
}

// t's mm: is sorted, q's mm is unsorted
//       | r1's minimizers | r2's minimizers | ... | rn's minimizers |
// mm_c: | 0 | n_r1_mm | n_r1..2_mm | ... | n_r1..n-1_mm | n_r1..n_mm |
// t is already in the graph, q is query sequence
// merge sort for t and q's minimizer buckets
int collect_anchors1(void *km, ab_u64_v *anchors, ab_u128_v mm, int *mm_c, int tid, int qid, int qlen, int k) {
    int i, j, _i, _j; uint64_t xi, xj, _xi, _xj, _yi, _yj, a;
    i = mm_c[tid], j = mm_c[qid];
    // t's mm is already sorted XXX
    radix_sort_ab_128x(mm.a + j, mm.a + mm_c[qid+1]);

    while (i < mm_c[tid+1] && j < mm_c[qid+1]) {
        xi = mm.a[i].x, xj = mm.a[j].x;
        if (xi == xj) {
            for (_i = i; _i < mm_c[tid+1]; ++_i) {
                _xi = mm.a[_i].x;
                if (_xi != xi) break;
                _yi = mm.a[_i].y;
                for (_j = j; _j < mm_c[qid+1]; ++_j) {
                    _xj = mm.a[_j].x;
                    if (_xj != xj) break;
                    _yj = mm.a[_j].y;
                    // t_strand<<63 | t_lastPos<<32 | q_lastPos
                    if ((_yi & 1) == (_yj & 1)) { // same strand
                        a = (uint64_t)((uint32_t)_yi>>1)<<32 | ((uint32_t)_yj>>1);
                    } else { // different strand
                        a = 1ULL<<63 | (uint64_t)((uint32_t)_yi>>1)<<32 | (qlen - (((uint32_t)_yj>>1)+1-k) - 1); // XXX qlen < pow(2,28)
                    }
                    kv_push(uint64_t, km, *anchors, a);
                }
            }
            i = _i, j = _j;
        } else if (xi < xj) ++i;
        else if (xi > xj) ++j;
    }
    // sort by tpos
    radix_sort_64(anchors->a, anchors->a + anchors->n);
    return anchors->n;
}

int get_local_chain_score(int j_end_tpos, int j_end_qpos, int i_end_anchor_i, ab_u64_v *anchors, int *pre_id, int *score) {
    int i = i_end_anchor_i, chain_score = 0;
    int i_tpos, i_qpos;
    do {
        i_tpos = (anchors->a[i] >> 32) & 0x7fffffff, i_qpos = (int32_t)anchors->a[i];

        if (i_tpos <= j_end_tpos && i_qpos <= j_end_qpos) break;
        i = pre_id[i];
    } while (i != -1);

    if (i == -1) chain_score = score[i_end_anchor_i];
    else chain_score = score[i_end_anchor_i] - score[i];
    return chain_score;
}


// local chains:
//   x: strand | end_tpos | end_qpos
//   y: end_anchor_i | start_anchor_i
int abpoa_dp_chaining_of_local_chains(void *km, ab_u128_t *local_chains, int n_local_chains, ab_u64_v *anchors, int *score, int *pre_id, ab_u64_v *par_anchors, int min_w, int tlen, int qlen) {
    int i, j, st, score1, global_max_score=INT32_MIN, global_max_i=-1;
    int *chain_score = (int*)kmalloc(km, n_local_chains * 4), *pre_chain_id = (int*)kmalloc(km, n_local_chains * 4);
    size_t _n = par_anchors->n;

    for (i = st = 0; i < n_local_chains; ++i) {
        uint64_t ix = local_chains[i].x, iy = local_chains[i].y;
        int istrand = ix >> 63, i_end_qpos = (int32_t)ix, i_end_anchor_i = iy >> 32, i_start_anchor_i = (int32_t)iy;
        int i_start_tpos = (anchors->a[i_start_anchor_i] >> 32) & 0x7fffffff, i_start_qpos = (int32_t)anchors->a[i_start_anchor_i];
        int max_j = -1, max_score = score[i_end_anchor_i];
        while (st < i) {
            if ((local_chains[st].x) >> 63 != istrand) ++st;
            else break;
        }
        for (j = i-1; j >= st; --j) {
            uint64_t jx = local_chains[j].x;
            int j_end_tpos = (jx >> 32) & 0x7fffffff, j_end_qpos = (int32_t)jx; //, j_end_anchor_i = iy >> 32;
            if (j_end_qpos >= i_end_qpos) continue;

            if (i_start_tpos > j_end_tpos && i_start_qpos > j_end_qpos) score1 = chain_score[j] + score[i_end_anchor_i];
            else score1 = chain_score[j] + get_local_chain_score(j_end_tpos, j_end_qpos, i_end_anchor_i, anchors, pre_id, score);

            if (score1 > max_score) {
                max_score = score1; max_j = j;
            }
        }
        chain_score[i] = max_score; pre_chain_id[i] = max_j;
        if (max_score > global_max_score) {
            global_max_score = max_score;
            global_max_i = i;
        }
    }
    if (global_max_i < 0) return 0;
    // collect anchors based on global_max_i
    int cur_i = global_max_i, pre_i = pre_chain_id[global_max_i];
    uint64_t cur_y = local_chains[cur_i].y, pre_x, pre_y;
    int last_tpos=tlen, last_qpos=qlen;
    while (pre_i != -1) { // collect valid anchors in local_chains[cur_i], constrained by local_chains[pre_i]
        pre_x = local_chains[pre_i].x, pre_y = local_chains[pre_i].y;
        int pre_end_tpos = (pre_x >> 32) & 0x7fffffff, pre_end_qpos = (int32_t)pre_x;
        i = cur_y >> 32;
        while (i != -1) {
            int cur_tpos = (anchors->a[i] >> 32) & 0x7fffffff, cur_qpos = (int32_t)anchors->a[i];
            if (cur_tpos > pre_end_tpos && cur_qpos > pre_end_qpos) {
                if (last_tpos - cur_tpos >= min_w && last_qpos - cur_qpos >= min_w) {
                    kv_push(uint64_t, 0, *par_anchors, anchors->a[i]);
                    last_tpos = cur_tpos, last_qpos = cur_qpos;
                }
            } else break;
            i = pre_id[i];
        }
        cur_i = pre_i, pre_i = pre_chain_id[pre_i], cur_y = pre_y;
    }
    // collect anchors of last chain: local_chains[cur_i]
    i = cur_y >> 32;
    while (i != -1) {
        int cur_tpos = (anchors->a[i] >> 32) & 0x7fffffff, cur_qpos = (int32_t)anchors->a[i];
        if (last_tpos - cur_tpos >= min_w && last_qpos - cur_qpos >= min_w) {
            kv_push(uint64_t, 0, *par_anchors, anchors->a[i]);
            last_tpos = cur_tpos, last_qpos = cur_qpos;
        }
        i = pre_id[i];
    }
    // reverse order of par_anchors
    for (i = 0; i < (int)(par_anchors->n-_n) >> 1; ++i) {
        uint64_t tmp = par_anchors->a[_n+i];
        par_anchors->a[_n+i] = par_anchors->a[par_anchors->n-i-1];
        par_anchors->a[par_anchors->n-i-1] = tmp;
    }

#ifdef __DEBUG__
    for (i = _n; i < par_anchors->n; ++i) {
        uint64_t ia = par_anchors->a[i];
        // strand, rpos, qpos
        fprintf(stderr, "%c\t%ld\t%d\n", "+-"[ia >> 63], (ia>>32) & 0x7fffffff, ((uint32_t)ia));
    }
#endif
    kfree(km, chain_score), kfree(km, pre_chain_id);
    return 0;
}

// for DP chaining
static int get_chain_score(int max_bw, int *score, int i_qpos, int i_tpos, int j_qpos, int j_tpos, int k) {
    int delta_q, delta_t, delta_tq, min_d;

    delta_q = i_qpos - j_qpos; delta_t = i_tpos - j_tpos;
    min_d = MIN_OF_THREE(delta_q, delta_t, k);
    *score = min_d;
    if (delta_q >= delta_t) {
        if ((delta_tq = delta_q - delta_t) > max_bw) return 0;
    } else {
        if ((delta_tq = delta_t - delta_q) > max_bw) return 0;
    }
    *score -= ((ilog2_32(delta_tq) >> 1) + delta_tq * 0.01 * k);
    return 1;
}

// Dynamic Programming-based Chaining for global alignment mode
// anchors:
//          strand<<63 | tpos<<32 | qpos
int abpoa_dp_chaining(void *km, ab_u64_v *anchors, ab_u64_v *par_anchors, abpoa_para_t *abpt, int tlen, int qlen) {
    int i, j, st, n_a = anchors->n;
    int *score = (int32_t*)kmalloc(km, n_a * 4), *pre_id = (int32_t*)kmalloc(km, n_a * 4), *end_pos = (int32_t*)kmalloc(km, n_a * 4);
    memset(end_pos, 0, n_a * 4);

    int max_bw = 100, max_dis = 100, max_skip_anchors = 25, max_non_best_anchors = 50, min_local_chain_score = 100;
    int min_w = abpt->min_w+abpt->k;
    int i_qpos, i_tpos, i_tstrand, j_qpos, j_tpos;
    for (i = st = 0; i < (int)anchors->n; ++i) {
        uint64_t ia = anchors->a[i];
        i_qpos = (int32_t)ia, i_tpos = (ia >> 32) & 0x7fffffff, i_tstrand = ia >> 63;
        int max_j = -1, n_skip=0, non_best_iter_n = 0, max_score=abpt->k, _score;
        while (st < i) {
            uint64_t st_a = anchors->a[st];
            if ((st_a >> 63) != i_tstrand || (int)((st_a >> 32) & 0x7fffffff) + max_dis < i_tpos) ++st;
            else break;
        }

        for (j = i-1; j >= st; --j) { // check if j is i's optimal pre anchor
            uint64_t ja = anchors->a[j];
            j_qpos = (uint32_t)ja; j_tpos = (ja >> 32) & 0x7fffffff;
            if (j_qpos >= i_qpos || j_qpos + max_dis < i_qpos) continue;
            if (!get_chain_score(max_bw, &_score, i_qpos, i_tpos, j_qpos, j_tpos, abpt->k)) continue;
            _score += score[j];
            if (_score > max_score) {
                max_score = _score; max_j = j;
                non_best_iter_n = 0;
                if (n_skip > 0) --n_skip;
            } else if (end_pos[j] == i) {
                if (++n_skip > max_skip_anchors) break;
            } else if (++non_best_iter_n > max_non_best_anchors) break;

            if (pre_id[j] >= 0) end_pos[pre_id[j]] = i;
        }
#ifdef __DEBUG__
        fprintf(stderr, "%d pre_id: %d, score: %d, tpos: %d, qpos: %d\n", i, max_j, max_score, i_tpos, i_qpos);
#endif
        score[i] = max_score, pre_id[i] = max_j;
    }

    memset(end_pos, 0, n_a * 4);
    int n_local_chains = 0;
    for (i = n_a-1; i >= 0; --i) {
        if (pre_id[i] >= 0) end_pos[pre_id[i]] = 1;
        if (end_pos[i] == 0 && score[i] >= min_local_chain_score) {
            end_pos[i] = 2;
            ++n_local_chains;
        }
    }
    // collect local chains
    // x: score
    // y: e_a_i
    ab_u128_t *local_chains = (ab_u128_t*)kmalloc(km, n_local_chains * sizeof(ab_u128_t));
    for (i = n_local_chains = 0; i < n_a; ++i) {
        if (end_pos[i] == 2) {
            local_chains[n_local_chains].x = score[i];
            local_chains[n_local_chains++].y = i;
        }
    }
    radix_sort_ab_128x(local_chains, local_chains + n_local_chains);

    // collect local chains
    // x: strand | endpos | score
    // y: s_a_i | e_a_i
    int32_t *anchor_map = end_pos; memset(anchor_map, 0, n_a * 4);
    int start_id, end_id, tot_chain_i; uint64_t strand, tpos, qpos;

    for (i = tot_chain_i = n_local_chains-1; i >=0; --i) {
        j = local_chains[i].y; end_id = j; strand = anchors->a[i] >> 63; tpos = (anchors->a[j] >> 32) & 0x7fffffff, qpos = (int32_t)anchors->a[j];
        do {
            start_id = j;
            anchor_map[j] = 1;
            j = pre_id[j];
        } while (j >= 0 && anchor_map[j] == 0);

        if (j < 0) { // reach the start of the chain
            local_chains[tot_chain_i].x = strand << 63 | tpos << 32 | qpos;
            local_chains[tot_chain_i--].y = (uint64_t)end_id << 32 | start_id;
        } 
        // not keep branched chains
        /*else if ((int32_t)local_chains[i].x - score[j] >= min_local_chain_score) { // anchor_map == 1, anchor was already used in other chain
            local_chains[tot_chain_i].x = strand << 63 | tpos << 32 | qpos;
            local_chains[tot_chain_i--].y = (uint64_t)end_id << 32 | ((int32_t)local_chains[i].x - score[j]);
            pre_id[start_id] = -1;
        }*/
    }

    radix_sort_ab_128x(local_chains+tot_chain_i+1, local_chains + n_local_chains);
    abpoa_dp_chaining_of_local_chains(km, local_chains+tot_chain_i+1, n_local_chains-1-tot_chain_i, anchors, score, pre_id, par_anchors, min_w, tlen, qlen);

    kfree(km, score); kfree(km, pre_id); kfree(km, end_pos); kfree(km, local_chains);
    return 0;
}

int bin_search_min_larger(int *lis, int left, int right, int key) {
    int mid;
    while (right - left > 1) {
        mid = ((right - left) >> 1) + left;

        if (lis[mid] >= key) right = mid;
        else left = mid;
    }
    return right;
}

// rank: qpos<<32 | tpos_rank
int LIS(void *km, int tot_n, uint64_t *rank, int n) {
    int *pre_rank = (int*)kcalloc(km, tot_n+1, sizeof(int));
    int *lis = (int*)kmalloc(km, n * sizeof(int));
    int i, n_lis, irank, idx;
    lis[0] = (uint32_t)rank[0]; n_lis = 1;

    // calculate LIS length
    for (i = 1; i < n; ++i) {
        irank = (uint32_t)rank[i];
        if (irank < lis[0]) {
            lis[0] = irank;
        } else if (irank > lis[n_lis-1]) {
            lis[n_lis] = irank;
            pre_rank[irank] = lis[n_lis-1];
            ++n_lis;
        } else {
            idx = bin_search_min_larger(lis, -1, n_lis-1, irank);
            lis[idx] = irank;
            if (idx > 0) pre_rank[irank] = lis[idx-1];
        }
    }
    // collect LIS, store ids in rank
    int r = lis[n_lis-1]; i = n_lis-1;
    while (r != 0) {
        if (i < 0) err_fatal_simple("Error in LIS.");
        rank[i--] = r;
        r = pre_rank[r];
    }
    kfree(km, pre_rank); kfree(km, lis);
    return n_lis;
}

// XXX TODO use dp-based chaining
// XXX TODO remove q_span

// Longest Increasing Subsequence-based Chaining (only works for global alignment mode)
// input:
//   anchors: (sorted by tpos)
//          strand<<63 | tpos<<32 | qpos
// output:
//   anchor list size: n
//   list of anchors: anchors
int LIS_chaining(void *km, ab_u64_v *anchors, ab_u64_v *par_anchors, int min_w) {
    size_t i, j, n_a = anchors->n, n_for, n_rev;
    uint64_t *for_rank = (uint64_t*)kmalloc(km, sizeof(uint64_t) * n_a);
    uint64_t *rev_rank = (uint64_t*)kmalloc(km, sizeof(uint64_t) * n_a);
    uint64_t qpos;

    n_for = 0, n_rev = 0;
    for (i = 0; i < n_a; ++i) {
        uint64_t ia = anchors->a[i];
        qpos = (uint32_t)ia;
        if (ia >> 63) { // reverse
            rev_rank[n_rev++] = qpos << 32 | (i+1);
        } else { // forward
            for_rank[n_for++] = qpos << 32 | (i+1);
        }
    }

    if (n_for > 0) {
        radix_sort_64(for_rank, for_rank + n_for);
        n_for = LIS(km, n_a, for_rank, n_for);
    }
    if (n_rev > 0) {
        radix_sort_64(rev_rank, rev_rank + n_rev);
        n_rev = LIS(km, n_a, rev_rank, n_rev);
    }

    size_t n; uint64_t *rank;
    if (n_for > n_rev) {
        n = n_for; rank = for_rank; kfree(km, rev_rank);
    } else {
        n = n_rev; rank = rev_rank; kfree(km, for_rank);
    }
    // filter anchors
    int last_tpos = -1, last_qpos = -1, cur_tpos, cur_qpos;
#ifdef __DEBUG__
    size_t _n = par_anchors->n;
#endif
    for (i = 0; i < n; ++i) {
        j = (int)rank[i]-1;
        cur_tpos = (anchors->a[j] >> 32) & 0x7fffffff;
        if (cur_tpos - last_tpos < min_w) continue;
        cur_qpos = (uint32_t)anchors->a[j];
        if (cur_qpos - last_qpos < min_w) continue;

        kv_push(uint64_t, 0, *par_anchors, anchors->a[j]); // store LIS anchors into par_anchors
        last_tpos = cur_tpos; last_qpos = cur_qpos;
    }
#ifdef __DEBUG__
    for (i = _n; i < par_anchors->n; ++i) {
        uint64_t ia = par_anchors->a[i];
        // strand, rpos, qpos
        fprintf(stderr, "%c\t%ld\t%d\n", "+-"[ia >> 63], (ia>>32) & 0x7fffffff, ((uint32_t)ia));
    }
#endif
    return 0;
}

int abpoa_collect_mm(void *km, uint8_t **seqs, int *seq_lens, int n_seq, abpoa_para_t *abpt, ab_u128_v *mm, int *mm_c) {
    if (abpt->verbose > 0) fprintf(stderr, "[%s] Collecting minimizers ... ", __func__);
    int i;
    mm_c[0] = 0;
    for (i = 0; i < n_seq; ++i) { // collect minimizers
        if (abpt->m > 5) mm_aa_sketch(km, seqs[i], seq_lens[i], abpt->w, abpt->k, i, 0, mm);
        else mm_sketch(km, seqs[i], seq_lens[i], abpt->w, abpt->k, i, 0, abpt->amb_strand, mm);
        mm_c[i+1] = mm->n;
    }
    if (abpt->verbose > 0) fprintf(stderr, "done!\n");
    return mm->n;
}

// split guide tree and seeding and partition
int abpoa_build_guide_tree_partition(uint8_t **seqs, int *seq_lens, int n_seq, abpoa_para_t *abpt, int *read_id_map, ab_u64_v *par_anchors, int *par_c) {
    int i; void *km = km_init();
    for (i = 0; i < n_seq; ++i) read_id_map[i] = i;
    ab_u128_v mm1 = {0, 0, 0}; int *mm_c = (int*)_err_malloc((n_seq+1) * sizeof(int));
    abpoa_collect_mm(km, seqs, seq_lens, n_seq, abpt, &mm1, mm_c);

    if (abpt->progressive_poa && n_seq > 2) {
        // copy mm1 to mm2
        ab_u128_v mm2 = {0, 0, 0};
        for (i = 0; i < (int)mm1.n; ++i) kv_push(ab_u128_t, km, mm2, mm1.a[i]);
        // use mm2 to build guide tree
        abpoa_build_guide_tree(abpt, n_seq, &mm2, read_id_map);
        kfree(km, mm2.a);
    }
    if (abpt->disable_seeding || n_seq < 2) {
        kfree(km, mm1.a); free(mm_c); km_destroy(km);
        return 0; // no anchor
    }
    // partition into small windows
    int qid, tid;
    tid = read_id_map[0];
    radix_sort_ab_128x(mm1.a + mm_c[tid], mm1.a + mm_c[tid+1]);

    par_c[0] = 0;
    for (i = 1; i < n_seq; ++i) {
        tid = read_id_map[i-1]; qid = read_id_map[i];
        ab_u64_v anchors = {0, 0, 0};
        // collect minimizer hit anchors between t and q
        collect_anchors1(km, &anchors, mm1, mm_c, tid, qid, seq_lens[qid], abpt->k);
        // filtering and only keep LIS anchors
#ifdef __DEBUG__
        fprintf(stderr, "%d vs %d (tot_n: %ld)\n", tid, qid, anchors.n);
#endif
        // alignment mode: different chaining result for global/local/extend
        abpoa_dp_chaining(km, &anchors, par_anchors, abpt, seq_lens[tid], seq_lens[qid]);
        par_c[i] = par_anchors->n;
        kfree(km, anchors.a);
    }

    kfree(km, mm1.a); free(mm_c); km_destroy(km);
    return 0; // par_anchors->n;
}
