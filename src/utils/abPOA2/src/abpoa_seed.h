#ifndef _ABPOA_SEED_H
#define _ABPOA_SEED_H

#include <stdint.h>
#include <stdlib.h>
#include "abpoa.h"

// emulate 128-bit integers and arrays
typedef struct { uint64_t x, y; } ab_u128_t;
typedef struct { size_t n, m; ab_u128_t *a; } ab_u128_v;

typedef struct { size_t n, m; uint64_t *a; } ab_u64_v;

#ifdef __cplusplus
extern "C" {
#endif

int abpoa_build_guide_tree_partition(uint8_t **seqs, int *seq_lens, int n_seq, abpoa_para_t *abpt, int *read_id_map, ab_u64_v *par_anchors, int *par_c);

#ifdef __cplusplus
}
#endif


#endif
