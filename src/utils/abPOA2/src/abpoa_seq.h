#ifndef _ABPOA_SEQ_H
#define _ABPOA_SEQ_H
#include <zlib.h>
#include "abpoa.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

#ifdef __cplusplus
extern "C" {
#endif

abpoa_seq_t *abpoa_realloc_seq(abpoa_seq_t *abs);
void abpoa_cpy_str(abpoa_str_t *str, char *s, int l);
abpoa_seq_t *abpoa_init_seq(void);
void abpoa_free_seq(abpoa_seq_t *abs);
int abpoa_read_seq(abpoa_seq_t *abs, kseq_t *kseq);
#ifdef __cplusplus
}
#endif


#endif
