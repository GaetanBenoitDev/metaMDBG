#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "abpoa_seq.h"
#include "abpoa_align.h"
#include "abpoa_graph.h"
#include "utils.h"
#include "kstring.h"
#include "khash.h"

KHASH_MAP_INIT_STR(abstr, uint32_t)

// for nt
// AaCcGgTtNn ==> 0,1,2,3,4
unsigned char ab_nt4_table[256] = {
       0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// 65,97=>A, 67,99=>C, 71,103=>G, 84,85,116,117=>T, else=>N
const char ab_nt256_table[256] = {
       'A', 'C', 'G', 'T',  'N', '-', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', '-',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
};

// for aa
// AaCcGgTtNn ... ==> 0,1,2,3,4 ...
// BbDdEeFf   ... ==> 5,6,7,8 ...
unsigned char ab_aa26_table[256] = {
	 0,  1,  2,  3,   4,  5,  6,  7,   8,  9, 10, 11,  12, 13, 14, 15, 
	16, 17, 18, 19,  20, 21, 22, 23,  24, 25, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26,  0,  5,  1,   6,  7,  8,  2,   9, 10, 11, 12,  13, 14,  4, 15, 
	16, 17, 18, 19,   3, 20, 21, 22,  23, 24, 25, 26,  26, 26, 26, 26, 
	26,  0,  5,  1,   6,  7,  8,  2,   9, 10, 11, 12,  13, 14,  4, 15, 
	16, 17, 18, 19,   3, 20, 21, 22,  23, 24, 25, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26
};

// 0/1/2/3/4=>ACGTN
// 5/6/7/8=>BDEF ...
const char ab_aa256_table[256] = {
	'A', 'C', 'G', 'T',  'N', 'B', 'D', 'E',  'F', 'H', 'I',  'J', 'K', 'L', 'M', 'O',
	'P', 'Q', 'R', 'S',  'U', 'V', 'W', 'X',  'Y', 'Z', '*', '-',  '*', '*', '*', '*',
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', 'A', 'B', 'C',  'D', 'E', 'F', 'G',  'H', 'I', 'J', 'K',  'L', 'M', 'N', 'O', 
	'P', 'Q', 'R', 'S',  'T', 'U', 'V', 'W',  'X', 'Y', 'Z', '*',  '*', '*', '*', '*', 
	'*', 'A', 'B', 'C',  'D', 'E', 'F', 'G',  'H', 'I', 'J', 'K',  'L', 'M', 'N', 'O', 
	'P', 'Q', 'R', 'S',  'T', 'U', 'V', 'W',  'X', 'Y', 'Z', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*'
};

char ab_char26_table[256];
char ab_char256_table[256];

abpoa_seq_t *abpoa_init_seq(void) {
    abpoa_seq_t *abs = (abpoa_seq_t*)_err_malloc(sizeof(abpoa_seq_t));
    abs->n_seq = 0; abs->m_seq = CHUNK_READ_N;
    abs->seq = (abpoa_str_t*)_err_calloc(abs->m_seq, sizeof(abpoa_str_t));
    abs->name = (abpoa_str_t*)_err_calloc(abs->m_seq, sizeof(abpoa_str_t));
    abs->comment = (abpoa_str_t*)_err_calloc(abs->m_seq, sizeof(abpoa_str_t));
    abs->qual = (abpoa_str_t*)_err_calloc(abs->m_seq, sizeof(abpoa_str_t));
    abs->is_rc = (uint8_t*)_err_calloc(abs->m_seq, sizeof(uint8_t));
    return abs;
}

void abpoa_free_seq(abpoa_seq_t *abs) {
    int i;
    for (i = 0; i < abs->m_seq; ++i) {
        if (abs->seq[i].m > 0) free(abs->seq[i].s);
        if (abs->name[i].m > 0) free(abs->name[i].s);
        if (abs->comment[i].m > 0) free(abs->comment[i].s);
        if (abs->qual[i].m > 0) free(abs->qual[i].s);
    }
    free(abs->seq); free(abs->name); free(abs->comment); free(abs->qual); 
    free(abs->is_rc); free(abs);
}

void abpoa_cpy_str(abpoa_str_t *str, char *s, int l) {
    if (l > 0) {
        str->l = l; str->m = l + 1;
        str->s = (char*)_err_malloc(str->m * sizeof(char));
        memcpy(str->s, s, l);
        str->s[str->l] = 0;
    }
}

void abpoa_cpy_seq(abpoa_seq_t *abs, int seq_i, kseq_t *kseq) {
    abpoa_cpy_str(abs->seq+seq_i, kseq->seq.s, kseq->seq.l);
    abpoa_cpy_str(abs->name+seq_i, kseq->name.s, kseq->name.l);
    abpoa_cpy_str(abs->comment+seq_i, kseq->comment.s, kseq->comment.l);
    abpoa_cpy_str(abs->qual+seq_i, kseq->qual.s, kseq->qual.l);
}

abpoa_seq_t *abpoa_realloc_seq(abpoa_seq_t *abs) {
    if (abs->n_seq >= abs->m_seq) {
        int m_seq = MAX_OF_TWO(abs->n_seq, abs->m_seq << 1);
        abs->seq = (abpoa_str_t*)_err_realloc(abs->seq, m_seq * sizeof(abpoa_str_t));
        abs->name = (abpoa_str_t*)_err_realloc(abs->name, m_seq * sizeof(abpoa_str_t));
        abs->comment = (abpoa_str_t*)_err_realloc(abs->comment, m_seq * sizeof(abpoa_str_t));
        abs->qual = (abpoa_str_t*)_err_realloc(abs->qual, m_seq * sizeof(abpoa_str_t));
        abs->is_rc = (uint8_t*)_err_realloc(abs->is_rc, m_seq * sizeof(uint8_t));
        int i;
        for (i = abs->m_seq; i < m_seq; ++i) {
            abs->seq[i].l = abs->seq[i].m = 0;
            abs->name[i].l = abs->name[i].m = 0;
            abs->comment[i].l = abs->comment[i].m = 0;
            abs->qual[i].l = abs->qual[i].m = 0;
            abs->is_rc[i] = 0;
        }
        abs->m_seq = m_seq;
    }
    return abs;
}

int abpoa_read_nseq(abpoa_seq_t *abs, kseq_t *kseq, int chunk_read_n) {
    int n = 0;
    while (n < chunk_read_n && kseq_read(kseq) >= 0) {
        abpoa_realloc_seq(abs);
        // copy kseq to abs->seq
        abpoa_cpy_seq(abs, abs->n_seq, kseq);
        abs->n_seq++; n++;
    }
    return n;
}

int abpoa_read_seq(abpoa_seq_t *abs, kseq_t *kseq) {
    int n = 0;
    while (kseq_read(kseq) >= 0) {
        abpoa_realloc_seq(abs);
        // copy kseq to abs->seq
        abpoa_cpy_seq(abs, abs->n_seq, kseq);
        abs->n_seq++; n++;
    }
    return n;
}

static long int _strtol10(const char *str, char **endptr) {
    long int res = 0; unsigned d;
    char *s;
    for (s = (char*)str, d = s[0]-'0'; d < 10; ++s, d = s[0]-'0')
        res = res * 10 + d;
    if (endptr) *endptr = s;
    return res;
}

static unsigned long int _strtoul10(const char *str, char **endptr) {
    unsigned long int res = 0; unsigned d;
    char *s;
    for (s = (char*)str, d = s[0]-'0'; d < 10; ++s, d = s[0]-'0')
        res = res * 10 + d;
    if (endptr) *endptr = s;
    return res;
}

int gfa_aux_parse(char *s, uint8_t **data, int *max)
{
	char *q, *p;
	kstring_t str;
	if (s == 0) return 0;
	str.l = 0, str.m = *max, str.s = (char*)*data;
	if (*s == '\t') ++s;
	for (p = q = s;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (p - q >= 5 && q[2] == ':' && q[4] == ':' && (q[3] == 'I' || q[3] == 'A' || q[3] == 'i' || q[3] == 'f' || q[3] == 'Z' || q[3] == 'B')) {
				int type = q[3];
				kputsn_(q, 2, &str);
				q += 5;
				if (type == 'A') {
					kputc_('A', &str);
					kputc_(*q, &str);
                } else if (type == 'I') {
					uint32_t x;
					// x = strtol(q, &q, 10);
                    x = _strtoul10(q, &q);
					kputc_(type, &str); kputsn_((char*)&x, 4, &str);
				} else if (type == 'i') {
					int32_t x;
					// x = strtol(q, &q, 10);
                    x = _strtol10(q, &q);
					kputc_(type, &str); kputsn_((char*)&x, 4, &str);
				} else if (type == 'f') {
					float x;
					x = strtod(q, &q);
					kputc_('f', &str); kputsn_(&x, 4, &str);
				} else if (type == 'Z') {
					kputc_('Z', &str); kputsn_(q, p - q + 1, &str); // note that this include the trailing NULL
				} else if (type == 'B') {
					type = *q++; // q points to the first ',' following the typing byte
					if (p - q >= 2 && (type == 'c' || type == 'C' || type == 's' || type == 'S' || type == 'i' || type == 'I' || type != 'f')) {
						int32_t n;
						char *r;
						for (r = q, n = 0; *r; ++r)
							if (*r == ',') ++n;
						kputc_('B', &str); kputc_(type, &str); kputsn_(&n, 4, &str);
						// TODO: to evaluate which is faster: a) aligned array and then memmove(); b) unaligned array; c) kputsn_()
						if (type == 'c')      while (q + 1 < p) { int8_t   x = strtol(q + 1, &q, 0); kputc_(x, &str); }
						else if (type == 'C') while (q + 1 < p) { uint8_t  x = strtol(q + 1, &q, 0); kputc_(x, &str); }
						else if (type == 's') while (q + 1 < p) { int16_t  x = strtol(q + 1, &q, 0); kputsn_(&x, 2, &str); }
						else if (type == 'S') while (q + 1 < p) { uint16_t x = strtol(q + 1, &q, 0); kputsn_(&x, 2, &str); }
						// else if (type == 'i') while (q + 1 < p) { int32_t  x = strtol(q + 1, &q, 0); kputsn_(&x, 4, &str); }
						else if (type == 'i') while (q + 1 < p) { int32_t  x = _strtol10(q + 1, &q); kputsn_(&x, 4, &str); }
						// else if (type == 'I') while (q + 1 < p) { uint32_t x = strtol(q + 1, &q, 0); kputsn_(&x, 4, &str); }
						else if (type == 'I') while (q + 1 < p) { uint32_t x = _strtoul10(q + 1, &q); kputsn_(&x, 4, &str); }
						else if (type == 'f') while (q + 1 < p) { float    x = strtod(q + 1, &q);    kputsn_(&x, 4, &str); }
					}
				} // should not be here, as we have tested all types
			}
			q = p + 1;
			if (c == 0) break;
		}
	}
	if (str.l > 0 && str.l == str.m) ks_resize(&str, str.l + 1);
	if (str.s) str.s[str.l] = 0;
	*max = str.m, *data = (uint8_t*)str.s;
	return str.l;
}

static inline int gfa_aux_type2size(int x)
{
	if (x == 'C' || x == 'c' || x == 'A') return 1;
	else if (x == 'S' || x == 's') return 2;
	else if (x == 'I' || x == 'i' || x == 'f') return 4;
	else return 0;
}

#define __skip_tag(s) do { \
		int type = toupper(*(s)); \
		++(s); \
		if (type == 'Z') { while (*(s)) ++(s); ++(s); } \
		else if (type == 'B') (s) += 5 + gfa_aux_type2size(*(s)) * (*(int32_t*)((s)+1)); \
		else (s) += gfa_aux_type2size(type); \
	} while(0)

uint8_t *gfa_aux_get(int l_data, const uint8_t *data, const char tag[2])
{
	const uint8_t *s = data;
	int y = tag[0]<<8 | tag[1];
	while (s < data + l_data) {
		int x = (int)s[0]<<8 | s[1];
		s += 2;
		if (x == y) return (uint8_t*)s;
		__skip_tag(s);
	}
	return 0;
}

int gfa_aux_del(int l_data, uint8_t *data, uint8_t *s)
{
	uint8_t *p;
	p = s - 2;
	__skip_tag(s);
	memmove(p, s, l_data - (s - data));
	return l_data - (s - p);
}

int abpoa_gfa_parse_H(abpoa_graph_t *abg, int *n_s, int *n_l, int *n_p, char *s) {
    if (s[1] != '\t' || s[2] == '0') return -1;
    int l_aux, m_aux = 0; uint8_t *aux = 0, *info;
    l_aux = gfa_aux_parse(s + 2, &aux, &m_aux);

    info = gfa_aux_get(l_aux, aux, "NS");
    if (info == 0 || info[0] != 'i') err_fatal_simple("Error: no \"NS\" tag in GFA header.");
    *n_s = *(int32_t*)(info+1);
    abg->node_m = *n_s + 2;
    abg->node = (abpoa_node_t*)_err_realloc(abg->node, abg->node_m * sizeof(abpoa_node_t*));
    l_aux = gfa_aux_del(l_aux, aux, info);

    info = gfa_aux_get(l_aux, aux, "NL");
    if (info == 0 || info[0] != 'i') err_fatal_simple("Error: no \"NL\" tag in GFA header.");
    *n_l = *(int32_t*)(info+1);
    l_aux = gfa_aux_del(l_aux, aux, info);

    info = gfa_aux_get(l_aux, aux, "NP");
    if (info == 0 || info[0] != 'i') err_fatal_simple("Error: no \"NP\" tag in GFA header.");
    *n_p = *(int32_t*)(info+1);
    l_aux = gfa_aux_del(l_aux, aux, info);

    if (aux) free(aux);
    return 0;
}

typedef struct {
    int n, m;
    kstring_t *seq, *name;
    khash_t(abstr) *h;
} seg_seq_t;

seg_seq_t *seg_seq_init(void) {
    seg_seq_t *s = (seg_seq_t*)_err_malloc(sizeof(seg_seq_t));
    s->n = s->m = 0; s->seq = 0, s->name = 0;
    s->h = kh_init(abstr);
    return s;
}

seg_seq_t *seg_seq_realloc(seg_seq_t *r) {
    if (r->n >= r->m) {
        int m;
        if (r->m == 0) m = 1;
        else m = MAX_OF_TWO(r->n, (r->m) << 1);
        r->seq = (kstring_t*)_err_realloc(r->seq, m * sizeof(kstring_t));
        r->name = (kstring_t*)_err_realloc(r->name, m * sizeof(kstring_t));
        int i;
        for (i = r->m; i < m; ++i) {
            r->seq[i] = (kstring_t){0,0,0};
            r->name[i] = (kstring_t){0,0,0};
        }
        r->m = m;
    }
    return r;
}

void seg_seq_free(seg_seq_t *s) {
    if (s->m > 0) {
        int i;
        for (i = 0; i < s->m; ++i) {
            if (s->seq[i].m) free(s->seq[i].s);
            if (s->name[i].m) free(s->name[i].s);
        }
        free(s->seq); free(s->name);
    }
    kh_destroy(abstr, s->h);
    free(s);
}

int abpoa_gfa_parse_S(seg_seq_t *segs,  char *s) {
    if (s[1] != '\t' || s[2] == '\0') return -1;
    char *deli_s, *info_s, *seq = 0;
    int i, seq_len, seg_name_len, is_ok = 0;
    char *seg_name=0;

    for (i = 0, deli_s = info_s = s + 2;; ++deli_s) {
        if (*deli_s == 0 || *deli_s == '\t') {
            int c = *deli_s;
            *deli_s = 0;
            if (i == 0) {
                seg_name = info_s;
                seg_name_len = deli_s - info_s;
            } else if (i == 1) {
                seq = info_s;
                seq_len = deli_s - info_s;
                is_ok = 1;
                break;
            }
            if (c == 0) break;
            ++i, info_s = deli_s + 1;
        }
    }

    if (is_ok) {
        seg_seq_realloc(segs);
        kputsn(seg_name, seg_name_len, segs->name+segs->n);
        kputsn(seq, seq_len, segs->seq+segs->n);
        int absent;
        khint_t pos = kh_put(abstr, segs->h, segs->name[segs->n].s, &absent);
        if (absent) kh_val(segs->h, pos) = segs->n;
        else err_fatal(__func__, "Duplicated chromosome: \"%s\".", seg_name);
        ++segs->n;
    } else err_fatal(__func__, "Error: no seq in GFA segment line (%s).", seg_name);
    return 0;
}

/*int abpoa_gfa_parse_S(abpoa_graph_t *abg, char *s) {
    if (s[1] != '\t' || s[2] == '\0') return -1;
    char *deli_s, *info_s, *seq = 0;
    int i, seq_len, is_ok = 0;
    char *seg_name=0;

    for (i = 0, deli_s = info_s = s + 2;; ++deli_s) {
        if (*deli_s == 0 || *deli_s == '\t') {
            int c = *deli_s;
            *deli_s = 0;
            if (i == 0) {
                seg_name = info_s;
                abpoa_realloc_seq(seg_names);
                abpoa_cpy_str(seg_names->name+seg_names->n_seq, seg_name, deli_s - info_s);
                seg_names->n_seq++;
            } else if (i == 1) {
                seq = info_s;
                seq_len = deli_s - info_s;
                is_ok = 1;
                break;
            }
            if (c == 0) break;
            ++i, info_s = deli_s + 1;
        }
    }

    if (is_ok) {
        int seg_id, absent;
        for (i = 0; i < seq_len; ++i) {
            seg_id = abpoa_add_graph_node(abg, ab_char26_table[(int)(seq[i])]);
            if (i == 0) {
                khint_t pos = kh_put(str, seg_name2in_id, seg_names->name[seg_names->n_seq-1].s, &absent);
                if (absent) kh_val(seg_name2in_id, pos) = seg_id;
                else err_fatal(__func__, "Error: duplicated seg name (%s).", seg_name);
            }
            if (i == seq_len-1) {
                khint_t pos = kh_put(str, seg_name2out_id,  seg_names->name[seg_names->n_seq-1].s, &absent);
                if (absent) kh_val(seg_name2out_id, pos) = seg_id;
                else err_fatal(__func__, "Error: duplicated seg name (%s).", seg_name);
            }
        }
    } else err_fatal(__func__, "Error: no seq in GFA segment line (%s).", seg_name);
    return 0;
}*/

int abpoa_gfa_parse_P(abpoa_graph_t *abg, abpoa_seq_t *abs, seg_seq_t *segs, int add_read_id, int p_i, int p_n, khash_t(abstr) *seg_name2in_id, khash_t(abstr) *seg_name2out_id, char *s) {
    if (s[1] != '\t' || s[2] == '\0') return -1;
    char *deli_s, *info_s, *path = 0;
    int i, is_ok = 0, is_rc = -1;
    char *path_name=0; int path_name_len=0;
    kstring_t *seg_seq, *seg_name; int read_ids_n = 1 + ((p_n-1) >> 6);

    for (i = 0, deli_s = info_s = s + 2;; ++deli_s) {
        if (*deli_s == 0 || *deli_s == '\t') {
            int c = *deli_s;
            *deli_s = 0;
            if (i == 0) {
                path_name = info_s;
                path_name_len = deli_s - info_s;
            } else if (i == 1) {
                path = info_s;
                is_ok = 1;
                break;
            }
            if (c == 0) break;
            ++i, info_s = deli_s + 1;
        }
    }

    if (is_ok) {
        char *deli_s, *info_s, *_seg_name; khint_t pos, seg_pos; int absent;
        int id, in_id=-1, out_id=-1, last_id = ABPOA_SRC_NODE_ID, next_id = ABPOA_SINK_NODE_ID;
        for (deli_s = info_s = path; ; ++deli_s) {
            if (*deli_s == '+') {
                if (is_rc == 1) err_fatal(__func__, "Error: path has both \'+\' and \'-\' seg. (%s)", path_name);
                is_rc = 0; *deli_s = 0; _seg_name = info_s;
                seg_pos = kh_get(abstr, segs->h, _seg_name);
                if (seg_pos == kh_end(segs->h)) err_fatal(__func__, "Error: seg (%s) not exist.", info_s);
                seg_name = segs->name + kh_val(segs->h, seg_pos);
                seg_seq = segs->seq + kh_val(segs->h, seg_pos);

                // check if seg already exist
                pos = kh_put(abstr, seg_name2in_id, seg_name->s, &absent);
                if (absent) { // add node for seg_seq
                    for (i = 0; i < (int)seg_seq->l; ++i) {
                        id = abpoa_add_graph_node(abg, ab_char26_table[(int)(seg_seq->s[i])]);
                        if (i == 0) in_id = id;
                        if (i == (int)seg_seq->l-1) out_id = id;
                    }
                    kh_val(seg_name2in_id, pos) = in_id;
                    pos = kh_put(abstr, seg_name2out_id, seg_name->s, &absent);
                    kh_val(seg_name2out_id, pos) = out_id;
                } else {
                    in_id = kh_val(seg_name2in_id, pos);
                    pos = kh_put(abstr, seg_name2out_id, seg_name->s, &absent);
                    out_id = kh_val(seg_name2out_id, pos);
                }
                // add edge
                abpoa_add_graph_edge(abg, last_id, in_id, 1, 1, add_read_id, p_i, read_ids_n);
                if (in_id < out_id) {
                    for (i = 0; i < out_id - in_id; ++i)
                        abpoa_add_graph_edge(abg, in_id+i, in_id+i+1, 1, 1, add_read_id, p_i, read_ids_n);
                } else if (in_id > out_id) err_fatal(__func__, "Error: in_id (%d) > out_id (%d).", in_id, out_id);

                last_id = out_id;
                info_s = deli_s + 2;
            } else if (*deli_s == '-') {
                if (is_rc == 0) err_fatal(__func__, "Error: path has both \'+\' and \'-\' seg. (%s)", path_name);
                is_rc = 1; *deli_s = 0; _seg_name = info_s;
                seg_pos = kh_get(abstr, segs->h, _seg_name);
                if (seg_pos == kh_end(segs->h)) err_fatal(__func__, "Error: seg (%s) not exist.", info_s);
                seg_name = segs->name + kh_val(segs->h, seg_pos);
                seg_seq = segs->seq + kh_val(segs->h, seg_pos);

                // check if seg exist
                pos = kh_put(abstr, seg_name2in_id, seg_name->s, &absent);
                if (absent) { // add node for seg_seq
                    for (i = 0; i < (int)seg_seq->l; ++i) {
                        id = abpoa_add_graph_node(abg, ab_char26_table[(int)(seg_seq->s[i])]);
                        if (i == 0) in_id = id;
                        if (i == (int)seg_seq->l-1) out_id = id;
                    }
                    kh_val(seg_name2in_id, pos) = in_id;
                    pos = kh_put(abstr, seg_name2out_id, seg_name->s, &absent);
                    kh_val(seg_name2out_id, pos) = out_id;
                } else {
                    in_id = kh_val(seg_name2in_id, pos); out_id = kh_val(seg_name2out_id, pos);
                }

                // add edge
                abpoa_add_graph_edge(abg, out_id, next_id, 1, 1, add_read_id, p_i, read_ids_n);
                if (in_id < out_id) {
                    for (i = 0; i < out_id - in_id; ++i)
                        abpoa_add_graph_edge(abg, in_id+i, in_id+i+1, 1, 1, add_read_id, p_i, read_ids_n);
                } else if (in_id > out_id) err_fatal(__func__, "Error: in_id (%d) > out_id (%d).", in_id, out_id);

                next_id = in_id;
                info_s = deli_s + 2;
            } else if (*deli_s == 0 || *deli_s == '\t') break;
        }
        if (is_rc) abpoa_add_graph_edge(abg, ABPOA_SRC_NODE_ID, next_id, 1, 1, add_read_id, p_i, read_ids_n);
        else abpoa_add_graph_edge(abg, last_id, ABPOA_SINK_NODE_ID, 1, 1, add_read_id, p_i, read_ids_n);
        // set abs
        abpoa_realloc_seq(abs);
        abpoa_cpy_str(abs->name+abs->n_seq, path_name, path_name_len); 
        abs->is_rc[abs->n_seq] = is_rc; abs->n_seq++;
    } else err_fatal(__func__, "Error: no path in GFA path line (%s).", path_name);
    return 0;
}

int abpoa_fa_parse_seq(abpoa_graph_t *abg, abpoa_seq_t *abs, kstring_t *seq, kstring_t *name, int add_read_id, int p_i, int p_n, int **rank2node_id) {
    if (*rank2node_id == 0) {
        *rank2node_id = (int*)_err_calloc(seq->l, sizeof(int));
    } 
    char *s = seq->s;
    int32_t read_ids_n = 1 + ((p_n-1) >> 6);
    int32_t i, rank, last_id = ABPOA_SRC_NODE_ID, cur_id, aln_id; uint8_t base;
    for (i = 0; s[i]; ++i) {
        if (s[i] == '-') continue; // gap
        else {
            base = ab_char26_table[(int)(s[i])];
            rank = i;
            cur_id = (*rank2node_id)[rank];
            if (cur_id == 0) {
                cur_id = abpoa_add_graph_node(abg, base);
                (*rank2node_id)[rank] = cur_id;
            } else {
                if (abg->node[cur_id].base != base) {
                    aln_id = abpoa_get_aligned_id(abg, cur_id, base);
                    if (aln_id == -1) {
                        aln_id = abpoa_add_graph_node(abg, base);
                        abpoa_add_graph_aligned_node(abg, cur_id, aln_id);
                    }
                    cur_id = aln_id;
                }
            }
            abpoa_add_graph_edge(abg, last_id, cur_id, 1, 1, add_read_id, p_i, read_ids_n);
            last_id = cur_id;
        }
    }
    abpoa_add_graph_edge(abg, last_id, ABPOA_SINK_NODE_ID, 1, 1, add_read_id, p_i, read_ids_n);
    abpoa_realloc_seq(abs);
    abpoa_cpy_str(abs->name + abs->n_seq, name->s, name->l); abs->n_seq++;
    return 0;
}

abpoa_t *abpoa_restore_graph(abpoa_t *ab, abpoa_para_t *abpt) {
    char *fn = abpt->incr_fn;
    if (fn == NULL) return ab;
    gzFile fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r"); if (fp == 0) return NULL;
    kstream_t *ks = ks_init(fp); kstring_t s={0,0,0}; int dret, line_n=0, read_name_e; 
    seg_seq_t *seqs = seg_seq_init(); 
    khash_t(abstr) *seg_name2in_id = kh_init(abstr), *seg_name2out_id = kh_init(abstr);
    int add_read_id = abpt->use_read_ids;
    int p_i = -1, is_fa = 0, *rank2node_id=0;

    abpoa_graph_t *abg = ab->abg; abpoa_seq_t *abs = ab->abs;
    while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
        line_n++;
        int ret = 0; 
        if (is_fa) {
            if (s.l > 0 && s.s[0] == '>') { // name
                // parse_seq
                if (seqs->seq[seqs->n].l > 0) {
                    ret = abpoa_fa_parse_seq(abg, abs, seqs->seq+seqs->n, seqs->name+seqs->n, add_read_id, p_i, p_i+1, &rank2node_id);
                    seqs->n++;
                }
                // kputsn seqs->name

                read_name_e = 1;
                while (read_name_e < (int)s.l && !isspace(s.s[read_name_e])) read_name_e++;
                seg_seq_realloc(seqs);
                kputsn(s.s+1, read_name_e-1, seqs->name + seqs->n);
                p_i++;
            } else { // seq
                // kputsn seqs->seq
                kputsn(s.s, s.l, seqs->seq + seqs->n);
            }
        } else {
            if (s.l > 0 && s.s[0] == '>') {
                read_name_e = 1;
                while (read_name_e < (int)s.l && !isspace(s.s[read_name_e])) read_name_e++;
                seg_seq_realloc(seqs);
                kputsn(s.s+1, read_name_e-1, seqs->name + seqs->n);
                is_fa = 1; p_i++;
            } 
            // else if (s.l < 2 || s.s[0] == '#') continue; // comment
            // else if (s.s[0] == 'H') ret = abpoa_gfa_parse_H(abg, &s_n, &l_n, &p_n, s.s);
            else if (s.s[0] == 'S') ret = abpoa_gfa_parse_S(seqs, s.s); // include Link information
            // else if (s.s[0] == 'L') ret = abpoa_gfa_parse_L(abg, seg_name, s.s);
            else if (s.s[0] == 'P') {
                p_i++;
                ret = abpoa_gfa_parse_P(abg, abs, seqs, add_read_id, p_i, p_i+1, seg_name2in_id, seg_name2out_id, s.s);
            }
        }
        if (ret < 0) err_fatal(__func__, "Error in %c-line at line %ld (error code %d)", s.s[0], (long)line_n, ret);
    }
    if (is_fa) { // last seq
        abpoa_fa_parse_seq(abg, abs, seqs->seq+seqs->n, seqs->name+seqs->n, add_read_id, p_i, p_i+1, &rank2node_id);
        seqs->n++;
    }
    if (s.m) free(s.s);
    ks_destroy(ks); gzclose(fp);
    seg_seq_free(seqs); kh_destroy(abstr, seg_name2in_id); kh_destroy(abstr, seg_name2out_id);
    if (rank2node_id) free(rank2node_id);
    if (abs->n_seq == 0) {
        err_func_printf(__func__, "Warning: no graph/sequence restored from file \'%s\'.\n", fn);
        abg->node_n = 2;
    }
    abg->is_called_cons = abg->is_set_msa_rank = abg->is_topological_sorted = 0;
    return ab;
}
