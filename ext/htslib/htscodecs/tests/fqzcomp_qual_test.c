/* Tests for fqz codec */
/*
 * Copyright (c) 2019-2021 Genome Research Ltd.
 * Author(s): James Bonfield
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *    1. Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
 *       Institute nor the names of its contributors may be used to endorse
 *       or promote products derived from this software without specific
 *       prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
 * LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include "config.h"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>
#include <ctype.h>
#include <limits.h>

#include "htscodecs/fqzcomp_qual.h"
#include "htscodecs/varint.h"

#ifndef MAX_REC
#define MAX_REC 1000000
#endif

#ifndef MAX_SEQ
#  define MAX_SEQ 100000
#endif

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#  define MAX(a,b) ((a)>(b)?(a):(b))
#endif

static fqz_slice fixed_slice = {0};

fqz_slice *fake_slice(size_t buf_len, int *len, int *r2, int *sel, int nlen) {
    fixed_slice.num_records = (nlen == 1) ? (buf_len+len[0]-1) / len[0] : nlen;
    assert(fixed_slice.num_records <= MAX_REC);
    int i;
    if (!fixed_slice.len)
        fixed_slice.len = malloc(MAX_REC * sizeof(*fixed_slice.len));
    if (!fixed_slice.flags)
        fixed_slice.flags = malloc(MAX_REC * sizeof(*fixed_slice.flags));
    for (i = 0; i < fixed_slice.num_records; i++) {
        int idx = i < nlen ? i : nlen-1;
        fixed_slice.len[i] = len[idx];
        fixed_slice.flags[i] = r2 ? r2[idx]*FQZ_FREAD2 : 0;
        fixed_slice.flags[i] |= sel ? (sel[idx]<<16) : 0;
    }

    return &fixed_slice;
}

static uint64_t manual_strats[10] = {0};
static int manual_nstrat = 0;

/*
 * Manually specified strategies held in global manual_strats[].
 */
static inline
int fqz_manual_parameters(fqz_gparams *gp,
                          fqz_slice *s,
                          unsigned char *in,
                          size_t in_size) {
    int i, p;
    int dsqr[] = {
        0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
    };

    gp->vers = FQZ_VERS;
    gp->nparam = manual_nstrat;
    gp->gflags = GFLAG_MULTI_PARAM | GFLAG_HAVE_STAB;
    for (i = 0; i < 256; i++)
        gp->stab[i] = 0;

    // Fill these out later
    gp->max_sel = 0;
    gp->max_sym = 0;
    gp->p = malloc(gp->nparam * sizeof(*gp->p));

    for (p = 0; p < gp->nparam; p++) {
        fqz_param *pm = &gp->p[p];
        uint64_t st = manual_strats[p];

        pm->do_qa  = st & 15; st >>= 4;
        pm->do_r2  = st & 15; st >>= 4;
        pm->dloc   = st & 15; st >>= 4;
        pm->ploc   = st & 15; st >>= 4;
        pm->sloc   = st & 15; st >>= 4;
        pm->qloc   = st & 15; st >>= 4;
        pm->dshift = st & 15; st >>= 4;
        pm->dbits  = st & 15; st >>= 4;
        pm->pshift = st & 15; st >>= 4;
        pm->pbits  = st & 15; st >>= 4;
        pm->qshift = st & 15; st >>= 4;
        pm->qbits  = st & 15; st >>= 4;

        // Gather some stats, as per qual_stats func.
        // r in rec count.
        // i = index to in[]
        // j = index within this rec
        uint32_t qhist[256] = {0};

        // qual stats for seqs using this parameter only
        fqz_qual_stats(s, in, in_size, pm, qhist, p);
        int max_sel = pm->max_sel;

        // Update max_sel running total. Eg with 4 sub-params:
        //
        // sel    param no.   => new
        // 0      0              0
        // 0/1    1              1,2
        // 0/1    2              3,4
        // 0      3              5
        for (i = gp->max_sel; i < gp->max_sel + max_sel+1; i++)
            gp->stab[i] = p;
        gp->max_sel += max_sel+1;

        pm->fixed_len = pm->fixed_len > 0;
        pm->use_qtab = 0;  // unused by current encoder
        pm->store_qmap = pm->nsym <= 8;

        // Adjust parameters based on quality stats.
        // FIXME:  dup from fqz_pick_parameters.
        for (i = 0; i < sizeof(dsqr)/sizeof(*dsqr); i++)
            if (dsqr[i] > (1<<pm->dbits)-1)
                dsqr[i] = (1<<pm->dbits)-1;

        if (pm->store_qmap) {
            int j;
            for (i = j = 0; i < 256; i++)
                if (qhist[i])
                    pm->qmap[i] = j++;
                else
                    pm->qmap[i] = INT_MAX;
            pm->max_sym = pm->nsym;
        } else {
            pm->nsym = 255;
            for (i = 0; i < 256; i++)
                pm->qmap[i] = i;
        }
        if (gp->max_sym < pm->max_sym)
            gp->max_sym = pm->max_sym;

        // Produce ptab from pshift.
        if (pm->qbits) {
            for (i = 0; i < 256; i++) {
                pm->qtab[i] = i; // 1:1

                // Alternative mappings:
                //qtab[i] = i > 30 ? MIN(max_sym,i)-15 : i/2;  // eg for 9827 BAM
            }

        }
        pm->qmask = (1<<pm->qbits)-1;

        if (pm->pbits) {
            for (i = 0; i < 1024; i++)
                pm->ptab[i] = MIN((1<<pm->pbits)-1, i>>pm->pshift);

            // Alternatively via analysis of quality distributions we
            // may select a bunch of positions that are special and
            // have a non-uniform ptab[].
            // Manual experimentation on a NovaSeq run saved 2.8% here.
        }

        if (pm->dbits) {
            for (i = 0; i < 256; i++)
                pm->dtab[i] = dsqr[MIN(sizeof(dsqr)/sizeof(*dsqr)-1, i>>pm->dshift)];
        }

        pm->use_ptab = (pm->pbits > 0);
        pm->use_dtab = (pm->dbits > 0);

        pm->pflags =
            (pm->use_qtab   ?PFLAG_HAVE_QTAB :0)|
            (pm->use_dtab   ?PFLAG_HAVE_DTAB :0)|
            (pm->use_ptab   ?PFLAG_HAVE_PTAB :0)|
            (pm->do_sel     ?PFLAG_DO_SEL    :0)|
            (pm->fixed_len  ?PFLAG_DO_LEN    :0)|
            (pm->do_dedup   ?PFLAG_DO_DEDUP  :0)|
            (pm->store_qmap ?PFLAG_HAVE_QMAP :0);
    }

    for (i = gp->max_sel; i < 256; i++)
        gp->stab[i] = gp->stab[gp->max_sel-1];

    return 0;
}

#define BS 1024*1024
static unsigned char *load(char *fn, size_t *lenp) {
    unsigned char *data = NULL;
    uint64_t dsize = 0;
    uint64_t dcurr = 0;
    signed int len;

    //build_rcp_freq();

#ifndef _O_BINARY
#define _O_BINARY 0
#endif

    int fd = open(fn, O_RDONLY | _O_BINARY);
    if (!fd) {
        perror(fn);
        return NULL;
    }

    do {
        if (dsize - dcurr < BS) {
            dsize = dsize ? dsize * 2 : BS;
            data = realloc(data, dsize);
        }

        len = read(fd, data + dcurr, BS);
        if (len > 0)
            dcurr += len;
    } while (len > 0);

    if (len == -1) {
        perror("read");
    }
    close(fd);

    *lenp = dcurr;
    return data;
}

#ifndef BLK_SIZE
#    define BLK_SIZE 300*1000000
#endif

int count_lines(unsigned char *in, size_t len) {
    size_t i;
    int lines = 0;

    for (i = 0; i < len; i++)
        if (in[i] == '\n')
            lines++;

    return lines;
}

// QUAL [is_read2 [selector]]
void parse_lines(unsigned char *in, size_t len,
                 int *rec_len, int *rec_r2, int *rec_sel, size_t *new_len) {
    size_t i, j, start;
    int rec = 0;

    for (start = i = j = 0; i < len; i++) {
        if (in[i] == '\n' || in[i] == ' ' || in[i] == '\t') {
            rec_len[rec] = i-start;

            // Read2 marker
            while (i < len && in[i] != '\n' && isspace(in[i]))
                i++;

            if (in[i] != '\n')
                rec_r2[rec] = atoi((char *)&in[i]);
            else
                rec_r2[rec] = 0;

            while (i < len && !isspace(in[i]))
                i++;

            // selector
            while (i < len && in[i] != '\n' && isspace(in[i]))
                i++;

            if (in[i] != '\n')
                rec_sel[rec] = atoi((char *)&in[i]);
            else
                rec_sel[rec] = 0;

            while (i < len && in[i] != '\n')
                i++;

            start = i+1;
            rec++;
        } else {
            in[j++] = in[i]-33; // ASCII phred to qual
        }
    }
    *new_len = j;
}

int main(int argc, char **argv) {
    unsigned char *in, *out;
    size_t in_len, out_len;
    int decomp = 0, vers = 4;  // CRAM version 4.0 (4) or 3.1 (3)
    int strat = 0, raw = 0;
    fqz_gparams *gp = NULL, gp_local;
    uint32_t blk_size = BLK_SIZE; // MAX

#ifdef _WIN32
        _setmode(_fileno(stdin),  _O_BINARY);
        _setmode(_fileno(stdout), _O_BINARY);
#endif

    extern char *optarg;
    extern int optind;
    int opt;

    while ((opt = getopt(argc, argv, "ds:s:b:rx:")) != -1) {
        switch (opt) {
        case 'd':
            decomp = 1;
            break;

        case 'b':
            blk_size = atoi(optarg);
            if (blk_size > BLK_SIZE)
                blk_size = BLK_SIZE;
            break;

        case 's':
            strat = atoi(optarg);
            break;

        case 'x': {
            // Hex digits are:
            // qbits  qshift
            // pbits  pshift
            // dbits  dshift
            // qloc   sloc
            // ploc   dloc
            // do_r2  do_qavg
            //
            // Examples: -x 0x5570000d6e14 q40+dir =  3473340
            //           -x 0x8252120e8d04 q4      =  724989
            uint64_t x = strtol(optarg, NULL, 0);
            manual_strats[manual_nstrat++] = x;

            gp = &gp_local;
            break;
        }
        case 'r':
            raw = 1;
            break;
        }
    }

    in = load(optind < argc ? argv[optind] : "/dev/stdin", &in_len);
    if (!in)
        exit(1);

    if (raw)
        blk_size = in_len;

    // Block based, for arbitrary sizes of input
    if (decomp) {
        unsigned char *in2 = in;
        while (in_len > 0) {
            // Read sizes as 32-bit
            size_t in2_len, out_len;
            if (raw) {
                uint32_t u32;
                var_get_u32(in2, in2+in_len, &u32);
                out_len = u32;
                in2_len = in_len;
            } else {
                out_len = *(uint32_t *)in2;  in2 += 4;
                in2_len = *(uint32_t *)in2;  in2 += 4;
            }

            fprintf(stderr, "out_len %ld, in_len %ld\n", (long)out_len, (long)in2_len);

            int *lengths = malloc(MAX_REC * sizeof(int));
            out = (unsigned char *)fqz_decompress((char *)in2, in_len-(raw?0:8), &out_len, lengths, MAX_REC);
            if (!out) {
                fprintf(stderr, "Failed to decompress\n");
                return 1;
            }

            // Convert from binary back to ASCII with newlines
            int i = 0, j = 0;
            while (j < out_len) {
                int k;
                char seq[MAX_SEQ];
                for (k = 0; k < lengths[i]; k++)
                    seq[k] = out[j+k]+33;
                seq[k] = 0;
                puts(seq);
                j += lengths[i++];
            }
            free(out);
            in2 += in2_len;
            in_len -= in2_len+(raw?0:8);

            free(lengths);

            break; // One cycle only until we fix blocking to be \n based
        }
    } else {
        // Convert from ASCII newline separated file to binary block.
        // We return an array of line lengths and optionally param selectors.
        int nlines = count_lines(in, in_len);
        fprintf(stderr, "nlines=%d\n", nlines);
        int *rec_len = calloc(nlines, sizeof(*rec_len));
        int *rec_r2  = calloc(nlines, sizeof(*rec_r2));
        int *rec_sel = calloc(nlines, sizeof(*rec_sel));
        parse_lines(in, in_len, rec_len, rec_r2, rec_sel, &in_len);

        unsigned char *in2 = in;
        long t_out = 0;
        out = NULL;
        while (in_len > 0) {
            // FIXME: blk_size no longer working in test.  One cycle only!
            size_t in2_len = in_len <= blk_size ? in_len : blk_size;
            fqz_slice *s = fake_slice(in2_len, rec_len, rec_r2, rec_sel, nlines);
            if (gp == &gp_local)
                if (fqz_manual_parameters(gp, s, in2, in2_len) < 0)
                    return 1;
            out = (unsigned char *)fqz_compress(vers, s, (char *)in2, in2_len, &out_len, strat, gp);

            // Write out 32-bit sizes.
            if (!raw) {
                uint32_t u32;
                u32 = in2_len; if (write(1, &u32, 4) != 4) return 1;
                u32 = out_len; if (write(1, &u32, 4) != 4) return 1;
            }
            if (write(1, out, out_len) < 0) return 1;
            in_len -= in2_len;
            in2 += in2_len;
            t_out += out_len + (raw?0:8);

            break; // One cycle only until we fix blocking to be \n based
        }
        free(out);
        free(rec_len);
        free(rec_r2);
        free(rec_sel);
        fprintf(stderr, "Total output = %ld\n", t_out);
    }

    free(in);

    return 0;
}
