/* Arithmetic coder tests */
/*
 * Copyright (c) 2019,2020 Genome Research Ltd.
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

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>
#include <limits.h>
#include <fcntl.h>

#include "htscodecs/arith_dynamic.h"

#ifndef BLK_SIZE
// Divisible by 4 for X4
#  define BLK_SIZE 1039*251*4
#endif

// Room to allow for expanded BLK_SIZE on worst case compression.
#define BLK_SIZE2 ((105LL*BLK_SIZE)/100)

static unsigned char *in_buf;

// Max 4GB
static unsigned char *load(FILE *infp, uint32_t *lenp) {
    unsigned char *data = NULL;
    uint32_t dsize = 0;
    uint32_t dcurr = 0;
    signed int len;

    do {
        if (dsize - dcurr < BLK_SIZE) {
            dsize = dsize ? dsize * 2 : BLK_SIZE;
            data = realloc(data, dsize);
        }

        len = fread(data + dcurr, 1, BLK_SIZE, infp);
        if (len > 0)
            dcurr += len;
    } while (len > 0);

    if (len == -1) {
        perror("fread");
    }

    *lenp = dcurr;
    return data;
}

int main(int argc, char **argv) {
    int opt, order = 0;
    int decode = 0, test = 0;
    FILE *infp = stdin, *outfp = stdout;
    struct timeval tv1, tv2, tv3, tv4;
    size_t bytes = 0, raw = 0;

    in_buf = malloc(BLK_SIZE2+257*257*3);

#ifdef _WIN32
        _setmode(_fileno(stdin),  _O_BINARY);
        _setmode(_fileno(stdout), _O_BINARY);
#endif

    extern char *optarg;
    extern int optind;

    while ((opt = getopt(argc, argv, "o:dtr")) != -1) {
        switch (opt) {
        case 'o': {
            char *optend;
            order = strtol(optarg, &optend, 0);
            if (*optend == '.')
                order += atoi(optend+1)<<8;
            break;
        }

        case 'd':
            decode = 1;
            break;
            
        case 't':
            test = 1;
            break;

        case 'r':
            raw = 1;
            break;
        }
    }

    //order = order ? 1 : 0; // Only support O(0) and O(1)

    if (optind < argc) {
        if (!(infp = fopen(argv[optind], "rb"))) {
            perror(argv[optind]);
            return 1;
        }
        optind++;
    }

    if (optind < argc) {
        if (!(outfp = fopen(argv[optind], "wb"))) {
            perror(argv[optind]);
            return 1;
        }
        optind++;
    }

    gettimeofday(&tv1, NULL);

    if (test) {
        size_t len, in_sz = 0, out_sz = 0;
        typedef struct {
            unsigned char *blk;
            uint32_t sz;
        } blocks;
        blocks *b = NULL, *bc = NULL, *bu = NULL;
        int nb = 0, i;
        
        while ((len = fread(in_buf, 1, BLK_SIZE, infp)) != 0) {
            // inefficient, but it'll do for testing
            b = realloc(b, (nb+1)*sizeof(*b));
            bu = realloc(bu, (nb+1)*sizeof(*bu));
            bc = realloc(bc, (nb+1)*sizeof(*bc));
            b[nb].blk = malloc(len);
            b[nb].sz = len;
            memcpy(b[nb].blk, in_buf, len);
            bc[nb].sz = arith_compress_bound(BLK_SIZE, order);
            bc[nb].blk = malloc(bc[nb].sz);
            bu[nb].sz = len;
            bu[nb].blk = malloc(BLK_SIZE);
            nb++;
            in_sz += len;
        }
        fprintf(stderr, "Testing %d blocks\n", nb);

#ifndef NTRIALS
#define NTRIALS 10
#endif
        int trials = NTRIALS;
        while (trials--) {
            // Warmup
            for (i = 0; i < nb; i++) memset(bc[i].blk, 0, bc[i].sz);

            gettimeofday(&tv1, NULL);

            out_sz = 0;
            for (i = 0; i < nb; i++) {
                unsigned int csz = bc[i].sz;
                bc[i].blk = arith_compress_to(b[i].blk, b[i].sz, bc[i].blk, &csz, order);
                assert(csz <= bc[i].sz);
                out_sz += 5 + csz;
            }

            gettimeofday(&tv2, NULL);
            
            // Warmup
            for (i = 0; i < nb; i++) memset(bu[i].blk, 0, BLK_SIZE);

            gettimeofday(&tv3, NULL);

            for (i = 0; i < nb; i++)
                bu[i].blk = arith_uncompress_to(bc[i].blk, bc[i].sz, bu[i].blk, &bu[i].sz);

            gettimeofday(&tv4, NULL);

            for (i = 0; i < nb; i++) {
                if (b[i].sz != bu[i].sz || memcmp(b[i].blk, bu[i].blk, b[i].sz))
                    fprintf(stderr, "Mismatch in block %d, sz %d/%d\n", i, b[i].sz, bu[i].sz);
                //free(bc[i].blk);
                //free(bu[i].blk);
            }

            fprintf(stderr, "%5.1f MB/s enc, %5.1f MB/s dec\t %ld bytes -> %ld bytes\n",
                    (double)in_sz / ((long)(tv2.tv_sec - tv1.tv_sec)*1000000 +
                                     tv2.tv_usec - tv1.tv_usec),
                    (double)in_sz / ((long)(tv4.tv_sec - tv3.tv_sec)*1000000 +
                                     tv4.tv_usec - tv3.tv_usec),
                    (long)in_sz, (long)out_sz);
        }

        exit(0);
        
    }

    if (raw) {
        // One naked / raw block, to match the specification
        uint32_t in_size, out_size;
        unsigned char *in = load(infp, &in_size), *out;
        if (!in) exit(1);

        if (decode) {
            if (!(out = arith_uncompress(in, in_size, &out_size)))
                exit(1);

            fwrite(out, 1, out_size, outfp);
            bytes = out_size;
        } else {
            if (!(out = arith_compress(in, in_size, &out_size, order)))
                exit(1);

            fwrite(out, 1, out_size, outfp);
            bytes += in_size;
        }

        free(in);
        free(out);
    } else {
        // Block based, to permit arbitrarily large data sets.
        if (decode) {
            for (;;) {
                uint32_t in_size, out_size;
                unsigned char *out;

                if (4 != fread(&in_size, 1, 4, infp))
                    break;
                if (in_size > BLK_SIZE)
                    exit(1);

                if (in_size != fread(in_buf, 1, in_size, infp)) {
                    fprintf(stderr, "Truncated input\n");
                    exit(1);
                }
                out = arith_uncompress(in_buf, in_size, &out_size);
                if (!out)
                    exit(1);

                fwrite(out, 1, out_size, outfp);
                fflush(outfp);
                free(out);

                bytes += out_size;
            }
        } else {
            int loop = 0;
            for (;;loop++) {
                uint32_t in_size, out_size;
                unsigned char *out;

                in_size = fread(in_buf, 1, BLK_SIZE, infp);
                if (loop && in_size <= 0)
                    break;

                if (in_size < 4)
                    order &= ~1;

                out = arith_compress(in_buf, in_size, &out_size, order);

                fwrite(&out_size, 1, 4, outfp);
                fwrite(out, 1, out_size, outfp);
                free(out);

                bytes += in_size;
            }
        }
    }

    gettimeofday(&tv2, NULL);

    fprintf(stderr, "Took %ld microseconds, %5.1f MB/s\n",
            (long)(tv2.tv_sec - tv1.tv_sec)*1000000 +
            tv2.tv_usec - tv1.tv_usec,
            (double)bytes / ((long)(tv2.tv_sec - tv1.tv_sec)*1000000 +
                             tv2.tv_usec - tv1.tv_usec));

    free(in_buf);
    return 0;
}
