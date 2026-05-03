/* Tests for CRAM-3.0 rANS codec */
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

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <fcntl.h>
#include <sys/time.h>

#include "htscodecs/rANS_static.h"

#ifndef BLK_SIZE
// Divisible by 4 for X4
#  define BLK_SIZE 1039*251*4
#endif

// Room to allow for expanded BLK_SIZE on worst case compression.
#define BLK_SIZE2 (1.05*BLK_SIZE)

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

/*-----------------------------------------------------------------------------
 * Main.
 *
 * This is a simple command line tool for testing order-0 and order-1
 * compression using the rANS codec. Simply compile with
 *
 * gcc -DTEST_MAIN -O3 -I. cram/rANS_static.c -o cram/rANS_static
 *
 * Usage: cram/rANS_static -o0 < file    > file.o0
 *        cram/rANS_static -d  < file.o0 > file2
 *
 *        cram/rANS_static -o1 < file    > file.o1
 *        cram/rANS_static -d  < file.o1 > file2
 */
int main(int argc, char **argv) {
    int opt, order = 0;
    unsigned char *in_buf = malloc(BLK_SIZE2+257*257*3);
    int decode = 0, test = 0;
    FILE *infp = stdin, *outfp = stdout;
    struct timeval tv1, tv2, tv3;
    size_t bytes = 0, raw = 0;

#ifdef _WIN32
        _setmode(_fileno(stdin),  _O_BINARY);
        _setmode(_fileno(stdout), _O_BINARY);
#endif

    extern char *optarg;
    extern int optind;

    while ((opt = getopt(argc, argv, "o:dtr")) != -1) {
        switch (opt) {
        case 'o':
            order = atoi(optarg);
            break;

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

    order = order ? 1 : 0; // Only support O(0) and O(1)

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
        blocks *b = NULL, *bc, *bu;
        int nb = 0, i;
        
        while ((len = fread(in_buf, 1, BLK_SIZE, infp)) != 0) {
            // inefficient, but it'll do for testing
            b = realloc(b, (nb+1)*sizeof(*b));
            b[nb].blk = malloc(len);
            b[nb].sz = len;
            memcpy(b[nb].blk, in_buf, len);
            nb++;
            in_sz += len;
        }

        int trials = 5;
        while (trials--) {
            bc = malloc(nb*sizeof(*bc));
            bu = malloc(nb*sizeof(*bu));

            gettimeofday(&tv1, NULL);

            out_sz = 0;
            for (i = 0; i < nb; i++) {
                bc[i].blk = rans_compress(b[i].blk, b[i].sz, &bc[i].sz, order);
                out_sz += 5 + bc[i].sz;
                bc[i].blk = realloc(bc[i].blk, bc[i].sz);
            }
        
            gettimeofday(&tv2, NULL);

            for (i = 0; i < nb; i++) {
                bu[i].blk = rans_uncompress(bc[i].blk, bc[i].sz, &bu[i].sz);
            }

            gettimeofday(&tv3, NULL);

            for (i = 0; i < nb; i++) {
                if (b[i].sz != bu[i].sz || memcmp(b[i].blk, bu[i].blk, b[i].sz))
                    fprintf(stderr, "Mismatch in block %d\n", i);
                free(bc[i].blk);
                free(bu[i].blk);
            }
            free(bc);
            free(bu);

            fprintf(stderr, "%5.1f MB/s enc, %5.1f MB/s dec\t %ld bytes -> %ld bytes\n",
                    (double)in_sz / ((long)(tv2.tv_sec - tv1.tv_sec)*1000000 +
                                     tv2.tv_usec - tv1.tv_usec),
                    (double)in_sz / ((long)(tv3.tv_sec - tv2.tv_sec)*1000000 +
                                     tv3.tv_usec - tv2.tv_usec),
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
            if (!(out = rans_uncompress(in, in_size, &out_size)))
                exit(1);

            fwrite(out, 1, out_size, outfp);
            bytes = out_size;
        } else {
            if (!(out = rans_compress(in, in_size, &out_size, order)))
                exit(1);

            fwrite(out, 1, out_size, outfp);
            bytes += in_size;
        }

        free(in);
        free(out);
    } else {
        if (decode) {
            // Only used in some test implementations of RC_GetFreq()
            //RC_init();
            //RC_init2();

            for (;;) {
                uint32_t in_size, out_size;
                unsigned char *out;

                order = fgetc(infp);
                if (4 != fread(&in_size, 1, 4, infp))
                    break;
                if (in_size != fread(in_buf, 1, in_size, infp)) {
                    fprintf(stderr, "Truncated input\n");
                    exit(1);
                }
                out = rans_uncompress(in_buf, in_size, &out_size);
                if (!out)
                    abort();

                fwrite(out, 1, out_size, outfp);
                free(out);

                bytes += out_size;
            }
        } else {
            int loop=0;
            for (;;loop++) {
                uint32_t in_size, out_size;
                unsigned char *out;

                in_size = fread(in_buf, 1, BLK_SIZE, infp);
                if (loop && in_size <= 0)
                    break;

                out = rans_compress(in_buf, in_size, &out_size,
                                    order && in_size >= 4);

                fputc(order && in_size >= 4, outfp);
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
