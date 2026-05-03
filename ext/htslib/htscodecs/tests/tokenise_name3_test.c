/*
 * Copyright (c) 2016-2020 Genome Research Ltd.
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
 *    Institute nor the names of its contributors may be used to endorse
 *    or promote products derived from this software without specific
 *    prior written permission.
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
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <ctype.h>
#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <fcntl.h>
#include <errno.h>
#include <time.h>

#include "htscodecs/tokenise_name3.h"

//-----------------------------------------------------------------------------
// main() implementation for testing

// Large enough for whole file for now.
#ifndef BLK_SIZE
#define BLK_SIZE 1*1024*1024
#endif
static char *blk;

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

static int encode(int argc, char **argv) {
    FILE *fp;
    int len, level = 9;
    int use_arith = 0;
    int raw = 0;

    while (argc > 1 && argv[1][0] == '-') {
        if (strcmp(argv[1], "-r") == 0) {
            raw = 1;
            argc--;
            argv++;
        }

        else if (argv[1][1] >= '0' && argv[1][1] <= '9') {
            level = atoi(argv[1]+1);
            if (level > 10) {
                level -= 10;
                use_arith = 1;
            }
            argc--;
            argv++;
        }

        else
            exit(1);
    }

    if (argc > 1) {
        fp = fopen(argv[1], "r");
        if (!fp) {
            perror(argv[1]);
            return 1;
        }
    } else {
        fp = stdin;
    }

    if (raw) {
        // One naked / raw block, to match the specification
        uint32_t in_len;
        int out_len;
        unsigned char *in = load(fp, &in_len), *out;
        if (!in) exit(1);
        out = tok3_encode_names((char *)in, in_len, level, use_arith,
                                &out_len, NULL);
        if (!out || write(1, out, out_len) < out_len) exit(1);   // encoded data
        free(in);
        free(out);
    } else {
        // Block based, to permit arbitrarily large files for benchmarking
        int blk_offset = 0;
        for (;;) {
            int last_start = 0;

            len = fread(blk+blk_offset, 1, BLK_SIZE-blk_offset, fp);
            if (len <= 0)
                break;
            len += blk_offset;

            int out_len;
            uint8_t *out = tok3_encode_names(blk, len, level, use_arith,
                                             &out_len, &last_start);
            if (!out) {
                fprintf(stderr, "Couldn't encode names\n");
                exit(1);
            }
            if (write(1, &out_len, 4) < 4) exit(1);
            if (write(1, out, out_len) < out_len) exit(1);   // encoded data
            free(out);

            if (len > last_start)
                memmove(blk, &blk[last_start], len - last_start);
            blk_offset = len - last_start;
        }
    }

    if (fclose(fp) < 0) {
        perror("closing file");
        return 1;
    }

    return 0;
}

static int decode(int argc, char **argv) {
    uint32_t in_sz, out_sz;
    int raw = 0;

    if (argc > 1 && strcmp(argv[1], "-r") == 0) {
        raw = 1;
        argc--;
        argv++;
    }

    if (raw) {
        // One naked / raw block, to match the specification
        uint32_t in_len;
        unsigned char *in = load(stdin, &in_len), *out;
        if (!in) exit(1);

        if ((out = tok3_decode_names(in, in_len, &out_sz)) == NULL)
            exit(1);
        if (write(1, out, out_sz) != out_sz)
            exit(1);

        free(in);
        free(out);
    } else {
        while (fread(&in_sz, 1, 4, stdin) == 4) {
            uint8_t *in = malloc(in_sz), *out;
            if (!in)
                return -1;

            if (fread(in, 1, in_sz, stdin) != in_sz) {
                free(in);
                return -1;
            }

            if ((out = tok3_decode_names(in, in_sz, &out_sz)) == NULL) {
                free(in);
                return -1;
            }

            if (write(1, out, out_sz) < out_sz) exit(1);

            free(in);
            free(out);
        }
    }

    return 0;
}

int main(int argc, char **argv) {
#ifdef _WIN32
        _setmode(_fileno(stdin),  _O_BINARY);
        _setmode(_fileno(stdout), _O_BINARY);
#endif

    // temporary fix for decoder, which needs more space
    blk = malloc(BLK_SIZE*2);

    int ret;

    if (argc > 1 && strcmp(argv[1], "-d") == 0)
        ret = decode(argc-1, argv+1);
    else
        ret = encode(argc, argv);

    free(blk);
    return ret;
}
