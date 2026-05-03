/* Tests for hts_codecs */
/*
 * Copyright (c) 2022 Genome Research Ltd.
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

/*
 * This test aims to test all entropy codecs on an input file.
 * This therefore validates the pthead_once memory allocations to ensure
 * there are not unforseen initialisation interactions.
 *
 * We repeatedly compress and decompress a single input file,
 * validating the result.
 */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <fcntl.h>
#include <sys/time.h>
#ifndef _WIN32
#include <sys/resource.h>
#include <pthread.h>
#endif

#include "htscodecs/arith_dynamic.h"
#include "htscodecs/rANS_static.h"
#include "htscodecs/rANS_static4x16.h"

#ifndef BLK_SIZE
#   define BLK_SIZE 1024*1024
#endif

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
    FILE *infp = stdin;
    int result = EXIT_SUCCESS;
#ifdef _WIN32
        _setmode(_fileno(stdin),  _O_BINARY);
        _setmode(_fileno(stdout), _O_BINARY);
#endif

    extern void force_sw32_decoder(void);
    extern void rans_disable_avx512(void);
    extern void rans_disable_avx2(void);

    int benchmark = 0;
    while (argc > 1 && strcmp(argv[1], "-b") == 0) {
        benchmark++;
        argc--;
        argv++;
    }

    if (argc > 1) {
        if (!(infp = fopen(argv[1], "rb"))) {
            perror(argv[1]);
            return 1;
        }
    }

#ifndef _WIN32
    // Specify an extra small stack, eg as in Alpine linux threads.
    // This checks we're not accidentally needing high stack usage.
    struct rlimit r = {64*1024, 64*1024};
    setrlimit(RLIMIT_STACK, &r);
#endif

    uint32_t in_size, csize, usize;
    unsigned char *in = load(infp, &in_size);
    int order_a[] = {0,1,                            // r4x8
                     64,65, 128,129, 192,193,        // r4x16, arith
                     4,5, 68,69, 132,133, 196,197,   // r4x16 SIMD
                     };
    char *codec[] = {"r4x8", "r4x16", "r32x16", "arith"};
    int i, j;
    for (i = 0; i < sizeof(order_a) / sizeof(*order_a); i++) {
        int order = order_a[i];
        uint8_t *comp, *uncomp;
        for (j = 0; j < 4; j++) {
            int chigh = 4, clow = 0, c;
            uint8_t *comp0 = NULL;
            uint32_t csize0 = 0;
            for (c = 0; c < 4; c+=(j==2)?1:4) {
                struct timeval tv1, tv2, tv3, tv4;

                // Test combinations of SIMD implementations
                uint32_t chex = benchmark
                    ? (clow<<8) | clow
                    : (clow<<8) | chigh;
                clow  = 1<<c;
                chigh >>= 1;
                rans_set_cpu(chex);

                int bloop = benchmark;
            bloop:
                // encode
                gettimeofday(&tv1, NULL);
                switch (j) {
                case 0: // r4x8
                    if (i >= 2) continue;
                    comp = rans_compress(in, in_size, &csize, order);
                    break;

                case 1: // r4x16
                    if (i >= 8) continue;
                    comp = rans_compress_4x16(in, in_size, &csize, order);
                    break;

                case 2: // r32x16
                    if (i < 8) continue;
                    comp = rans_compress_4x16(in, in_size, &csize, order);
                    break;

                case 3: // arith
                    if (i >= 8) continue;
                    comp = arith_compress(in, in_size, &csize, order);
                    break;
                }
                gettimeofday(&tv2, NULL);

                if (j == 2)
                    printf("%10s-o%d-c%04x\t", codec[j], order, chex);
                else
                    printf("%10s-o%d      \t", codec[j], order);
                printf("%10d uncomp, %10d comp", in_size, csize);

                if (comp == NULL) {
                    printf("\tFAIL (comp)\n");
                    result = EXIT_FAILURE;
                    continue;
                }

                if (comp0) {
                    if (csize != csize0 || memcmp(comp, comp0, csize) != 0) {
                        printf("\tFAIL (comp cmp)\n");
                        result = EXIT_FAILURE;
                    }
                } else {
                    csize0 = csize;
                    comp0 = comp;
                }

                // decode
                gettimeofday(&tv3, NULL);
                switch (j) {
                case 0: // r4x8
                    if (i >= 2) continue;
                    uncomp = rans_uncompress(comp, csize, &usize);
                    break;

                case 1: // r4x16
                    if (i >= 8) continue;
                    uncomp = rans_uncompress_4x16(comp, csize, &usize);
                    break;

                case 2: // r32x16
                    if (i < 8) continue;
                    uncomp = rans_uncompress_4x16(comp, csize, &usize);
                    break;

                case 3: // arith
                    if (i >= 8) continue;
                    uncomp = arith_uncompress(comp, csize, &usize);
                    break;
                }
                gettimeofday(&tv4, NULL);

                if (usize != in_size || uncomp == NULL || memcmp(in, uncomp, usize) != 0) {
                    printf("\tFAIL\n");
                    result = EXIT_FAILURE;
                } else if (benchmark) {
                    printf(" %6.1f enc MB/s %6.1f dec MB/s\n",
                           (double)usize /
                             ((long)(tv2.tv_sec - tv1.tv_sec)*1000000 +
                              tv2.tv_usec - tv1.tv_usec),
                           (double)usize /
                             ((long)(tv4.tv_sec - tv3.tv_sec)*1000000 +
                              tv4.tv_usec - tv3.tv_usec));
                } else {
                    printf("\tpass\n");
                }

                if (comp != comp0)
                    free(comp);
                free(uncomp);
                if (--bloop > 0) goto bloop;
            }
            free(comp0);
        }
        printf("\n");
    }

    free(in);

    if (result != EXIT_SUCCESS)
        return result;

#ifndef _WIN32
    // We wouldn't normally exit this way, but we explicitly call it to
    // check htscodecs_tls_free_all has no leaks.  Note that this will
    // cause the program to return an exit status of zero.
    pthread_exit(NULL);
#endif

    return 0;
}
