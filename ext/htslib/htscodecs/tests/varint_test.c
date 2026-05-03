/* Variable-length encoding tests */
/*
 * Copyright (c) 2020 Genome Research Ltd.
 * Author(s): Rob Davies
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

#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

#include "htscodecs/varint.h"

typedef struct unsigned_test {
    uint64_t val;
    int      len;
    uint8_t  encoding[12];
} unsigned_test;

typedef struct signed_test {
    int64_t  val;
    int      len;
    uint8_t  encoding[12];
} signed_test;

void dump_encoding(size_t sz, const uint8_t *buffer) {
    size_t byte;
    for (byte = 0; byte < sz; byte++) {
        printf("%s0x%02x", byte ? " " : "", buffer[byte]);
    }
}

void dump_unsigned(int bits, int is_enc,
                   uint64_t val, int sz, const uint8_t *buffer) {
    printf("%d-bit ", bits);
    if (is_enc) printf("0x%0*"PRIx64" => ", bits / 4, val);
    dump_encoding(sz, buffer);
    if (!is_enc) printf(" => 0x%0*"PRIx64, bits / 4, val);
    printf("\n");
}

void dump_signed(int bits, int is_enc,
                 int64_t val, int sz, const uint8_t *buffer) {
    printf("%d-bit ", bits);
    if (is_enc) printf("%"PRId64" => ", val);
    dump_encoding(sz, buffer);
    if (!is_enc) printf(" => %"PRId64, val);
    printf("\n");
}

int check_put_unsigned(int bits, const unsigned_test *t,
                       int len, const uint8_t buffer[16], int verbose) {
    if (len != t->len || memcmp(t->encoding, buffer, len) != 0) {
        printf("var_put_u%d failed:\nExpected ", bits);
        dump_unsigned(bits, 1, t->val, t->len, t->encoding);
        printf("Got      ");
        dump_unsigned(bits, 1, t->val, len, buffer);
        return 1;
    }
    if (verbose) {
        dump_unsigned(bits, 1, t->val, len, buffer);
    }
    return 0;
}

int check_get_unsigned(int bits, const unsigned_test *t,
                       int len, uint64_t val) {
    if (val == t->val && len == t->len)
        return 0;

    printf("var_get_u%d failed:\nExpected ", bits);
    dump_unsigned(bits, 0, t->val, t->len, t->encoding);
    printf("Got      ");
    dump_unsigned(bits, 0, val, len, t->encoding);
    return 1;
}

int check_put_signed(int bits, const signed_test *t,
                     int len, const uint8_t buffer[16], int verbose) {
    if (len != t->len || memcmp(t->encoding, buffer, len) != 0) {
        printf("var_put_s%d failed:\nExpected ", bits);
        dump_signed(bits, 1, t->val, t->len, t->encoding);
        printf("Got      ");
        dump_signed(bits, 1, t->val, len, buffer);
        return 1;
    }
    if (verbose) {
        dump_signed(bits, 1, t->val, len, buffer);
    }
    return 0;
}

int check_get_signed(int bits, const signed_test *t,
                     int len, int64_t val) {
    if (val == t->val && len == t->len)
        return 0;

    printf("var_get_s%d failed:\nExpected ", bits);
    dump_signed(bits, 0, t->val, t->len, t->encoding);
    printf("Got      ");
    dump_signed(bits, 0, val, len, t->encoding);
    return 1;
}

#define NELE(X) (sizeof(X)/sizeof(X[0]))

int test_unsigned(int verbose) {
    uint8_t buffer[16] = { 0 };
    uint8_t *endp = buffer + sizeof(buffer);
    uint32_t v32;
    uint64_t v64;
    size_t i;
    int len;
    int res = 0;
    unsigned_test tests32[] = {
        {          0U, 1, { 0x00 } },
        {          1U, 1, { 0x01 } },
        {       0x7fU, 1, { 0x7f } },
        {       0x80U, 2, { 0x81, 0x00 } },
        {     0x1234U, 2, { 0xa4, 0x34 } },
        {   0x123456U, 3, { 0xc8, 0xe8, 0x56 } },
        { 0x12345678U, 5, { 0x81, 0x91, 0xd1, 0xac, 0x78 } },
        { 0x80000000U, 5, { 0x88, 0x80, 0x80, 0x80, 0x00 } },
        { 0xffffffffU, 5, { 0x8f, 0xff, 0xff, 0xff, 0x7f } }
    };
    unsigned_test tests64[] = {
        { 0x100000000ULL, 5, { 0x90, 0x80, 0x80, 0x80, 0x00 } },
        { 0x123456789abcULL, 7, { 0x84, 0xc6, 0xc5, 0xb3, 0xe2, 0xb5, 0x3c} },
        { 0x123456789abcdef0ULL, 9,
          { 0x92, 0x9a, 0x95, 0xcf, 0x89, 0xd5, 0xf3, 0xbd, 0x70 } },
        { 0x8000000000000000ULL, 10,
          { 0x81, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x00 } },
        { 0xffffffffffffffffULL, 10,
          { 0x81, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f } } 
    };

    for (i = 0; i < NELE(tests32); i++) {
        memset(buffer, 0x55, sizeof(buffer));
        len = var_put_u32(buffer, endp, (uint32_t) tests32[i].val);
        res |= check_put_unsigned(32, &tests32[i], len, buffer, verbose);
        memset(buffer, 0x55, sizeof(buffer));
        len = var_put_u64(buffer, endp, tests32[i].val);
        res |= check_put_unsigned(64, &tests32[i], len, buffer, verbose);
        len = var_get_u32(tests32[i].encoding,
                          tests32[i].encoding + tests32[i].len,
                          &v32);
        res |= check_get_unsigned(32, &tests32[i], len, v32);
        len = var_get_u64(tests32[i].encoding,
                          tests32[i].encoding + tests32[i].len,
                          &v64);
        res |= check_get_unsigned(64, &tests32[i], len, v64);
    }

    for (i = 0; i < NELE(tests64); i++) {
        memset(buffer, 0x55, sizeof(buffer));
        len = var_put_u64(buffer, endp, tests64[i].val);
        res |= check_put_unsigned(64, &tests64[i], len, buffer, verbose);
        len = var_get_u64(tests64[i].encoding,
                          tests64[i].encoding + tests64[i].len,
                          &v64);
        res |= check_get_unsigned(64, &tests64[i], len, v64);
    }

    return res;
}

int test_signed(int verbose) {
    uint8_t  buffer[16] = { 0 };
    uint8_t *endp = buffer + sizeof(buffer);
    int32_t  v32;
    int64_t  v64;
    size_t   i;
    int      len;
    int      res = 0;
    signed_test tests32[] = {
        {   0, 1, { 0x00 }, },
        {  -1, 1, { 0x01 }, },
        {   1, 1, { 0x02 }, },
        { -63, 1, { 0x7d }, },
        {  63, 1, { 0x7e }, },
        { -64, 1, { 0x7f, } },
        {  64, 2, { 0x81, 0x00 } },
        { -65, 2, { 0x81, 0x01 } },
        {  65, 2, { 0x81, 0x02 } },
        { -12345678, 4, { 0x8b, 0xe3, 0x85, 0x1b } },
        {  12345678, 4, { 0x8b, 0xe3, 0x85, 0x1c } },
        { -2147483647, 5, { 0x8f, 0xff, 0xff, 0xff, 0x7d} },
        {  2147483647, 5, { 0x8f, 0xff, 0xff, 0xff, 0x7e} },
        { -2147483647-1, 5, { 0x8f, 0xff, 0xff, 0xff, 0x7f} },
    };

    signed_test tests64[] = {
        { 2147483648LL, 5, { 0x90, 0x80, 0x80, 0x80, 0x00 } },
        {  -1234567890123456LL, 8,
           { 0x84, 0xb1, 0xb5, 0xa7, 0xc8, 0xd5, 0xea, 0x7f } },
        {  1234567890123456LL, 8,
           { 0x84, 0xb1, 0xb5, 0xa7, 0xc8, 0xd5, 0xeb, 0x00 } },
        { -9223372036854775807LL, 10,
          { 0x81, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7d } },
        { 9223372036854775807LL, 10,
          { 0x81, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7e } },
        { -9223372036854775807LL - 1LL, 10,
          { 0x81, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f } },
    };

    for (i = 0; i < NELE(tests32); i++) {
        memset(buffer, 0x55, sizeof(buffer));
        len = var_put_s32(buffer, endp, (int32_t) tests32[i].val);
        res |= check_put_signed(32, &tests32[i], len, buffer, verbose);
        memset(buffer, 0x55, sizeof(buffer));
        len = var_put_s64(buffer, endp, tests32[i].val);
        res |= check_put_signed(64, &tests32[i], len, buffer, verbose);
        len = var_get_s32(tests32[i].encoding,
                          tests32[i].encoding + tests32[i].len,
                          &v32);
        res |= check_get_signed(32, &tests32[i], len, v32);
        len = var_get_s64(tests32[i].encoding,
                          tests32[i].encoding + tests32[i].len,
                          &v64);
        res |= check_get_signed(64, &tests32[i], len, v64);
    }

    for (i = 0; i < NELE(tests64); i++) {
        memset(buffer, 0x55, sizeof(buffer));
        len = var_put_s64(buffer, endp, tests64[i].val);
        res |= check_put_signed(64, &tests64[i], len, buffer, verbose);
        len = var_get_s64(tests64[i].encoding,
                          tests64[i].encoding + tests64[i].len,
                          &v64);
        res |= check_get_signed(64, &tests64[i], len, v64);
    }

    return res;
}


int main(int argc, char **argv) {
    int opt;
    int verbose = 0;
    int res = 0;

    while ((opt = getopt(argc, argv, "v")) != -1) {
        switch (opt) {
        case 'v':
            verbose++;
            break;
        default:
            fprintf(stderr, "Unknown option '%c'\n", opt);
            return EXIT_FAILURE;
        }
    }

    res |= test_unsigned(verbose);
    res |= test_signed(verbose);
    return res;
}
