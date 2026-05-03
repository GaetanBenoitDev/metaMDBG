# Copyright (c) 2019 Genome Research Ltd.
# Author(s): James Bonfield
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# 
#    1. Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
# 
#    2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
# 
#    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
#    Institute nor the names of its contributors may be used to endorse
#    or promote products derived from this software without specific
#    prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
# LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 

# Standalone test programs
noinst_PROGRAMS = rans4x16pr tokenise_name3 arith_dynamic rans4x8 rans4x16pr fqzcomp_qual varint entropy

LDADD = $(top_builddir)/htscodecs/libhtscodecs.la
AM_CPPFLAGS = -I$(top_srcdir)

fqzcomp_qual_SOURCES = fqzcomp_qual_test.c
rans4x8_SOURCES = rANS_static_test.c
rans4x16pr_SOURCES = rANS_static4x16pr_test.c
arith_dynamic_SOURCES = arith_dynamic_test.c
tokenise_name3_SOURCES = tokenise_name3_test.c
varint_SOURCES = varint_test.c
entropy_SOURCES = entropy.c

test_scripts = \
	entropy.test \
	rans4x8.test \
	rans4x16.test \
	arith.test \
	tok3.test \
	fqzcomp.test

TESTS = $(test_scripts) \
	varint

EXTRA_DIST = $(test_scripts) dat names

test_outdir = test.out

distclean-local:
	-rm -rf $(test_outdir)

# Fuzz testing programs, for local use.
#
# Best configure the package in a build subdirectory
# ../configure --disable-shared CFLAGS="-g -gdwarf-2 -Wall -O3 -fsanitize=address,undefined"
# Then "make fuzz" to recompile these tools.
#
# We can get non-fuzzing code for valgrind use, having found a bug, with e.g.
# clang70 -I. -I../../tests -I..  -I../.. -DFUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION -Wall -g -Wall -gdwarf-2 ../../tests/tokenise_name3_fuzz.c -DNOFUZZ -lbz2 -lpthread


# Our fuzz tools are built against a version of the library compiled
# with -DFUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION.
# This limits the maximum size of some data blocks to speed up
# fuzz testing, avoiding small inputs uncompressing to huge
# outputs.

fuzz: $(EXTRA_PROGRAMS)

# Additional flags
fuzzer_cflags  = -fsanitize=fuzzer -DFUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
fuzzer_ldflags = -fsanitize=fuzzer
fuzzer_ldadd   = $(top_builddir)/htscodecs/libcodecsfuzz.a \
	$(top_builddir)/htscodecs/libcodecsfuzz_sse4.a \
	$(top_builddir)/htscodecs/libcodecsfuzz_avx2.a \
	$(top_builddir)/htscodecs/libcodecsfuzz_avx512.a

EXTRA_PROGRAMS = \
	rans4x8_fuzz \
	rans4x16pr_fuzz \
	arith_dynamic_fuzz \
	tokenise_name3_fuzz \
	tokenise_name3_fuzzrt \
	fqzcomp_qual_fuzz \
	fqzcomp_qual_fuzzrt \
	entropy_fuzz

rans4x8_fuzz_SOURCES = rANS_static_fuzz.c
rans4x8_fuzz_CFLAGS  = $(fuzzer_cflags)
rans4x8_fuzz_LDFLAGS = $(fuzzer_ldflags)
rans4x8_fuzz_LDADD   = $(fuzzer_ldadd)

rans4x16pr_fuzz_SOURCES = rANS_static4x16pr_fuzz.c
rans4x16pr_fuzz_CFLAGS  = $(fuzzer_cflags)
rans4x16pr_fuzz_LDFLAGS = $(fuzzer_ldflags)
rans4x16pr_fuzz_LDADD   = $(fuzzer_ldadd)

arith_dynamic_fuzz_SOURCES = arith_dynamic_fuzz.c
arith_dynamic_fuzz_CFLAGS  = $(fuzzer_cflags)
arith_dynamic_fuzz_LDFLAGS = $(fuzzer_ldflags)
arith_dynamic_fuzz_LDADD   = $(fuzzer_ldadd)

tokenise_name3_fuzz_SOURCES = tokenise_name3_fuzz.c
tokenise_name3_fuzz_CFLAGS  = $(fuzzer_cflags)
tokenise_name3_fuzz_LDFLAGS = $(fuzzer_ldflags)
tokenise_name3_fuzz_LDADD   = $(fuzzer_ldadd)

tokenise_name3_fuzzrt_SOURCES = tokenise_name3_fuzzrt.c
tokenise_name3_fuzzrt_CFLAGS  = $(fuzzer_cflags)
tokenise_name3_fuzzrt_LDFLAGS = $(fuzzer_ldflags)
tokenise_name3_fuzzrt_LDADD   = $(fuzzer_ldadd)

fqzcomp_qual_fuzz_SOURCES = fqzcomp_qual_fuzz.c
fqzcomp_qual_fuzz_CFLAGS  = $(fuzzer_cflags)
fqzcomp_qual_fuzz_LDFLAGS = $(fuzzer_ldflags)
fqzcomp_qual_fuzz_LDADD   = $(fuzzer_ldadd)

entropy_fuzz_SOURCES = entropy_fuzz.c
entropy_fuzz_CFLAGS  = $(fuzzer_cflags)
entropy_fuzz_LDFLAGS = $(fuzzer_ldflags)
entropy_fuzz_LDADD   = $(fuzzer_ldadd)

fqzcomp_qual_fuzzrt_SOURCES = fqzcomp_qual_fuzzrt.c
fqzcomp_qual_fuzzrt_CFLAGS  = $(fuzzer_cflags)
fqzcomp_qual_fuzzrt_LDFLAGS = $(fuzzer_ldflags)
fqzcomp_qual_fuzzrt_LDADD   = $(fuzzer_ldadd)
