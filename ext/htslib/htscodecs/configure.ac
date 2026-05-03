dnl Process this file with autoconf to produce a configure script.
AC_INIT(htscodecs, 1.6.6)

# Some functions benefit from -O3 optimisation, so if the user didn't
# explicitly set any compiler flags, we'll plump for O3.
_=${CFLAGS="-g -O3"}

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE
AM_EXTRA_RECURSIVE_TARGETS([fuzz])

# For my own benefit, here's a summary of how versioning should work.
#
# Library versioning, not to be confused with the release number:
# Update VERS_CURRENT when an interface changes in an incompatible manner.
#     eg. change of struct, removing a function, changing a prototype.
#     Also increment this when we add new functions in a backwards compatible
#     manner, but update VERS_AGE too (see below).
#
# Update VERS_REVISION whenever VERS_CURRENT hasn't changed, but the source
#     code has. (Otherwise reset revision to 0).
#
# Update VERS_AGE if new interfaces have been added, causing VERS_CURRENT to
#        be incremented too. The idea is that if the previous values are
#        VERS_CURRENT=2, VERS_AGE=0 then after adding a new function we
#        get VERS_CURRENT=3, VERS_AGE=1 (ie ABI compat with vers 2 or 3).
#
# Also see:
# http://www.gnu.org/software/libtool/manual/libtool.html#Updating-version-info

# Example: We have htscodecs 1.11.6 implementing ABI vers 0, rev 6, age 0.
# Our library has an soname of libhtscodecs.so.0 and files as follows:
#       (SONAME      libhtscodecs.so.0)
# 	libhtscodecs.so -> libhtscodecs.so.0.0.6
#       libhtscodecs.so.0 -> libhtscodecs.so.0.0.6
#       libhtscodecs.so.0.0.6
#
# We then release htscodecs 1.12.0 with modified structs, so it's incompatible.
# We bump to ABI vers 1, rev 0, age 0 - giving:
#       (SONAME      libhtscodecs.so.1)
# 	libhtscodecs.so -> libhtscodecs.so.1.0.0
#       libhtscodecs.so.1 -> libhtscodecs.so.1.0.0
#       libhtscodecs.so.1.0.0
#
# If in 1.12.1 to 1.12.5 we then perform several successive minor tweaks (bug
# fix, optimisation, etc) that do not change the ABI at all, we may
# end up with ABI vers 1, rev 5, age 0:
#       (SONAME      libhtscodecs.so.1)
# 	libhtscodecs.so -> libhtscodecs.so.1.0.5
#       libhtscodecs.so.1 -> libhtscodecs.so.1.0.5
#       libhtscodecs.so.1.0.5
#
# And assume in 1.13.0 we add an extra public function. Code linked
# against 1.12.x still works (aka library version 1) meaning our
# soname doesn't change, but we do update ABI version: vers 2, rev 0, age 1.
#       (SONAME      libhtscodecs.so.1)
# 	libhtscodecs.so -> libhtscodecs.so.1.1.0
#       libhtscodecs.so.1 -> libhtscodecs.so.1.1.0
#       libhtscodecs.so.1.1.0

VERS_CURRENT=3
VERS_REVISION=11
VERS_AGE=1
AC_SUBST(VERS_CURRENT)
AC_SUBST(VERS_REVISION)
AC_SUBST(VERS_AGE)

dnl Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
LT_INIT

# # by Marcelo Magallon <mmagallo@efis.ucr.ac.cr>
# # Turn around -rpath problem with libtool 1.0c
# # This define should be improbable enough to not conflict with anything.
# echo host=$host
# case ${host} in
#   *-linux-gnu | *-apple-darwin*)
#     AC_MSG_RESULT([Fixing libtool for -rpath problems.])
#     sed < libtool > libtool-2 \
#     's/^hardcode_libdir_flag_spec.*$'/'hardcode_libdir_flag_spec=""/'
#     mv libtool-2 libtool
#     chmod 755 libtool
#   ;;
# esac

dnl Turn on C compiler warnings.
VL_PROG_CC_WARNINGS

dnl Check for libraries
AC_ARG_ENABLE([bz2],
  [AS_HELP_STRING([--disable-bz2],
                  [omit support for BZ2 within Arith streams])],
  [], [enable_bz2=yes])

if test "$enable_bz2" != no; then
  bz2_devel=ok
  AC_CHECK_LIB([bz2], [BZ2_bzBuffToBuffCompress], [
	LIBS="-lbz2 $LIBS"
	AC_DEFINE([HAVE_LIBBZ2],1,[Define to 1 if you have the libbz2 library.])], [bz2_devel=missing])
  if test "$bz2_devel" != "ok"; then
    AC_MSG_ERROR([libbzip2 development files not found.

This is needed for some options of the arithmetic coder.  It can be disabled
with configure --disable-bz2, but be aware some data streams will become
unable to be decompressed and arith tests will fail.])
  fi
fi

dnl Check if __builtin_prefetch exists.
AC_CACHE_CHECK([for __builtin_prefetch], [ax_cv_builtin_prefetch],[
AC_LINK_IFELSE([AC_LANG_PROGRAM([], [__builtin_prefetch("")])],
     [ax_cv_builtin_prefetch=yes],
     [ax_cv_builtin_prefetch=no])
])
AS_IF([test "$ax_cv_builtin_prefetch" = "yes"],
        [AC_DEFINE(HAVE_BUILTIN_PREFETCH, 1,
           [Define to 1 if you have a __builtin_prefetch])],
	[])

dnl Check is cpuid works, needed by rANS_static4x16pr.c.
AC_CHECK_DECLS([__get_cpuid_max, __cpuid_count], [], [], [[#include <cpuid.h>]])

dnl AC_CHECK_LIB([lzma], [lzma_easy_buffer_encode], [
dnl 	LIBS="-llzma $LIBS"
dnl 	AC_DEFINE([HAVE_LIBLZMA],1,[Define to 1 if you have the liblzma library.])])
dnl 
dnl # Experimental
dnl AC_CHECK_LIB([bsc], [bsc_compress], [
dnl 	LIBS="-lbsc $LIBS"
dnl 	AC_DEFINE([HAVE_LIBBSC],1,[Define to 1 if you have the libbsc library.])])

dnl Check if we can use our SSE4.1 too.
dnl Our SSE4 codec uses SSE4.1, SSSE3 (shuffle) and POPCNT, so we check all 3
dnl together.  This helps Zig builds which don't work well if we test each
dnl individually.
HTS_CHECK_COMPILE_FLAGS_NEEDED([sse4.1], [-msse4.1 -mssse3 -mpopcnt], [AC_LANG_PROGRAM([[
	  #ifdef __x86_64__
	  #include "x86intrin.h"
	  #endif
	]],[[
	  #ifdef __x86_64__
	  __m128i a = _mm_set_epi32(1, 2, 3, 4), b = _mm_set_epi32(4, 3, 2, 1);
	  __m128i c = _mm_shuffle_epi8(_mm_max_epu32(a, b), b);
	  return _mm_popcnt_u32(*((char *) &c));
	  #endif
	]])], [
        MSSE4_1="$flags_needed"
	build_rans_sse4=yes
	AC_SUBST([MSSE4_1])
        AC_DEFINE([HAVE_SSE4_1],1,[Defined to 1 if rANS source using SSE4.1 can be compiled.])
        AC_DEFINE([HAVE_SSSE3],1,[Defined to 1 if rANS source using SSSE3 can be compiled.])
        AC_DEFINE([HAVE_POPCNT],1,[Defined to 1 if rANS source using popcnt can be compiled.])
])
AM_CONDITIONAL([RANS_32x16_SSE4],[test "$build_rans_sse4" = yes])

dnl Check if we can use our AVX2 implementations.
build_rans_avx2=no
HTS_CHECK_COMPILE_FLAGS_NEEDED([avx2], [-mavx2 -mpopcnt], [AC_LANG_PROGRAM([[
	  #ifdef __x86_64__
	  #include "x86intrin.h"
	  #endif
	]],[[
	  #ifdef __x86_64__
	  __m256i a = _mm256_set_epi32(1, 2, 3, 4, 5, 6, 7, 8);
	  __m256i b = _mm256_add_epi32(a, a);
	  long long c = _mm256_extract_epi64(b, 0);
	  return _mm_popcnt_u32((int)c);
	  #endif
	]])], [
        MAVX2="$flags_needed"
	build_rans_avx2=yes
	AC_SUBST([MAVX2])
        AC_DEFINE([HAVE_AVX2],1,[Defined to 1 if rANS source using AVX2 can be compiled.])
        AC_DEFINE([HAVE_POPCNT],1,[Defined to 1 if rANS source using popcnt can be compiled.])
])
AM_CONDITIONAL([RANS_32x16_AVX2],[test "$build_rans_avx2" = yes])

dnl Check also if we have AVX512.
build_rans_avx512=no
HTS_CHECK_COMPILE_FLAGS_NEEDED([avx512f], [-mavx512f -mpopcnt], [AC_LANG_PROGRAM([[
	  #ifdef __x86_64__
	  #include "x86intrin.h"
	  #endif
	]],[[
	  #ifdef __x86_64__
	  __m512i a = _mm512_set1_epi32(1);
	  __m512i b = _mm512_add_epi32(a, a);
          __m256i c = _mm512_castsi512_si256(b);
	  __m256i d = _mm512_extracti64x4_epi64(a, 1);
	  return _mm_popcnt_u32(*((char *) &c)) + (*(char *) &d);
	  #endif
	]])], [
        MAVX512="$flags_needed"
	build_rans_avx512=yes
	AC_SUBST([MAVX512])
        AC_DEFINE([HAVE_AVX512],1,[Defined to 1 if rANS source using AVX512F can be compiled.])
])
        AC_DEFINE([HAVE_POPCNT],1,[Defined to 1 if rANS source using popcnt can be compiled.])
AM_CONDITIONAL([RANS_32x16_AVX512],[test "$build_rans_avx512" = yes])

AC_SUBST([HTSCODECS_SIMD_SRC])

dnl Checks for header files.
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS(fcntl.h limits.h unistd.h malloc.h)
AC_CHECK_HEADER(zlib.h)

dnl Checks for library functions.
AX_SEARCH_LIBS_REV([pthread_join], [thr pthread])

AC_CONFIG_FILES([Makefile htscodecs/Makefile tests/Makefile])
AC_OUTPUT
