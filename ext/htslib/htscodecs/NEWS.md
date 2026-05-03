Release 1.6.6: 11th January 2026
--------------------------------

This is a bug fix release.

Bug fixes

- Prevent use of unaligned 16-bit integers from memmove.
  (PR #145 Robert Davies)

- Add missing support for X_CAT in arithmetic codec.
  (PR #146; Issue #144 Michael Macias)


Release 1.6.5: 9th December 2025
--------------------------------

This is a bug fix release.

Bug fixes

- Add cpuid checks for XSAVE, OSXSAVE and AVX.  Corrects auto-detection of
  SIMD version on machines that have but disable specific CPU features.
  (PR #140 Robert Davies, fixes samtools/samtools#2256 Ran Fan).

- Avoid undefined behaviour by replacing literal copies with memcpy
  (PR #142 James Bonfield, fixes Issue #141 Vasudeva Easwara Sarma)


Release 1.6.4: 9th July 2025
----------------------------

This is primarily a bug fix release.

Fixes

- Fixed a minor thread data race in first call of the rans4x16 codec.

- Protect against SIMD rANS encoding on small or highly compressible data
  sets.  This could fail when combined with the RLE method where one
  sub-component was very small (<32 bytes) and the other was large.

Changes

- UUID4 based read names are now compressed better with the name tokeniser.
  This also slightly improves name compression of mixed data sets.


Release 1.6.3: 22nd May 2025
----------------------------

A tiny bug fix to 1.6.2 to fix a memory leak in rans_compress_to_4x16,
detected by htslib's CI system. (#136)


Release 1.6.2: 22nd May 2025
----------------------------

This release has minor bug fixes and some continuous integration test
improvements.

Bug fixes

- Improved check of out_size in rans4x16 and arithmetic coder, plus
  better memory freeing on error. (#127]

- [CI] Bump FreeBSD release used to 14.2 and Ubuntu to 24.04 (#129 jkb,
  #133 from John Marshall).

- [CI] Remove GitHub workflow shell override. (#133, John Marshall)

- [JavaScript] Correct arithmetic coder, (Commit 9d3127d, with thanks to
  Colin Diesh)


Release 1.6.1: 22nd August 2024
-------------------------------

This release is primarily portability and minor bug fixes.

Changes

- Improve warning levels by the compiler in CI. (#125)

- Switch to GitHub actions for some CI builds. (#121, #123)

- Add configure check for cpuid systems. (#115, #116.  Reported by
  Ryan Carsten Schmidt)

Bug fixes

- Use unsigned chars for ctype macros in the name tokeniser.
  On many systems this was already mitigated against, but on some OSes
  a char > 128 could trigger a buffer underrun. (#124)

- Fix interaction between _XOPEN_SOURCE and FreeBSD.
  (#119, John Marshall)

- Improve AVX512 compiler support, notably MacOS El Capitan's XCode.
  (#118, Rob Davies)

- Fix -std=c99 -pendantic pedantry (#117)


Release 1.6.0: 7th December 2023
--------------------------------

This release is primarily bug fixes, mostly spotted through improved fuzz
testing.

One big change however is the SIMD rANS codecs are now performant on Intel
CPUs with the DownFall mitigation microcode applied.


Changes

- Replaced the rANS codec SIMD gathers with simulated gathers via scalar
  memory fetches.  This helps AMD Zen4, but importantly it also fixes a
  disastrous performance regression caused by Intel's DownFall microcode fix.

  There is an impact on pre-DownFall speeds, but we should focus on patched
  CPUs as a priority.

- A small speed up to the rans_F_to_s3 function used by order-0 rans decode.

- Small speed up to SIMD rans32x16 order-1 encoder by reducing cache misses.
  Also sped up the rans4x8 order-1 encoder, particularly on AMD Zen4.

- Now supports building with "zig cc"
  (Issue #109, reported by David Jackson)


Bug fixes

- Improve robustness of name tokeniser when given non 7-bit ASCII and on
  machines where "char" defaults to unsigned.
  (Issue #105, reported by Shubham Chandak)

- Also fixed a 1 byte buffer read-overrun in name tokeniser.

- Fix name tokeniser encoder failure with some duplicated streams.

- Fixed rans_set_cpu to work multiple times, as well as reinstating the
  ability to change decode and encode side independently (accidentally lost in
  commit 958032c).  No effect on usage, but it improves the test coverage.

- Added a round-trip fuzz tester to test the ability to encode.  The old fuzz
  testing was decode streams only.

- Fixed bounds checking in rans_uncompress_O0_32x16_avx2, fixing buffer read
  overruns.

- Removed undefined behaviour in transpose_and_copy(), fixing zig cc builds.


Release 1.5.2: 6th October 2023
-------------------------------

*** SECURITY FIXES ***

This release contains multiple bug fixes, including a couple
buffer overruns that could corrupt memory when used in specific
scenarios.  These have not been observed with real data, but could
represent an attack vector for a malicious user.  (We know of no
exploit.)


Changes

- The range coder has been extended to do bounds checking if the
  new RC_SetOutputEnd() is called.  This has a small performance hit
  for the encoder, depending on compiler, but tests showed within 10%
  at worst.

Bug fixes

- Fix write-buffer overruns in fqzcomp and name tokeniser.

  SECURITY ISSUE: FQZComp could overflow the computed maximum growth
  size, causing writes beyond the ends of the allocated memory.  This
  is triggered by many very small 1bp reads.  Fixed the maximum
  bounds for compressed data.

  SECURITY ISSUE: The name tokeniser using the maximum number of
  tokens (128) would erroneously write a 129th token.  This is a
  restricted overflow of a few bytes.

  (PR#97, reported by Shubham Chandak)

- Fix an maximum 8-byte read overflow in the AVX2 rans decoder.
  SECURITY ISSUE: This was only present when using gcc.
  (PR#100, reported by Rob Davies)

- The rANS Order-1 SSE4 decoder could decode incorrectly.
  When a single symbol only occurs and we're using 12-bit freqs, the
  frequency of 4096 was interpreted as freq 0.  This only happens in
  the non-SIMD tidy-up stage at the end of the decode, so at worst the
  final 31 bytes may be incorrect. (PR#102)

- Fixed a 1-byte heap read-buffer overflow. Existed since 6a87ead2
  (Oct 2021).  Low severity security due to size and high likelihood
  it's just malloc meta-data. (PR#95; OSS-Fuzz 62270)

- rans_compress_4x16 now works on zero length input.
  Previously this was giving divide-by-zero errors.
  (PR#101, reported by Shubham Chandak)

- Remove asserts which caused warnings about unused variables when
  building with -DNDEBUG.

- Fix ARM builds when HWCAP_ASIMD is missing (on Conda) (PR#91)

- Improve FreeBSD CI testing

- Fix undefined behaviour from signed bit-shifting (PR#90).


Release 1.5.1: 19th July 2023
-----------------------------

This release is mainly small updates and bug fixes focusing on
specific platforms, with no new features added.

Changes

- Be more selective in use of AVX512 on AMD Zen4 processors.  This can
  be faster (e.g. with 64-way unrolling), but in the current rANS codec
  implementations AVX2 is faster for certain operations (PR#85).

- Add config.h to test programs to help them pick up definitions such
  as XOPEN_SOURCE (PR#84)

- Add FreeBSD to CI testing (PR#83)

Bug fixes

- Trivial bug fix to the rans4x16pr test harness when given
  incompressible data (PR#86).

- Make ARM NEON checks specific to AArch64 and exclude AArch32 systems.
  (PR#82 to fix issue#81, reported by Robert Clausecker)


Release 1.5.0: 14th April 2023
------------------------------

Changes

- Significant speed ups to the fqzcomp codec via code restructuring
  and use of memory prefetch instructions.  Encode is 30-40% faster
  and decode 5-8% faster. (PR#75 James Bonfield)

- Improve multiarch builds on MacOS, fixing issues with getting the
  various SIMD implementations integrated. (Issue#76 John Marshall,
  PR#77/#78 Rob Davies)

- Remove unused ax_with_libdeflate.m4 file from build system.


Release 1.4.0: Februrary 2023
-----------------------------

This is almost entirely minor bug fixing with a few small updates.

Changes

- Optimise compression / speed of the name tokeniser.
  - In arithmetic coding mode, it can now utilise bzip2 at higher levels.
  - For both rans / arith entropy encoders, the choice of method / order
    is now optimised per token type, giving faster compression.
  - Culled a pointless zlib check in the configure script.
  - Made lack of bzip2 a hard failure in configure, unless an explicit
    --disable-bz2 option is given.
  (#72, #73)

- Switch CI to use ARM for MacOS builds
  (#69, thanks to Rob Davies)


Bug fixes

- Remove some newer compiler warnings (#61)

- Improvements for Intel -m32 builds, including better AVX2 validation
  (m32 misses _mm256_extract_epi64) and improved data alignment.
  (#62. See also samtools/htslib#1500)

- Detect Neon capability at runtime via operating system APIs.
  (#63, thanks to John Marshall)

- Improve FreeBSD diagnostics when neglecting to use -lpthread / -lthr. 
  Plus additional extra error checking too.
  (#68, #64, thanks to John Marshall)

- Update hts_pack to operate in line with CRAMcodecs spec, where the
  number of symbols > 16.
  (#65/#66, reported by Michael Macias)

- Fixed too-stringent buffer overflow checking in O1 rans decoder.
  (#71, reported by Divon Lan)


Release 1.3.0: 9th August 2022
------------------------------

The primary change in this release is a new SIMD enabled rANS codec.

Changes

- There is a 32-way unrolled rANS implementation.  This is accessed
  using the existing rans 4x16 API with the RANS_ORDER_X32 bit set.
  Implementations exist for SSE4.1, AVX2, AVX512 and ARM Neon, as
  well as traditional non-SIMD scalar code in C and JavaScript. See
  the commit logs for benchmarks.

- Improved memory allocation via a new htscodecs_tls_alloc function.
  This uses Thread Local Storage (TLS) to avoid multiple malloc/free
  calls, reducing system CPU time.

- Some external functions have been renamed, with the old ones still
  existing in a deprecated fashion.  Every symbol should now start
  hts_, rans_, arith_, fqz_ or tok3_*.

- Improved test framework with an "entropy" tool that iterates over
  all entropy encoders.

- Updated the Appveyor CI image to user a newer gcc.  Also added ARM
  to the list of processors to test on.

- Tab vs space code changes.  Use "git diff -w" to see through these.

- Reworked fuzzing infrastructure.

- Small speed improvements to various rANS encoders and decoders.
  These were tested on a broad range of compilers, versions and
  systems.  The new code may be slightly slower with some combinations,
  but is faster overall and removes a few outliers with considerably
  degraded performance.

- Substantial memory reduction to the name tokeniser (tok3).

Bug fixes

- Fixed undefined behaviour in our use of _builtin_clz().

- Fixed a few redundant #includes.

- Work around strict aliasing bugs, uncovered with gcc -O2.

- Fixed an issue with encoding data blocks close to 2GB in size.
  (Additionally blocks above 2GB now error, rather than crashing or
  returning incorrect results.)

- Fix encode error with large blocks using RANS_ORDER_STRIPE.


Release 1.2.2: 1st April 2022
-----------------------------

This release contains some fixes found during fuzzing with Clang's
memory-sanitizer.  None of these are involving writing memory so there
is no possibility for code execution vulnerabilities.  However some do
could access uninitialised elements in locally allocated memory, which
could leak private data if the library was used in conjunction with
other tools which don't zero sensitive data before freeing.

Bug fixes:

- The name tokeniser now validates the stored length in the data
  stream matches the actual decoded length.  Discovered by Taotao Gu.

- Fixed an endless loop in arith_dynamic and rans4x16pr involving
  X_STRIPE with 0 stripes.

- Avoid a harmless (and wrong?) undefined behaviour sanitizer error
  when calling memcpy(ptr, NULL, 0) in the name tokeniser.

- Fixed possible uninitialised memory access in
  rans_uncompress_O1_4x16.  If the frequency table didn't add up to
  the correct amount, parts of the "fb" table were left unpopulated.
  It was then possible to use these array elements in some of the rANS
  calculations.

- Similarly rans_uncompress_O0 could access an uninitialised element
  4095 of the decoder tables if the frequencies summed to 4095 instead
  of the expected 4096.

- Improved error detection from fqzcomp's read_array function.

- Reject fqzcomp parameters with inconsistent "sel" parameters, which
  could lead to uninitialised access to the model.sel range coder.


Release 1.2.1: 15th February 2022
---------------------------------

The only change in this release is a minor adjustment to the histogram
code so it works on systems with small stacks.  This was detected on
Windows Mingw builds.


Release 1.2: 10th February 2022
-------------------------------

This release contains the following minor changes.
Please see the "git log" for the full details.

Improvements / changes:

- Speed up of rANS4x16 order-0.  We now use a branchless encoder
  renormalisation step.  For complex data it's between 13 and 50%
  speed up depending on compiler.

- Improve rANS4x16 compute_shift estimates.  The entropy calculation
  is now more accurate.  This leads to more frequent use of the 10-bit
  frequency mode, at an expense of up to 1% size growth.

- Speed improvements to the striped rANS mode, both encoding and
  decoding.  Encoder gains ~8% and decoder ~5%, but varies
  considerably by compiler and data.

- Added new var_put_u64_safe and var_put_u32_safe interfaces.
  These are automatically used by var_put_u64 and var_put_u32 when
  near the end of the buffer, but may also be called directly.

- Small speed ups to the hist8 and hist1_4 functions.

- Minor speed up to RLE decoding.

Bug fixes:

- Work around an icc-2021 compiler bug, but also speed up the varint
  encoding too (#29).

- Fix an off-by-one error in the initial size check in arith_dynamic.
  This meant the very smallest of blocks could fail to decode.
  Reported by Divon Lan.

- Fixed hist1_4 to also count the last byte when computing T0[].

- Fixed overly harsh bounds checking in the fqzcomp read_array
  function, which meant it failed to decode some configurations.


Release 1.1.1: 6th July 2021
----------------------------

This release contains the following minor changes.
Please see the "git log" for the full details.

Improvements / changes:

- Modernised autoconf usage to avoid warnings with newer versions.
  (John Marshall)

- Avoid using awk with large records, due to some systems
  (e.g. Solaris / OpenIndiana) with line length limits .
  (John Marshall)

- Applied Debian patch to make the library link against -lm.

Bug fixes:

- Fixed an issue with the name tokeniser when a slice (name_context)
  has exactly 1 more name than the previous call. (James Bonfield)

- Removed access to an uninitialised variable in the name tokeniser
  decode when given malformed data.  This occurs when we use delta
  encoding for the very first name. (James Bonfield, OSS-Fuzz)

- Minor fixes to distcheck and distclean targets


Release 1.0: 23rd Feb 2021
--------------------------

This marks the first non-beta release of htscodecs, following a
perioid of integration with Htslib and automated fuzzing by Google's
OSS-Fuzz program.

[Note this testing only applies to the C implementation.  The
JavaScript code should still be considered as examples of the codecs,
more for purposes of understanding and clarity than as a fully
optimised and tested release.]

Since the last release (0.5) the key changes are:

- Improved support for big endian platforms

- Speed improvements to CRAM 3.0 4x8 rANS order-1 encoding.
  It's between 10 and 50% faster at encoding, based on input data.

- Improved autoconf bzip2 checks and tidy up "make test" output.

- Added some more files into "make install", so that "make distcheck"
  now passes.

- Replaced Travis with Cirrus-CI testing.

- Removed various C undefined behaviour, such as left shifting of
  negative values and integer overflows.  As far as we know these were
  currently harmless on the supported platforms, but may break future
  compiler optimisations.

- Fixed numerous OSS-Fuzz identified flaws.  Some of these were
  potential security issues such as small buffer overruns.

- Tidied up some code to prevent warnings.

- The name tokeniser now has a limit on the size of data it can encode
  (10 million records).  This may still be too high given the memory
  it will require, so it may be reduced again.

