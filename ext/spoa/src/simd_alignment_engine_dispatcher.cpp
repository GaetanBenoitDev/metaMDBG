// Copyright (c) 2020 Mario Brcic, Robert Vaser

#include "simd_alignment_engine_implementation.hpp"

#if defined(SPOA_GENERATE_DISPATCH)

#include "cpuinfo_x86.h"  // NOLINT

static const cpu_features::X86Features features =
    cpu_features::GetX86Info().features;

#elif defined(SPOA_GENERATE_DISPATCH_CPUIDEX)

#  if defined(_MSC_VER)

#include <intrin.h>

#  endif

// adapted from https://github.com/intel/linux-sgx/blob/master/common/inc/internal/linux/cpuid_gnu.h
void RunCpuidex(int cpu_info[4], int function_id, int subfunction_id) {
#  if defined(_MSC_VER)
  __cpuidex(cpu_info, function_id, subfunction_id);
#  elif defined(__X86_64__)
  __asm__ volatile ("cpuid")
                  : "=a" (cpu_info[0]), "=b" (cpu_info[1]), "=c" (cpu_info[2]), "=d" (cpu_info[3])
                  : "0" (function_id), "2" (subfunction_id));
#  else
  __asm__ volatile ("xchgl %%ebx, %1; cpuid; xchgl %%ebx, %1"
                  : "=a" (cpu_info[0]), "=r" (cpu_info[1]), "=c" (cpu_info[2]), "=d" (cpu_info[3])
                  : "0" (function_id), "2" (subfunction_id));
#  endif
}

constexpr int SPOA_SSE2 = 0x1;
constexpr int SPOA_SSE4_1 = 0x2;
constexpr int SPOA_AVX2 = 0x4;

static int GetX86Info() {
  int cpu_info[4] = {0};

  RunCpuidex(cpu_info, 0, 0);

  const int n = cpu_info[0];
  if (n == 0) {
    return 0;
  }

  RunCpuidex(cpu_info, 1, 0);

  int features = 0;
  if (cpu_info[3] >> 26) {
    features |= SPOA_SSE2;
  }
  if (cpu_info[2] >> 19) {
    features |= SPOA_SSE4_1;
  }

  if (n >= 7) {
    RunCpuidex(cpu_info, 7, 0);

    if (cpu_info[1] >> 5) {
      features |= SPOA_AVX2;
    }
  }

  return features;
}

static const int features = GetX86Info();

#endif

namespace spoa {
 
#if !defined(SPOA_GENERATE_DISPATCH) && !defined(SPOA_GENERATE_DISPATCH_CPUIDEX)

template class SimdAlignmentEngine<Architecture::kAutomatic>;

#endif

std::unique_ptr<AlignmentEngine> CreateSimdAlignmentEngine(
    AlignmentType type,
    AlignmentSubtype subtype,
    std::int8_t m,
    std::int8_t n,
    std::int8_t g,
    std::int8_t e,
    std::int8_t q,
    std::int8_t c) {
#if defined(SPOA_GENERATE_DISPATCH)
  if (features.avx2) {
    return SimdAlignmentEngine<Architecture::kAVX2>::Create(
        type, subtype, m, n, g, e, q, c);
  } else if (features.sse4_1) {
    return SimdAlignmentEngine<Architecture::kSSE4_1>::Create(
        type, subtype, m, n, g, e, q, c);
  } else {
    return SimdAlignmentEngine<Architecture::kSSE2>::Create(
        type, subtype, m, n, g, e, q, c);
  }
#elif defined(SPOA_GENERATE_DISPATCH_CPUIDEX)
  if (features & SPOA_AVX2) {
    return SimdAlignmentEngine<Architecture::kAVX2>::Create(
        type, subtype, m, n, g, e, q, c);
  } else if (features & SPOA_SSE4_1) {
    return SimdAlignmentEngine<Architecture::kSSE4_1>::Create(
        type, subtype, m, n, g, e, q, c);
  } else {
    return SimdAlignmentEngine<Architecture::kSSE2>::Create(
        type, subtype, m, n, g, e, q, c);
  }
#else
  return SimdAlignmentEngine<Architecture::kAutomatic>::Create(
      type, subtype, m, n, g, e, q, c);
#endif
}

}  // namespace spoa
