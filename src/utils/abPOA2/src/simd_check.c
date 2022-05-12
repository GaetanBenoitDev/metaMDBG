#include "simd_instruction.h"


//https://stackoverflow.com/questions/152016/detecting-cpu-architecture-compile-time                                                         
#if MSVC
#ifdef _M_X86
#define ARCH_X86
#endif
#endif

#if GCC
#ifdef __i386__
#define ARCH_X86
#endif
#endif

#ifndef ARCH_X86

int simd_check(void) {
  return SIMD_AVX2;
}
#else

#ifndef _MSC_VER
// adapted from https://github.com/01org/linux-sgx/blob/master/common/inc/internal/linux/cpuid_gnu.h
void __cpuidex(int cpuid[4], int func_id, int subfunc_id)
{
#if defined(__x86_64__)
	__asm__ volatile ("cpuid"
			: "=a" (cpuid[0]), "=b" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
			: "0" (func_id), "2" (subfunc_id));
#else // on 32bit, ebx can NOT be used as PIC code
	__asm__ volatile ("xchgl %%ebx, %1; cpuid; xchgl %%ebx, %1"
			: "=a" (cpuid[0]), "=r" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
			: "0" (func_id), "2" (subfunc_id));
#endif
}
#endif

int simd_check(void) {
	int flag = 0, cpuid[4], max_id;
	__cpuidex(cpuid, 0, 0);
    // int i;
    // for (i = 0 ; i < 4; ++i) printf("%d\t", cpuid[i]); printf("\n");
	max_id = cpuid[0];
	if (max_id == 0) return 0;
	__cpuidex(cpuid, 1, 0);
    // for (i = 0 ; i < 4; ++i) printf("%d\t", cpuid[i]); printf("\n");
	if (cpuid[3]>>25&1) flag |= SIMD_SSE;
	if (cpuid[3]>>26&1) flag |= SIMD_SSE2;
	if (cpuid[2]>>0 &1) flag |= SIMD_SSE3;
	if (cpuid[2]>>9 &1) flag |= SIMD_SSSE3;
	if (cpuid[2]>>19&1) flag |= SIMD_SSE41;
	if (cpuid[2]>>20&1) flag |= SIMD_SSE42;
	if (cpuid[2]>>28&1) flag |= SIMD_AVX;
	if (max_id >= 7) {
		__cpuidex(cpuid, 7, 0);
        // for (i = 0 ; i < 4; ++i) printf("%d\t", cpuid[i]); printf("\n");
		if (cpuid[1]>>5 &1) flag |= SIMD_AVX2;
		if (cpuid[1]>>16&1) flag |= SIMD_AVX512F;
		if (cpuid[1]>>30&1) flag |= SIMD_AVX512BW;
	}

	return flag;
}

#ifdef __CHECK_SIMD_MAIN__
int main(void) {
    char simd_label[6][20] = {"No SIMD", "SSE2 (128 bits)", "SSE4.1 (128 bits)", "AVX2 (256 bits)", "AVX512F (512 bits)", "AVX512BW (512 bits)"};
    int simd_flag = simd_check(), t=0;
    if (simd_flag & SIMD_AVX512BW) printf("__AVX512BW__\n"), t = 5;
    else if (simd_flag & SIMD_AVX512F) printf("__AVX512F__\n"), t = 4;
    else if (simd_flag & SIMD_AVX2) printf("__AVX2__\n"), t = 3;
    else if (simd_flag & SIMD_SSE41) printf("__SSE4_1__\n"), t = 2;
    else if (simd_flag & SIMD_SSE2) printf("__SSE2__\n"), t = 1;
    else printf("NO SIMD\n"), t = 0;

    char msg[100], i;
    fprintf(stderr, "\n");
    sprintf(msg, "==== %s will be used. ====", simd_label[t]);
    for (i = 0; msg[i]; ++i) fprintf(stderr, "="); fprintf(stderr, "\n");
    fprintf(stderr, "%s\n",msg);
    for (i = 0; msg[i]; ++i) fprintf(stderr, "="); fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    return simd_flag;
}
#endif

#endif
