import os, platform, sys

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

cmdclass = {}

try:
    from Cython.Build import build_ext
except ImportError: # without Cython
    #module_src = 'python/pyabpoa.c'
    sys.stderr.write('Error: \'cython\' is required to install pyabpoa\n')
    sys.exit(0)
else: # with Cython
    module_src = 'python/pyabpoa.pyx'
    cmdclass['build_ext'] = build_ext


simde = '-DUSE_SIMDE -DSIMDE_ENABLE_NATIVE_ALIASES'
sys.path.append('python')

if platform.machine() in ["aarch64", "arm64"]:
    simd_flag = '-march=armv8-a+simd -D__AVX2__'
elif platform.machine() in ["aarch32"]:
    simd_flag = '-march=armv8-a+simd -mfput=auto -D__AVX2__'
else:
    simd_flag='-march=native'
    if os.getenv('SSE4', False):
        simd_flag='-msse4.1'
    elif os.getenv('SSE2', False):
        simd_flag='-msse2'
    elif os.getenv('AVX2', False):
        simd_flag='-mavx2'
    #elif os.getenv('AVX512F', False):
    #    simd_flag='-mavx512f'
    #elif os.getenv('AVX512BW', False):
    #    simd_flag='-mavx512bw'

src_dir='src/'
inc_dir='include/'

src=[module_src, src_dir+'abpoa_align.c', src_dir+'abpoa_graph.c', src_dir+'abpoa_output.c', src_dir+'abpoa_plot.c', src_dir+'abpoa_seed.c', src_dir+'abpoa_seq.c', src_dir+'kalloc.c', src_dir+'kstring.c', src_dir+'simd_abpoa_align.c', src_dir+'simd_check.c', src_dir+'utils.c']

long_description = open('python/README.md').read()

setup(
    # Information
    name = "pyabpoa",
    description = "pyabpoa: SIMD-based partial order alignment using adaptive band",
    long_description = long_description,
    long_description_content_type="text/markdown",
    version = "1.4.0",
    url = "https://github.com/yangao07/abPOA",
    author = "Yan Gao",
    author_email = "gaoy286@mail.sysu.edu.cn",
    license = "MIT",
    keywords = "multiple-sequence-alignment  partial-order-graph-alignment",
    # Build instructions
    ext_modules = [Extension("pyabpoa",
                    sources=src,
                    include_dirs=[inc_dir],
                    depends=[src_dir+'abpoa.h', src_dir+'abpoa_align.h', src_dir+'abpoa_graph.h', src_dir+'abpoa_output.h', src_dir+'abpoa_seed.h', src_dir+'abpoa_seq.h', src_dir+'kalloc.h', src_dir+'khash.h', src_dir+'kdq.h', src_dir+'kseq.h', src_dir+'ksort.h', src_dir+'kstring.h', src_dir+'kvec.h', src_dir+'simd_abpoa_align.h', src_dir+'simd_instruction.h', src_dir+'utils.h', 'python/cabpoa.pxd'],
                    libraries = ['z', 'm', 'pthread'],
                    extra_compile_args=['-O3', '-Wno-error=declaration-after-statement', simde, simd_flag])],
    install_requires=['cython'],
    cmdclass = cmdclass
)
