spoa_config = configuration_data()
spoa_config.set('SPOA_VERSION', meson.project_version())

###########
# Library #
###########

spoa_config_headers = [
  configure_file(
    input : 'spoa_config.h.in',
    output : 'spoa_config.h',
    configuration : spoa_config)
]

spoa_sources = files([
  'alignment_engine.cpp',
  'graph.cpp',
  'simd_alignment_engine_dispatcher.cpp',
  'sisd_alignment_engine.cpp',
  'version.cpp',
]) + spoa_config_headers

spoa_simd_lib = []

if get_option('dispatch')
  simd = import('unstable_simd')

  spoa_simd_sources = files([
    'simd_alignment_engine_dispatch.cpp',
  ])

  spoa_simd_lib += simd.check(
    'spoa',
    sse2 : spoa_simd_sources,
    sse41 : spoa_simd_sources,
    avx2 : spoa_simd_sources,
    dependencies : spoa_lib_deps,
    include_directories : spoa_include_directories,
    cpp_args : spoa_lib_flags,
    compiler : meson.get_compiler('cpp'))[0]
endif

spoa_lib_install = (not meson.is_subproject()) or (get_option('default_library') != 'static')

spoa_lib = library(
  'spoa',
  spoa_sources,
  soversion : meson.project_version(),
  version : meson.project_version(),
  install : spoa_lib_install,
  dependencies : spoa_lib_deps,
  include_directories : spoa_include_directories,
  link_with : spoa_simd_lib,
  cpp_args : spoa_lib_flags)

##############
# Executable #
##############

if (not meson.is_subproject()) and get_option('exe')
  spoa_exe = executable(
    'spoa',
    files(['main.cpp']),
    install : true,
    dependencies : [spoa_lib_deps, spoa_biosoup_dep, spoa_bioparser_dep],
    include_directories : spoa_include_directories,
    link_with : spoa_lib)
endif
