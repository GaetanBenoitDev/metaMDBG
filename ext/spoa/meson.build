project(
  'spoa',
  ['cpp'],
  version : '4.1.5',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++17',
    'b_ndebug=if-release'],
  license : 'MIT',
  meson_version : '>=0.60.0'
)

################
# Dependencies #
################

spoa_lib_flags = []
spoa_lib_deps = []

if get_option('cereal')
  # cereal
  spoa_cereal_dep = dependency('cereal', version : '>= 1.3.0', fallback : ['cereal', 'cereal_dep'])

  spoa_lib_flags += '-DSPOA_USE_CEREAL'
  spoa_lib_deps += spoa_cereal_dep
endif

if (not meson.is_subproject()) and (get_option('exe') or get_option('tests'))
  # biosoup
  spoa_biosoup_dep = dependency('biosoup', version : '>= 0.11.0', fallback : ['biosoup', 'biosoup_dep'])

  # bioparser
  spoa_bioparser_dep = dependency('bioparser', version : '>= 3.1.0', fallback : ['bioparser', 'bioparser_dep'])
endif

if not get_option('dispatch')
  if get_option('sse41')
    spoa_lib_flags += '-msse4.1'
  elif get_option('avx2')
    spoa_lib_flags += '-mavx2'
  endif
endif

if (get_option('simde') or
    get_option('simde_nonvec') or
    get_option('simde_openmp') or
    get_option('dispatch'))
  # simde
  spoa_simde_dep = dependency('simde', version : '>= 0.7.6', fallback : ['simde', 'simde_dep'])

  add_project_arguments('-DSIMDE_ENABLE_NATIVE_ALIASES', language : ['c', 'cpp'])
  spoa_lib_flags += '-DSPOA_USE_SIMDE'
  spoa_lib_deps += spoa_simde_dep

  if get_option('simde_nonvec')
    add_project_arguments('-DSIMDE_NO_NATIVE', language : ['c', 'cpp'])
  endif

  if get_option('simde_openmp')
    add_project_arguments('-DSIMDE_ENABLE_OPENMP', language : ['c', 'cpp'])
    spoa_lib_flags += '-fopenmp-simd'
  endif

  if get_option('dispatch')
    spoa_lib_flags += '-DSPOA_GENERATE_DISPATCH_CPUIDEX'
  endif
endif

###########
# Headers #
###########

subdir('include')

###########
# Sources #
###########

subdir('src')

#########
# Tests #
#########

if (not meson.is_subproject()) and get_option('tests')
  # gtest
  spoa_gtest_dep = dependency('gtest', version : '>= 1.10.0', main : true, fallback : ['gtest', 'gtest_main_dep'])

  subdir('test')
endif

###################
# Dependency info #
###################

if (not meson.is_subproject())
  import('pkgconfig').generate(
    spoa_lib,
    name : 'spoa',
    version : meson.project_version(),
    filebase : 'spoa',
    description : 'C++ tool/library for SIMD vectorized partial order alignment.')
endif

spoa_dep = declare_dependency(
  include_directories : spoa_include_directories,
  link_with : spoa_lib,
  dependencies : spoa_lib_deps,
  version : meson.project_version()
)
