spoa_test_config = configuration_data()
spoa_test_config.set('SPOA_TEST_DATA',  meson.project_source_root() + '/test/data/sample.fastq.gz')

###########
# Sources #
###########

spoa_test_config_headers = [
  configure_file(
    input : 'spoa_test_config.h.in',
    output : 'spoa_test_config.h',
    configuration : spoa_test_config)
]

spoa_test_sources = files([
  'spoa_test.cpp',
]) + spoa_test_config_headers

spoa_test = executable(
  'spoa_test',
  spoa_test_sources,
  dependencies : [spoa_lib_deps, spoa_biosoup_dep, spoa_bioparser_dep, spoa_gtest_dep],
  include_directories : spoa_include_directories,
  link_with : spoa_lib,
  install : false)

#########
# Tests #
#########

test('spoa gtest unit tests', spoa_test)
