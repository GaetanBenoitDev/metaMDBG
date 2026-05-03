###########
# Headers #
###########

spoa_include_directories = include_directories(['.'])

if meson.is_subproject()
  subdir_done()
endif

install_headers(
  files([
    'spoa/alignment_engine.hpp',
    'spoa/graph.hpp',
    'spoa/spoa.hpp',
    'spoa/version.hpp',
  ]),
  subdir : 'spoa')
