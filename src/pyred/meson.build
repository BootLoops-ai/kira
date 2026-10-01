
config_h_in = configure_file(output : 'config.h.in', configuration : conf_data)

config_h = vcs_tag(
  command : ['git', 'describe', '--tags', '--broken', '--always'],
  input : config_h_in,
  output : 'config.h',
  fallback : 'no')

pyred_src = [
  'coeff_helper.cpp',
  'coeff_int.cpp',
  'coeff_vec.cpp',
  'integrals.cpp',
  'interface.cpp',
  'keyvaluedb.cpp',
  'parser.cpp',
  'relations.cpp']

sqlite_src = ['../sqlite3/sqlite3.c']

gzstream_src = ['../gzstream/gzstream.cpp']

pyred = static_library(
  'pyred',
  pyred_src + sqlite_src + gzstream_src + [config_h],
  include_directories : include_dirs,
  dependencies : deps,
  install : false)
