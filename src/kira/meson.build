
kira_src = [
  'main.cpp',
  'connect2kira.cpp',
  'dataBase.cpp',
  'dgl.cpp',
  'ibp.cpp',
  'integral.cpp',
  'kira.cpp',
  'ReadYamlFiles.cpp',
  'reduction.cpp',
  'threadpool.cpp',
  'tools.cpp',
  'trivial_sym.cpp']

if get_option('firefly')
  kira_src += [
  'reductionFF.cpp',
  'mpi_work.cpp']
endif

kira_exe_name = 'kira'
if get_option('weight_width') != '64'
  kira_exe_name = kira_exe_name + get_option('weight_width')
endif

kira_exe = executable(
  kira_exe_name,
  kira_src + [config_h],
  include_directories : include_dirs,
  link_with : pyred,
  dependencies : deps + deps_kira,
  install_rpath : rpath,
  install : true)
