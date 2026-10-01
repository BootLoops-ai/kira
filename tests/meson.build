
# meson test -v --suite kira

symcompare = executable(
  'symcompare',
  ['symcompare.cpp'] + [config_h],
  include_directories : include_dirs,
  link_with : pyred,
  dependencies : deps,
  install_rpath : rpath,
  install : false)

# mathematica = find_program('math', required : false)
# if mathematica.found()
#   message('FOUND MATHEMATICA')
# else
#   message('MATHEMATICA NOT FOUND')
# endif
# # mathematica.full_path() # Meson >= 0.55

end2end_test_script = find_program('run.sh')
input_dir = meson.current_source_dir() + '/input_data'
reference_dir = meson.current_source_dir() + '/reference_data'

if get_option('test-level') == 'sym'
  test_level = 0
elif get_option('test-level') == 'masters'
  test_level = 1
elif get_option('test-level') == 'full'
  test_level = 2
else
  error('Invalid test level: ' + get_option('test-level'))
endif

test_projects = [
  # [project_name,
  #  job_id,
  #  [test_level_0_suffix, test_level_1_suffix, test_level_2_suffix]]
  # --> run job_<job_id>_<test_level_suffix>.yaml
  # Results are compared to project_name/results_<job_id>.
  ['hqet2L_III_H', '1', ['1','2','3']],
  ['prop2L', '1', ['1','2','3']]
]

foreach test_project : test_projects
  test(
    test_project[0] + '__' + test_project[1],
    end2end_test_script,
    args : [kira_exe.full_path(),
            test_project[1],
            symcompare.full_path(),
            input_dir,
            reference_dir,
            test_project[0],
            test_project[2][test_level]],
    env : 'MALLOC_PERTURB_=0',
    suite : ['kira'],
    is_parallel : true)
endforeach
