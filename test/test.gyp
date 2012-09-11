# -*- mode: python -*-
{
  'targets': [
    {
      'target_name': 'unit_test',
      'type': 'executable',
      'sources': [
        'cluto_matrix_test.cc',
        'kernel_k_means_test.cc',
        'kernel_matrix_test.cc',
      ],
      'include_dirs': [
        '<(DEPTH)', 
      ],
      'link_settings': {
        'libraries': ['-lgtest', '-lgtest_main',],
      },
    },
  ],
}
