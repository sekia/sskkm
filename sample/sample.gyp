# -*- mode: python -*-
{
  'targets': [
    {
      'target_name': 'kernel_k_means',
      'type': 'executable',
      'sources': ['kernel_k_means.cc',],
    },
    {
      'target_name': 'graph_clustering',
      'type': 'executable',
      'sources': ['graph_clustering.cc',],
      'link_settings': {
        'libraries': ['-lboost_graph',],
      },
    },
  ],
}
