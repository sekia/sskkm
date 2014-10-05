# -*- mode: python -*-
{
  'target_defaults': {
    'include_dirs': [
      '<(DEPTH)/include',
      '<(DEPTH)/eigen-3.1.4',
    ],
    'includes': [
      'configure.gypi',
    ],
    'configurations': {
      'Debug': {
        'xcode_settings': {
          'GCC_OPTIMIZATION_LEVEL': '0',
        },
      },
      'Release': {
        'defines': ['NDEBUG',],
        'xcode_settings': {
          'GCC_OPTIMIZATION_LEVEL': '3',
        },
      },
    },
  },
}
