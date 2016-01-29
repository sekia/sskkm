sskkm
===

sskkm is a C++ library for graph/vector clustering.

This library implements a semi-supervised graph clustering algorithm employing a vector clustering method, which is described in "Semi-supervised graph clustering: a kernel approach" (Brian Kulis, et al., 2009).

Requirements
---

sskkm depends on 2 external C++ libraries:

- boost 1.50 -- http://www.boost.org/
- Eigen 3.1 (bundled with this distribution) -- http://eigen.tuxfamily.org/

Installation
---

Download the distribution and run `tar xzf sskkm-x.x.x.tar.gz`. Then you should have an extracted directory. Its structure is:

```
sskkm/
|---README          - This file.
|---include/        - C++ headers.
|---sample/         - Sample programs.
|---test/           - Unit tests.
|---{sample,test}_* - Working directories for building sample programs/tests.
```

Since sskkm is a header-only library, no installation instruction is needed. You can just copy the header files under `/include` directory to wherever compiler knows then get ready.

If you want to build a unit test and sample programs, Install [OMake](http://omake.metaprl.org/) and run `omake` on the distrinbution's root directory. And you should have executables at `{sample,test}_*` directories.

Documentation
---

Currently we have no documentation. Sorry!

You can learn usage of the library via sample applications in `/sample` directory.

Copyrigth and License
---

Copyright (c) 2012-2016 Koichi SATOH, All rights reserved.

The MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
