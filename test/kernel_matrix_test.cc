#include "gtest/gtest.h"
#include "sskkm/kernel_matrix.h"

using namespace sskkm;

TEST(KernelMatrixTest, ComputePolynomialKernelMatrix) {
  DenseMatrix dense_vectors(4, 3);
  dense_vectors <<
      1, 1, 1,
      1, 1, 0,
      0, 1, 1,
      0, 1, 1;

  SparseMatrix sparse_vectors(4, 3);
  sparse_vectors.insert(0, 0) = 1;
  sparse_vectors.insert(1, 0) = 1;
  sparse_vectors.insert(0, 1) = 1;
  sparse_vectors.insert(1, 1) = 1;
  sparse_vectors.insert(2, 1) = 1;
  sparse_vectors.insert(3, 1) = 1;
  sparse_vectors.insert(0, 2) = 1;
  sparse_vectors.insert(2, 2) = 1;
  sparse_vectors.insert(3, 2) = 1;
  sparse_vectors.finalize();

  KernelMatrix kernels_expected(3, 3);
  kernels_expected <<
      9, 9, 4,
      9, 25, 16,
      4, 16, 16;
  EXPECT_EQ(kernels_expected,
            ComputePolynomialKernelMatrix(dense_vectors, 1, 2));
  EXPECT_EQ(kernels_expected,
            ComputePolynomialKernelMatrix(sparse_vectors, 1, 2));
}
