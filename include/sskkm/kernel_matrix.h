#ifndef SSKKM_KERNEL_MATRIX_H_
#define SSKKM_KERNEL_MATRIX_H_

#include <cmath>
#include <functional>

#include "sskkm/base.h"
#include "sskkm/kernel_k_means.h"

namespace sskkm {

inline KernelMatrix ComputePolynomialKernelMatrix(
    const DenseMatrix& vectors, double linear_shift, double exponent) {
  return KernelMatrix(
      ((vectors.transpose() * vectors).array() + linear_shift).pow(exponent));
}

inline KernelMatrix ComputePolynomialKernelMatrix(
    const SparseMatrix& vectors, double linear_shift, double exponent) {
  DenseMatrix products = (vectors.transpose() * vectors).eval();
  return KernelMatrix((products.array() + linear_shift).pow(exponent));
}

template <typename Matrix>
inline KernelMatrix ComputeGaussianKernelMatrix(
    const Matrix& vectors, double deviation) {
  DenseMatrix products = (vectors.transpose() * vectors).eval();
  KernelMatrix exponents =
      (products.diagonal().replicate(1, vectors.cols()) +
        (-2 * products) +
        products.diagonal().transpose().replicate(vectors.cols(), 1)) /
      -(deviation * deviation);
  return exponents.unaryExpr(std::ptr_fun<double, double>(std::exp));
}

}  // namespace sskkm

#endif  // SSKKM_KERNEL_MATRIX_H_
