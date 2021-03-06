#ifndef SSKKM_BASE_H_
#define SSKKM_BASE_H_

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace sskkm {

using DenseVector = Eigen::VectorXd;

using DenseMatrix = Eigen::MatrixXd;

using SparseVector = Eigen::SparseVector<double>;

using SparseMatrix = Eigen::SparseMatrix<double>;

using SparseMatrixCoefficient = Eigen::Triplet<double>;

/**
   N-size vector where N is the number of nodes. Each element indicates whether
   node having corresponding order belongs to a cluster (value 1.0) or not
   (0.0).
 */
using ClusterIndicatorVector = SparseVector;

/**
   N x C matrix where N is the number of nodes and C is the number of clusters.
   The element at (n, c) is 1.0 if the n-th node belongs to the c-th cluster,
   Otherwise 0.0. i.e., the c-th column vector of the matrix is a
   ClusterIndicatorVector for the c-th cluster.
 */
using ClusterIndicatorMatrix = SparseMatrix;

/**
   Reads entire bytes of a stream.
 */
inline std::string SlurpStream(std::istream& in) {
  if (!in.good()) { throw std::invalid_argument("The stream is not readable"); }
  std::ostringstream buf;
  buf << in.rdbuf();
  return buf.str();
}

/**
   Reads whole content in a file.
 */
inline std::string SlurpFile(const char *filename) {
  std::ifstream file(filename);
  return SlurpStream(file);
}

}  // namespace sskkm

#endif  // SSKKM_BASE_H_
