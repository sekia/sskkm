#ifndef SSKKM_BASE_H_
#define SSKKM_BASE_H_

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace sskkm {

typedef Eigen::VectorXd DenseVector;

typedef Eigen::MatrixXd DenseMatrix;

typedef Eigen::SparseVector<double> SparseVector;

typedef Eigen::SparseMatrix<double> SparseMatrix;

typedef Eigen::Triplet<double> SparseMatrixCoefficient;

/**
   N-size vector where N is the number of nodes. Each element indicates whether
   node having corresponding order belongs to a cluster (value 1.0) or not
   (0.0).
 */
typedef SparseVector ClusterIndicatorVector;

/**
   N x C matrix where N is the number of nodes and C is the number of clusters.
   The element at (n, c) is 1.0 if the n-th node belongs to the c-th cluster,
   Otherwise 0.0. i.e., the c-th column vector of the matrix is a
   ClusterIndicatorVector for the c-th cluster.
 */
typedef SparseMatrix ClusterIndicatorMatrix;

class InvalidArgument : public std::invalid_argument {
public:
  explicit InvalidArgument(const std::string &what_arg)
      : std::invalid_argument(what_arg) {}
};

class RuntimeError : public std::runtime_error {
public:
  explicit RuntimeError(const std::string &what_arg)
      : std::runtime_error(what_arg) {}
};

/**
   Reads entire bytes of a stream.
 */
inline std::string SlurpStream(std::istream &in) throw (InvalidArgument) {
  if (!in.good()) { throw InvalidArgument("The stream is not readable"); }
  std::ostringstream buf;
  buf << in.rdbuf();
  return buf.str();
}

/**
   Reads whole content in a file.
 */
inline std::string SlurpFile(const char *filename) throw (InvalidArgument) {
  std::ifstream file(filename);
  return SlurpStream(file);
}

}  // namespace sskkm

#endif  // SSKKM_BASE_H_
