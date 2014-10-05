#ifndef TEST_TEST_UTIL_H_
#define TEST_TEST_UTIL_H_

#include <random>
#include <vector>
#include <Eigen/Sparse>

namespace sskkm {

inline bool operator==(
    const SparseMatrix& matrix1, const SparseMatrix& matrix2) {
  if (matrix1.rows() != matrix2.rows() || matrix1.cols() != matrix2.cols()) {
    return false;
  }
  for (int i = 0; i < matrix1.outerSize(); ++i) {
    for (SparseMatrix::InnerIterator iter1(matrix1, i), iter2(matrix2, i);
         iter1 && iter2;
         ++iter1, ++iter2) {
      if (iter1.row() != iter2.row()) { return false; }
      if (iter1.col() != iter2.col()) { return false; }
      if (iter1.value() != iter2.value()) { return false; }
    }
  }
  return true;
}

inline ClusterIndicatorMatrix InitializeRandomClusters(
    unsigned num_vectors, unsigned num_clusters) {
  static std::mt19937 rng;
  std::uniform_int_distribution<> dist(0, num_clusters - 1);
  std::vector<SparseMatrixCoefficient> cluster_indicators;
  cluster_indicators.reserve(num_vectors);
  while (true) {
    std::vector<bool> cluster_has_vertices(num_vectors);
    for (unsigned i = 0; i < num_vectors; ++i) {
      unsigned idx_cluster = dist(rng);
      cluster_indicators.push_back(
          SparseMatrixCoefficient(i, idx_cluster, 1.0));
      cluster_has_vertices[idx_cluster] = true;
    }
    ClusterIndicatorMatrix clusters(num_vectors, num_clusters);
    for (unsigned i = 0; i < num_clusters; ++i) {
      if (!cluster_has_vertices[i]) { goto CONTINUE_OUTER_LOOP; }
    }
    clusters.setFromTriplets(
        cluster_indicators.begin(), cluster_indicators.end());
    return clusters;
 CONTINUE_OUTER_LOOP:
    cluster_indicators.clear();
  }
}

}  // namespace sskkm

#endif
