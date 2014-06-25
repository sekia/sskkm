#ifndef SSKKM_KERNEL_K_MEANS_H_
#define SSKKM_KERNEL_K_MEANS_H_

#include <string>
#include <vector>

#include "sskkm/base.h"

namespace sskkm {

/**
   N x N matrix where N is the number of vectors. The element at (n, m) is 
   K(nv, mv) where K is the kernel function, and nv, mv are the n-th and m-th 
   vectors respectively.
 */
typedef DenseMatrix KernelMatrix;

/**
   N x C matrix where N is the number of vectors and C is the number of 
   clusters. The element at (n, c) is the distance between the n-th vector and 
   centroid of the c-th cluster in a kernel space.
 */
typedef DenseMatrix NormMatrix;

/**
   N-length vector where N is the number of vectors. Each scalar represents the 
   weight of corresponding vectors.
 */
typedef DenseVector WeightVector;

namespace internal {

/*
  A distance between the i-th vector (a_i) and centroid of the c-th cluster 
  (\pi_c) is computed by an expression:

  K_{ii} 
  - \frac {2 \sum_{a_j     \in \pi_c} K_{ij}} {|\pi_c|} 
  + \frac {  \sum_{a_j,a_k \in \pi_c} K_{jk}} {|\pi_c|^2}

  Reference: Brian Kulis, et al. 2009. Semi-supervised graph clustering: a 
             kernel approach. p.4.

  For optimization, we can omit the 1st term because it is fixed when comparing 
  distances between a vector and each cluster centroid.
*/
inline NormMatrix ComputeNormMatrix(
    const ClusterIndicatorMatrix &clusters,
    const KernelMatrix &kernels) {
  // Each element is the number of vectors in corresponding cluster.
  Eigen::ArrayXd cluster_sizes(clusters.cols());
  for (ClusterIndicatorMatrix::Index i = 0; i < clusters.cols(); ++i) {
    cluster_sizes(i) = clusters.col(i).nonZeros();
  }

  /*
    K * C, where K is a kernel matrix and C is a cluster indicator matrix, is a 
    N x C matrix. An element at (n, c) is sum of the kernels between the n-th 
    vector and each vector which is in the c-th cluster.
  */
  DenseMatrix kc = kernels * clusters;
  DenseMatrix term2s = -2 * kc;
  for (DenseMatrix::Index i = 0; i < term2s.cols(); ++i) {
    term2s.col(i) /= cluster_sizes(i);
  }

  /*
    diag((K*C)^t * C) = diag(C^t * K * C), where K is a kernel matrix and C is 
    a cluster indicator matrix, is a vector which contains sums of the kernels 
    between every pair of vectors which are in the corresponding cluster.
  */
  DenseVector term3s =
      Eigen::ArrayXd((kc.transpose() * clusters).diagonal()) /
      (cluster_sizes * cluster_sizes);

  return term2s + term3s.transpose().replicate(clusters.rows(), 1);
}

/*
   Weighted version of |ComputeNormMatrix|. The distance calculation is 
   modified to:

   K_{ii} 
   - \frac {2 \sum_{a_j     \in \pi_c} w_j     K_{ij}} 
           {  \sum_{a_j     \in \pi_c} w_j} 
   + \frac {  \sum_{a_j,a_k \in \pi_c} w_j w_k K_{jk}} 
           { (\sum_{a_j     \in \pi_c} w_j)^2}

   Where w_j, w_k are the weights for a_j, a_k respectively. We can omit the 
   1st term as well as the unweighted version.
*/
inline NormMatrix ComputeWeightedNormMatrix(
    const ClusterIndicatorMatrix &clusters,
    const KernelMatrix &kernels,
    const WeightVector &weights) {
  DenseVector cluster_total_weights = clusters.transpose() * weights;

  DenseMatrix kc = kernels.cwiseProduct(
      weights.transpose().replicate(kernels.rows(), 1)) * clusters;
  DenseMatrix term2s = -2 * kc.cwiseQuotient(
      cluster_total_weights.transpose().replicate(kc.rows(), 1));

  /*
    Although this expression is equivalent to another one used in an overloaded 
    function below, makes a slightly different result because of rounding 
    errors due to difference in the order of evaluation.
  */
  DenseVector term3s = Eigen::ArrayXd(
      (kc.cwiseProduct(weights.replicate(1, kc.cols())).transpose() *
       clusters).diagonal()) /
      (cluster_total_weights.array() * cluster_total_weights.array());

  return term2s + term3s.transpose().replicate(clusters.rows(), 1);
}

inline NormMatrix ComputeWeightedNormMatrix(
    const ClusterIndicatorMatrix &clusters,
    const KernelMatrix &kernels,
    const KernelMatrix &weighted_kernels,
    const KernelMatrix &weighted_weighted_kernels,
    const WeightVector &weights) {
  DenseVector cluster_total_weights = clusters.transpose() * weights;

  DenseMatrix term2s = -2 * (weighted_kernels * clusters).cwiseQuotient(
      cluster_total_weights.transpose().replicate(kernels.rows(), 1));

  DenseVector term3s = Eigen::ArrayXd(
      ((weighted_weighted_kernels * clusters).transpose() *
       clusters).diagonal()) /
      (cluster_total_weights.array() * cluster_total_weights.array());

  return term2s + term3s.transpose().replicate(clusters.rows(), 1);
}

/**
   Re-assigns vectors based on updated |norms|. Note that number of clusters may
   decrease.
 */
inline ClusterIndicatorMatrix ComputeClusterIndicatorMatrix(
    const NormMatrix &norms) {
  std::vector<SparseMatrixCoefficient> triplets;
  for (NormMatrix::Index i = 0; i < norms.rows(); ++i) {
    NormMatrix::Index cluster_index;
    norms.row(i).minCoeff(&cluster_index);
    triplets.push_back(SparseMatrixCoefficient(i, cluster_index, 1.0));
  }

  ClusterIndicatorMatrix temp_clusters(norms.rows(), norms.cols());
  temp_clusters.setFromTriplets(triplets.begin(), triplets.end());

  int num_empty_clusters = 0;
  for (ClusterIndicatorMatrix::Index i = 0; i < temp_clusters.cols(); ++i) {
    if (temp_clusters.col(i).nonZeros() == 0) { ++num_empty_clusters; }
  }
  ClusterIndicatorMatrix clusters(
      temp_clusters.rows(), temp_clusters.cols() - num_empty_clusters);
  for (ClusterIndicatorMatrix::Index i = 0, j = 0;
       i < clusters.cols();
       ++i, ++j) {
    while (temp_clusters.col(j).nonZeros() == 0) { ++j; }
    clusters.col(i) = temp_clusters.col(j);
  }
  return clusters;
}

inline void ComputeWeightedKernelMatrix(
    const KernelMatrix &kernels,
    const WeightVector &weights,
    KernelMatrix *weighted_kernels,
    KernelMatrix *weighted_weighted_kernels) {
  weighted_kernels->rowwise() = weights.transpose();
  weighted_kernels->array() *= kernels.array();
  weighted_weighted_kernels->colwise() = weights;
  weighted_weighted_kernels->array() *= weighted_kernels->array();
}

}  // namespace internal

class TooFewClustersLeft : public std::runtime_error {
 public:
  explicit TooFewClustersLeft(const std::string &what_arg)
      : std::runtime_error(what_arg) {}
};

template <typename ConvergencePredicator>
inline ClusterIndicatorMatrix ExecuteKernelKMeans(
    const ClusterIndicatorMatrix &initial_clusters,
    int k_min,
    const KernelMatrix &kernels,
    ConvergencePredicator &converged) {
  if (kernels.rows() != kernels.cols()) {
    throw std::invalid_argument("Kernel matrix must be a square matrix");
  }
  if (initial_clusters.rows() != kernels.rows()) {
    throw std::invalid_argument(
        "Cluster indicator matrix doesn't match for kernel matrix in size");
  }

  ClusterIndicatorMatrix clusters = initial_clusters;
  NormMatrix norms = internal::ComputeNormMatrix(clusters, kernels);
  double prev_norms_total;
  double norms_total = (norms.transpose() * clusters).trace();
  do {
    clusters = internal::ComputeClusterIndicatorMatrix(norms);
    if (clusters.cols() < k_min) {
      throw TooFewClustersLeft("Left clusters are less than wanted");
    }
    norms = internal::ComputeNormMatrix(clusters, kernels);
    prev_norms_total = norms_total;
    norms_total = (norms.transpose() * clusters).trace();
  } while (!converged(prev_norms_total, norms_total));
  return clusters;
}

template <typename ConvergencePredicator>
inline ClusterIndicatorMatrix ExecuteWeightedKernelKMeans(
    const ClusterIndicatorMatrix &initial_clusters,
    int k_min,
    const KernelMatrix &kernels,
    const WeightVector &weights,
    ConvergencePredicator &converged,
    bool less_memory = false) {
  if (kernels.rows() != kernels.cols()) {
    throw std::invalid_argument("Kernel matrix must be a square matrix");
  }
  if (initial_clusters.rows() != kernels.rows()) {
    throw std::invalid_argument(
        "Cluster indicator matrix doesn't match for kernel matrix in size");
  }
  if (weights.size() != initial_clusters.rows()) {
    throw std::invalid_argument("Not enough or too many weights");
  }

  ClusterIndicatorMatrix clusters = initial_clusters;
  double norms_total;
  double prev_norms_total;
  if (less_memory) {
    NormMatrix norms =
        internal::ComputeWeightedNormMatrix(clusters, kernels, weights);
    norms_total = (norms.transpose() * clusters).trace();
    do {
      clusters = internal::ComputeClusterIndicatorMatrix(norms);
      if (clusters.cols() < k_min) {
        throw TooFewClustersLeft("Left clusters are less than wanted");
      }
      norms = internal::ComputeWeightedNormMatrix(clusters, kernels, weights);
      prev_norms_total = norms_total;
      norms_total = (norms.transpose() * clusters).trace();
    } while (!converged(prev_norms_total, norms_total));
  } else {
    KernelMatrix weighted_kernels(kernels.rows(), kernels.cols());
    KernelMatrix weighted_weighted_kernels(kernels.rows(), kernels.cols());
    internal::ComputeWeightedKernelMatrix(
        kernels, weights, &weighted_kernels, &weighted_weighted_kernels);
    NormMatrix norms = internal::ComputeWeightedNormMatrix(
        clusters, kernels, weighted_kernels, weighted_weighted_kernels,
        weights);
    norms_total = (norms.transpose() * clusters).trace();
    do {
      clusters = internal::ComputeClusterIndicatorMatrix(norms);
      if (clusters.cols() < k_min) {
        throw TooFewClustersLeft("Left clusters are less than wanted");
      }
      norms = internal::ComputeWeightedNormMatrix(
          clusters, kernels, weighted_kernels, weighted_weighted_kernels,
          weights);
      prev_norms_total = norms_total;
      norms_total = (norms.transpose() * clusters).trace();
    } while (!converged(prev_norms_total, norms_total));
  }
  
  return clusters;
}

}  // namespace sskkm

#endif  // SSKKM_KERNEL_K_MEANS_H_
