#include <boost/foreach.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <cstdlib>  // for std::atai()
#include <iostream>
#include <string>
#include <vector>

#include "sskkm/convergence_predicator.h"
#include "sskkm/format/cluto.h"
#include "sskkm/kernel_k_means.h"
#include "sskkm/kernel_matrix.h"

using namespace sskkm;
using namespace sskkm::cluto;

namespace {

boost::random::mt19937 rng_;

// TODO(sekia): Move to library.
ClusterIndicatorMatrix InitializeRandomClusters(
    ClusterIndicatorMatrix::Index num_vectors,
    ClusterIndicatorMatrix::Index num_clusters) {
  boost::random::uniform_int_distribution<> dist(0, num_clusters - 1);

  std::vector<SparseMatrixCoefficient> coeffs(num_vectors);
  std::vector<int> cluster_sizes(num_clusters);
  while (true) {
    for (ClusterIndicatorMatrix::Index i = 0; i < num_vectors; ++i) {
      ClusterIndicatorMatrix::Index cluster_index = dist(rng_);
      coeffs[i] = SparseMatrixCoefficient(i, cluster_index, 1.0);
      ++cluster_sizes[cluster_index];
    }
    BOOST_FOREACH (const int cluster_size, cluster_sizes) {
      if (cluster_size == 0) { goto CONTINUE_OUTER_LOOP; }
    }
    break;
 CONTINUE_OUTER_LOOP:
    cluster_sizes.clear();
    cluster_sizes.resize(num_clusters);
  }

  ClusterIndicatorMatrix clusters(num_vectors, num_clusters);
  clusters.setFromTriplets(coeffs.begin(), coeffs.end());
  return clusters;
}

}  // namespace

int main(int argc, const char **argv) {
  if (argc != 3) {
    std::cerr << "Usage: kernel_k_means num_clusters cluto_matrix_file.mat"
              << std::endl;
    return -1;
  }

  int k = std::atoi(argv[1]);
  std::string cluto_matrix = SlurpFile(argv[2]);
  KernelMatrix kernels;
  ClusterIndicatorMatrix clusters;
  switch (DetermineMatrixType(cluto_matrix)) {
    case kDenseMatrix: {
      /*
        Usually, vector data are serialized as *row* vectors. So we need to
        transpose it by setting on |transpose| flag as optional 2nd argument.
      */
      DenseMatrix dense_vectors = ParseDenseMatrix(cluto_matrix, true);
      kernels = ComputeGaussianKernelMatrix(dense_vectors, 1);
      clusters = InitializeRandomClusters(dense_vectors.cols(), k);
      break;
    }
    case kSparseMatrix: {
      // Ditto.
      SparseMatrix sparse_vectors = ParseSparseMatrix(cluto_matrix, true);
      kernels = ComputePolynomialKernelMatrix(sparse_vectors, 0, 2);
      clusters = InitializeRandomClusters(sparse_vectors.cols(), k);
      break;
    }
    default:
      std::cerr << "Invalid input" << std::endl;
      return -1;
  }

  /*
    Classifies the vectors into 6 clusters. The clustering repeats at lest 20
    times, and terminates when it exceeds 500 times or eps converged within
    0.1 %.
  */
  ConvergencePredicator converged(500, 20, 0.001);
  try {
    clusters = ExecuteKernelKMeans(clusters, k, kernels, converged);
  } catch (const TooFewClustersLeft &) {
    std::cerr << "Number of clusters became less than wanted." << std::endl;
    return 1;
  }

  for (ClusterIndicatorMatrix::Index i = 0; i < clusters.cols(); ++i) {
    std::cerr << "Cluster #" << i << ": " << clusters.col(i).nonZeros()
              << " vector(s)" << std::endl;
  }

  std::vector<ClusterIndicatorMatrix::Index> cluster_numbers(clusters.rows());
  for (ClusterIndicatorMatrix::Index i = 0; i < clusters.cols(); ++i) {
    for (SparseMatrix::InnerIterator iter(clusters, i); iter; ++iter) {
      cluster_numbers[iter.row()] = iter.col();
    }
  }
  BOOST_FOREACH (
      const ClusterIndicatorMatrix::Index cluster_id, cluster_numbers) {
    std::cout << cluster_id << std::endl;
  }

  return 0;
}
