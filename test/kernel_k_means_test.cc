#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <vector>

#include "gtest/gtest.h"
#include "sskkm/kernel_k_means.h"
#include "sskkm/kernel_matrix.h"
#include "test/test_util.h"

using namespace sskkm;

namespace {

class Counter {
 public:
  Counter(int count) : count_(count), current_count_(0) {}
  bool operator()(
      double /* previous_total_norms */, double current_total_norms) {
    return ++current_count_ >= count_;
  }
  void Reset() { current_count_ = 0; }

 private:
  Counter(const Counter&);
  Counter& operator=(const Counter&);
  int count_;
  int current_count_;
};

}  // namespace

TEST(KernelKMeansTest, WeightedKkmIsSameAsKkmWhenWeightsAreOnes) {
  DenseMatrix vectors(3, 8);
  vectors <<
      3, 2,  4, 1, 10, 3, 9, 1,
      5, 4,  8, 1,  4, 5, 7, 1,
      7, 6, 16, 1,  6, 7, 4, 1;
  // Simple inner products.
  KernelMatrix kernels = ComputePolynomialKernelMatrix(vectors, 0, 1);
  WeightVector weights = WeightVector::Ones(vectors.cols(), 1);
  ClusterIndicatorMatrix clusters = InitializeRandomClusters(vectors.cols(), 3);
  Counter counter(10);

  ClusterIndicatorMatrix unweighted_result =
      ExecuteKernelKMeans(clusters, clusters.cols(), kernels, counter);

  counter.Reset();
  ClusterIndicatorMatrix weighted_result1 = ExecuteWeightedKernelKMeans(
      clusters, clusters.cols(), kernels, weights, counter);

  counter.Reset();
  ClusterIndicatorMatrix weighted_result2 = ExecuteWeightedKernelKMeans(
      clusters, clusters.cols(), kernels, weights, counter, true);

  EXPECT_TRUE(unweighted_result == weighted_result1);
  EXPECT_TRUE(unweighted_result == weighted_result2);
}

TEST(KernelKMeansTest, TwoImplementationsOfWeightedKkmMakeIdenticalResults) {
  boost::random::mt19937 rng;
  for (int i = 10; i < 20; ++i) {
    boost::random::uniform_int_distribution<> cluster_selector(0, i - 1);
    DenseMatrix vectors = DenseMatrix::Random(3, 500);
    KernelMatrix kernels = ComputePolynomialKernelMatrix(vectors, 0, 2);
    WeightVector weights = WeightVector::Random(vectors.cols());
    ClusterIndicatorMatrix clusters =
        InitializeRandomClusters(vectors.cols(), i);

    Counter counter(10);
    ClusterIndicatorMatrix result1 = ExecuteWeightedKernelKMeans(
        clusters, 1, kernels, weights, counter);

    counter.Reset();
    ClusterIndicatorMatrix result2 = ExecuteWeightedKernelKMeans(
        clusters, 1, kernels, weights, counter, true);

    EXPECT_TRUE(result1 == result2);
  }
}
