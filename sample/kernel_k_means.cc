#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "sskkm/convergence_predicator.h"
#include "sskkm/format/cluto.h"
#include "sskkm/kernel_k_means.h"
#include "sskkm/kernel_matrix.h"

using namespace sskkm;
using namespace sskkm::cluto;

namespace {

template <typename Matrix>
KernelMatrix ComputeKernelMatrix(
    const Matrix &vectors, const std::vector<std::string> &kernel_parameters) {
  if (kernel_parameters.size() == 0) {
    throw std::invalid_argument("No kernel parameters are given.");
  }

  KernelMatrix kernels;
  if (kernel_parameters[0] == "gaussian") {
    if (kernel_parameters.size() != 2) {
      std::ostringstream error_message;
      error_message << "Wrong number of kernel parameters; "
                    << "Expected 1, but got " << kernel_parameters.size() - 1
                    << ".";
      throw std::invalid_argument(error_message.str());
    }
    kernels = ComputeGaussianKernelMatrix(
        vectors, boost::lexical_cast<double>(kernel_parameters[1]));
  } else if (kernel_parameters[0] == "polynomial") {
    if (kernel_parameters.size() != 3) {
      std::ostringstream error_message;
      error_message << "Wrong number of kernel parameters; "
                    << "Expected 2, but got " << kernel_parameters.size() - 1
                    << ".";
      throw std::invalid_argument(error_message.str());
    }
    kernels = ComputePolynomialKernelMatrix(
        vectors,
        boost::lexical_cast<double>(kernel_parameters[1]),
        boost::lexical_cast<double>(kernel_parameters[2]));
  } else {
    throw std::invalid_argument(
        "Unknown kernel function: " + kernel_parameters[0] + ".");
  }
  return kernels;
}

void ExitWithHelpMessage(
    int status,
    const boost::program_options::options_description &options_description) {
  (status == 0 ? std::cout : std::cerr) << options_description << std::endl;
  std::exit(status);
}

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

std::vector<std::string> SplitString(
    const std::string &original_string, const std::string &separator) {
  std::vector<std::string> separated_strings;
  std::size_t offset = 0;
  for (std::size_t found_position = original_string.find(separator, offset);
       found_position != std::string::npos;
       found_position = original_string.find(separator, offset)) {
    separated_strings.push_back(
        original_string.substr(offset, found_position - offset));
    offset = found_position + separator.size();
  }
  separated_strings.push_back(original_string.substr(offset));
  return separated_strings;
}

}  // namespace

int main(int argc, const char **argv) {
  // Parse command line options.
  namespace opts = boost::program_options;

  opts::options_description
      clustering_options_description("Clustering options");
  clustering_options_description.add_options()
      ("input-file,i",
       opts::value<std::string>(),
       "Input matrix file (in CLUTO format.)")
      ("kernel,k",
       opts::value<std::string>()->default_value("gaussian,1"),
       "Kernel function and its parameter(s.)")
      ("num-clusters,n",
       opts::value<int>()->required(),
       "The number of clusters.")
      ("num-min-clusters,m",
       opts::value<int>()->default_value(0),
       "The number of clusters wanted at least. "
       "If set to 0, same number as --num-clusters is used.");

  opts::options_description other_options_description("Other options");
  other_options_description.add_options()
      ("help,h", "Show this help message.");

  opts::options_description options_description;
  options_description.add(clustering_options_description);
  options_description.add(other_options_description);

  opts::variables_map clustering_options, other_options;

  opts::store(
      opts::command_line_parser(argc, argv).
      options(other_options_description).allow_unregistered().run(),
      other_options);
  opts::notify(other_options);
  if (!other_options["help"].empty()) {
    ExitWithHelpMessage(0, options_description);
  }

  try {
    opts::store(
        opts::parse_command_line(argc, argv, options_description),
        clustering_options);
    opts::notify(clustering_options);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    ExitWithHelpMessage(-1, options_description);
  }

  // Set up clustering parameters.
  int k = clustering_options["num-clusters"].as<int>();
  int k_min = clustering_options["num-min-clusters"].as<int>() || k;
  std::string cluto_matrix =
      clustering_options["input-file"].empty() ?
      SlurpStream(std::cin) :
      SlurpFile(clustering_options["input-file"].as<std::string>().c_str());
  std::vector<std::string> kernel_parameters =
      SplitString(clustering_options["kernel"].as<std::string>(), ",");
  KernelMatrix kernels;
  ClusterIndicatorMatrix clusters;
  try {
    switch (DetermineMatrixType(cluto_matrix)) {
      case kDenseMatrix: {
        /*
          Usually, vector data are serialized as *row* vectors. So we need to
          transpose it by setting on |transpose| flag as optional 2nd argument.
        */
        DenseMatrix dense_vectors = ParseDenseMatrix(cluto_matrix, true);
        kernels = ComputeKernelMatrix(dense_vectors, kernel_parameters);
        clusters = InitializeRandomClusters(dense_vectors.cols(), k);
        break;
      }
      case kSparseMatrix: {
        // Ditto.
        SparseMatrix sparse_vectors = ParseSparseMatrix(cluto_matrix, true);
        kernels = ComputeKernelMatrix(sparse_vectors, kernel_parameters);
        clusters = InitializeRandomClusters(sparse_vectors.cols(), k);
        break;
      }
      default: {
        throw std::invalid_argument("Input matrix is invalid.");
      }
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    ExitWithHelpMessage(-1, options_description);
  }

  /*
    Execute clustering the vectors into at least k_min clusters. The clustering
    repeats at lest 20 times, and terminates when it exceeds 500 times or eps
    converged within 0.1 %.
  */
  ConvergencePredicator converged(500, 20, 0.001);
  try {
    clusters = ExecuteKernelKMeans(clusters, k_min, kernels, converged);
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
