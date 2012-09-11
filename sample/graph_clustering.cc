#include <boost/foreach.hpp>
/*
  This header is not used directly in this source file but
  |sskkm/ss_kernel_k_means.h|, a header included below, does. And it seems that
  the header needs to be loaded before other Boost.Graph headers (I'm not sure
  the reason :<).
*/
#include <boost/graph/transitive_closure.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <string>

#include "sskkm/convergence_predicator.h"
#include "sskkm/format/graphviz.h"
#include "sskkm/ss_kernel_k_means.h"

using namespace sskkm;
using namespace sskkm::graphviz;

int main(int argc, const char **argv) {
  if (argc != 6) {
    std::cerr << "Usage: " << argv[0]
              << " num_clusters num_least_clusters graph_file must_links_file"
              << " cannot_links_file"
              << std::endl;
    return -1;
  } 

  int k;
  int k_min;
  UndirectedGraph graph;
  MustLinks must_links;
  CannotLinks cannot_links;
  try {
    k = boost::lexical_cast<int>(argv[1]);
    k_min = boost::lexical_cast<int>(argv[2]);
    std::string graph_src = SlurpFile(argv[3]);
    std::string must_links_src = SlurpFile(argv[4]);
    std::string cannot_links_src = SlurpFile(argv[5]);
    graph = ReadUndirectedGraph(graph_src, 1.0);
    must_links = ReadUndirectedGraph(must_links_src, 5.0);
    cannot_links = ReadUndirectedGraph(must_links_src, 5.0);
  } catch (...) {
    std::cerr << "An error occured during input graph load." << std::endl;
    return -1;
  }

  try {
    ConvergencePredicator converged(500, 20, 0.01);
    ClusterIndicatorMatrix clusters = ExecuteSSKernelKMeans(
        k, k_min, kRatioCut, graph, must_links, cannot_links, converged, 10);
    std::vector<ClusterIndicatorMatrix::Index> cluster_numbers(clusters.rows());
    for (ClusterIndicatorMatrix::Index i = 0; i < clusters.cols(); ++i) {
      for (ClusterIndicatorMatrix::InnerIterator iter(clusters, i);
           iter;
           ++iter) {
        cluster_numbers[iter.row()] = iter.col();
      }
    }
    BOOST_FOREACH (
        const ClusterIndicatorMatrix::Index cluster_id, cluster_numbers) {
      std::cout << cluster_id << std::endl;
    }
  } catch (...) {
    std::cerr << "An error occured during clustering." << std::endl;
    return -1;
  }
  return 0;
}
