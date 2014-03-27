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
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

#include "sskkm/convergence_predicator.h"
#include "sskkm/format/graphviz.h"
#include "sskkm/ss_kernel_k_means.h"

using namespace sskkm;
using namespace sskkm::graphviz;

namespace {

void ExitWithHelpMessage(
    int status,
    const boost::program_options::options_description &options_description) {
  (status == 0 ? std::cout : std::cerr) << options_description << std::endl;
  std::exit(status);
}

ClusteringObjective DetermineObjectiveFunction(
    const std::string &objective_name) {
  if (objective_name == "ratio_cut") {
    return kRatioCut;
  } else if (objective_name == "ratio_association") {
    return kRatioAssociation;
  } else if (objective_name == "normalized_cut") {
    return kNormalizedCut;
  } else {
    throw std::invalid_argument(
        "Not supported objective function: " + objective_name);
  }
}

}  // namespace

int main(int argc, const char **argv) {
  namespace opts = boost::program_options;

  opts::options_description
      clustering_options_description("Clustering options");
  clustering_options_description.add_options()
      ("objective,o",
       opts::value<std::string>()->default_value("ratio_cut"),
       "Objective function to be minimize/maximize.")
      ("input-file,i",
       opts::value<std::string>()->required(),
       "Input graph file (in graphviz format.)")
      ("must-link-file,M",
       opts::value<std::string>(),
       "Must link constraint graph file.")
      ("cannot-link-file,c",
       opts::value<std::string>(),
       "Cannot link constraint graph file.")
      ("default-penalty-weight,p",
       opts::value<double>()->default_value(5.0),
       "Default penalty given when the constrants are violated.")
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

  int k = clustering_options["num-clusters"].as<int>();
  int k_min = clustering_options["num-min-clusters"].as<int>() || k;
  double default_penalty =
      clustering_options["default-penalty-weight"].as<double>();
  ClusteringObjective objective;
  UndirectedGraph graph;
  MustLinks must_links;
  CannotLinks cannot_links;
  try {
    objective = DetermineObjectiveFunction(
        clustering_options["objective"].as<std::string>());

    std::string graph_source =
        SlurpFile(clustering_options["input-file"].as<std::string>().c_str());
    graph = ReadUndirectedGraph(graph_source, 1.0);

    if (!clustering_options["must-link-file"].empty()) {
      std::string must_links_source =
          SlurpFile(
              clustering_options["must-link-file"].as<std::string>().c_str());
      must_links = ReadUndirectedGraph(must_links_source, default_penalty);
    } else {
      must_links = MustLinks(boost::num_vertices(graph));
    }

    if (!clustering_options["cannot-link-file"].empty()) {
      std::string cannot_links_source =
          SlurpFile(
              clustering_options["cannot-link-file"].as<std::string>().c_str());
      cannot_links = ReadUndirectedGraph(cannot_links_source, default_penalty);
    } else {
      cannot_links = CannotLinks(boost::num_vertices(graph));
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return -1;
  }

  try {
    ConvergencePredicator converged(500, 20, 0.01);
    ClusterIndicatorMatrix clusters = ExecuteSSKernelKMeans(
        k, k_min, objective, graph, must_links, cannot_links, converged, 10);
    std::vector<ClusterIndicatorMatrix::Index> cluster_numbers(clusters.rows());
    for (ClusterIndicatorMatrix::Index i = 0; i < clusters.cols(); ++i) {
      for (ClusterIndicatorMatrix::InnerIterator iter(clusters, i);
           iter;
           ++iter) {
        cluster_numbers[iter.row()] = iter.col();
      }
    }
    for (const auto cluster_id: cluster_numbers) {
      std::cout << cluster_id << std::endl;
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return -1;
  }
  return 0;
}
