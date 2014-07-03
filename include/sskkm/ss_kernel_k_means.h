#ifndef SSKKM_SS_K_MEANS_H_
#define SSKKM_SS_K_MEANS_H_

#include <algorithm>
#include <boost/foreach.hpp>
// XXX: Why transitive_closure.hpp needs to be loaded prior to other BGL
// headers?
#include <boost/graph/transitive_closure.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <boost/property_map/property_map.hpp>
#include <limits>
#include <sstream>
#include <vector>

#include "sskkm/base.h"
#include "sskkm/kernel_k_means.h"

namespace sskkm {

namespace internal {

typedef std::size_t VertexIndex;

typedef boost::property<
  boost::vertex_name_t,
  std::string,
  boost::property<boost::vertex_index_t, VertexIndex> > VertexProperty;

typedef boost::property<boost::edge_weight_t, double> EdgeProperty;

typedef SparseMatrix AdjacencyMatrix;  // aka. Similarity matrix.

typedef SparseMatrix ConstraintPenaltyMatrix;

typedef boost::adjacency_list<
  boost::hash_setS,
  boost::vecS,
  boost::directedS,
  internal::VertexProperty,
  internal::EdgeProperty> DirectedGraph;

}  // namespace internal

typedef boost::adjacency_list<
  boost::hash_setS,
  boost::vecS,
  boost::undirectedS,
  internal::VertexProperty,
  internal::EdgeProperty> UndirectedGraph;

typedef UndirectedGraph MustLinks;

typedef UndirectedGraph CannotLinks;

enum ClusteringObjective { kRatioCut, kRatioAssociation, kNormalizedCut };

class InconsistentConstraints : public std::runtime_error {
 public:
  explicit InconsistentConstraints(const std::string& what_arg)
      : std::runtime_error(what_arg) {}
};

namespace internal {

typedef boost::graph_traits<UndirectedGraph>::vertices_size_type ComponentIndex;

typedef std::map<
  boost::graph_traits<UndirectedGraph>::vertex_descriptor,
  ComponentIndex> ComponentIndices;

typedef SparseMatrix ComponentIndicatorMatrix;

inline ComponentIndicatorMatrix ComputeComponentIndicatorMatrix(
    const ComponentIndices& components) {
  ComponentIndex max_component_index = 0;
  std::vector<SparseMatrixCoefficient> component_indicator_coeffs;
  component_indicator_coeffs.reserve(components.size());
  for (const auto& vertex_component_pair: components) {
    const auto vertex = vertex_component_pair.first;
    const auto component = vertex_component_pair.second;
    if (max_component_index < component) { max_component_index = component; }
    component_indicator_coeffs.push_back(
        SparseMatrixCoefficient(vertex, component, 1.0));
  }
  ComponentIndicatorMatrix component_indicator(
      component_indicator_coeffs.size(), max_component_index + 1);
  component_indicator.setFromTriplets(
      component_indicator_coeffs.begin(), component_indicator_coeffs.end());
  return component_indicator;
}

inline DirectedGraph Undirected2Directed(const UndirectedGraph& undirected) {
  DirectedGraph directed;
  boost::copy_graph(undirected, directed);
  BOOST_FOREACH (
      const boost::graph_traits<DirectedGraph>::edge_descriptor edge,
      boost::edges(directed)) {
    boost::add_edge(
        boost::target(edge, directed),
        boost::source(edge, directed),
        directed);
  }  
  return directed;
}

inline CannotLinks ComputeTransitiveCannotLinks(
    const CannotLinks& given_cannot_links,
    const ComponentIndices& component_indices,
    const ComponentIndicatorMatrix& components) {
  // Lists up component pairs in which there is at least one cannot link across
  // them.
  std::vector<
    std::pair<ComponentIndex, ComponentIndex> > cannot_link_component_pairs;
  BOOST_FOREACH (
      const boost::graph_traits<CannotLinks>::edge_descriptor given_edge,
      boost::edges(given_cannot_links)) {
    ComponentIndex component1 = component_indices.find(
        boost::source(given_edge, given_cannot_links))->second;
    ComponentIndex component2 = component_indices.find(
        boost::target(given_edge, given_cannot_links))->second;
    // TODO(sekia): Inconsistent constraint. Probably should raise an error.
    if (component1 == component2) { continue; }
    if (component1 > component2) { std::swap(component1, component2); }
    cannot_link_component_pairs.push_back(
        std::make_pair(component1, component2));
  }
  std::sort(
      cannot_link_component_pairs.begin(), cannot_link_component_pairs.end());
  std::unique(
      cannot_link_component_pairs.begin(), cannot_link_component_pairs.end());

  UndirectedGraph cannot_links(boost::num_vertices(given_cannot_links));
  for (const auto& cannot_link_component_pair: cannot_link_component_pairs) {
    const auto component1 = cannot_link_component_pair.first;
    const auto component2 = cannot_link_component_pair.second;
    for (ComponentIndicatorMatrix::InnerIterator iter1(components, component1);
         iter1;
         ++iter1) {
      const boost::graph_traits<UndirectedGraph>::vertex_descriptor
          vertex1 = boost::vertex(iter1.row(), cannot_links);
      for (ComponentIndicatorMatrix::InnerIterator
               iter2(components, component2);
           iter2;
           ++iter2) {
        const boost::graph_traits<UndirectedGraph>::vertex_descriptor
            vertex2 = boost::vertex(iter2.row(), cannot_links);
        boost::add_edge(vertex1, vertex2, cannot_links);
      }      
    }
  }
  return cannot_links;
}

inline UndirectedGraph ComputeTransitiveMustLinks(
    const UndirectedGraph& given_must_links) {
  DirectedGraph inferred_must_links;
  boost::transitive_closure(
      Undirected2Directed(given_must_links), inferred_must_links);
  UndirectedGraph must_links;
  boost::copy_graph(inferred_must_links, must_links);
  return must_links;
}

inline AdjacencyMatrix ComputeAdjacencyMatrix(const UndirectedGraph& graph) {
  std::vector<SparseMatrixCoefficient> coeffs;
  coeffs.reserve(2 * boost::num_edges(graph));
  BOOST_FOREACH (
      const boost::graph_traits<UndirectedGraph>::edge_descriptor edge,
      boost::edges(graph)) {
    double weight = boost::get(boost::edge_weight, graph, edge);
    VertexIndex
        source_index = boost::get(
            boost::vertex_index, graph, boost::source(edge, graph)),
        target_index = boost::get(
            boost::vertex_index, graph, boost::target(edge, graph));
    coeffs.push_back(
        SparseMatrixCoefficient(source_index, target_index, weight));
    coeffs.push_back(
        SparseMatrixCoefficient(target_index, source_index, weight));
  }
  AdjacencyMatrix matrix(
      boost::num_vertices(graph), boost::num_vertices(graph));
  matrix.setFromTriplets(coeffs.begin(), coeffs.end());
  return matrix;
}

inline ConstraintPenaltyMatrix ComputeConstraintPenaltyMatrix(
    const MustLinks& must_links,
    const CannotLinks& cannot_links) {
  if (boost::num_vertices(must_links) != boost::num_vertices(cannot_links)) {
    throw std::invalid_argument(
        "The numbers of vertices in must-links and cannot-links must be"
        " equal.");
  }

  std::vector<SparseMatrixCoefficient> coeffs;
  coeffs.reserve(boost::num_edges(must_links) + boost::num_edges(cannot_links));

  const boost::property_map<MustLinks, boost::vertex_index_t>::const_type&
      vertex_indices = boost::get(boost::vertex_index, must_links);

  BOOST_FOREACH (
      const boost::graph_traits<MustLinks>::edge_descriptor edge,
      boost::edges(must_links)) {
    double weight = boost::get(boost::edge_weight, must_links, edge);
    internal::VertexIndex
        source_index = boost::get(
            vertex_indices, boost::source(edge, must_links)),
        target_index = boost::get(
            vertex_indices, boost::target(edge, must_links));
    coeffs.push_back(
        SparseMatrixCoefficient(source_index, target_index, weight));
    coeffs.push_back(
        SparseMatrixCoefficient(target_index, source_index, weight));
  }

  BOOST_FOREACH (
      const boost::graph_traits<CannotLinks>::edge_descriptor edge,
      boost::edges(cannot_links)) {
    double weight = - boost::get(boost::edge_weight, cannot_links, edge);
    internal::VertexIndex
        source_index = boost::get(
            vertex_indices, boost::source(edge, cannot_links)),
        target_index = boost::get(
            vertex_indices, boost::target(edge, cannot_links));
    coeffs.push_back(
        SparseMatrixCoefficient(source_index, target_index, weight));
    coeffs.push_back(
        SparseMatrixCoefficient(target_index, source_index, weight));
  }

  ConstraintPenaltyMatrix penalty_matrix(
      boost::num_vertices(must_links), boost::num_vertices(must_links));
  penalty_matrix.setFromTriplets(coeffs.begin(), coeffs.end());
  return penalty_matrix;
}

inline DenseVector ComputeVertexDegreeVector(const UndirectedGraph& graph) {
  DenseVector degrees = DenseVector::Zero(boost::num_vertices(graph));
  BOOST_FOREACH (
      const boost::graph_traits<UndirectedGraph>::vertex_descriptor vertex,
      boost::vertices(graph)) {
    degrees(boost::get(boost::vertex_index, graph, vertex)) =
        boost::degree(vertex, graph);
  }
  return degrees;
}

inline KernelMatrix ComputeNormalizedCutKernelMatirx(
    const UndirectedGraph& graph,
    const MustLinks& must_links,
    const CannotLinks& cannot_links,
    const WeightVector& degrees,
    double diagonal_shift) {
  DenseMatrix inversed_degrees = degrees.cwiseInverse().asDiagonal();
  KernelMatrix kernel = inversed_degrees *
      (ComputeAdjacencyMatrix(graph) +
       ComputeConstraintPenaltyMatrix(must_links, cannot_links)).toDense() *
      inversed_degrees;
  if (diagonal_shift != 0.0) { kernel += diagonal_shift * inversed_degrees; }
  return kernel;
}

inline KernelMatrix ComputeRatioAssociationKernelMatirx(
    const UndirectedGraph& graph,
    const MustLinks& must_links,
    const CannotLinks& cannot_links,
    double diagonal_shift) {
  KernelMatrix kernel = ComputeAdjacencyMatrix(graph) +
      ComputeConstraintPenaltyMatrix(must_links, cannot_links);
  if (diagonal_shift != 0.0) {
    kernel += diagonal_shift * DenseMatrix::Identity(
        boost::num_vertices(graph), boost::num_vertices(graph));
  }
  return kernel;
}

inline KernelMatrix ComputeRatioCutKernelMatrix(
    const UndirectedGraph& graph,
    const MustLinks& must_links,
    const CannotLinks& cannot_links,
    double diagonal_shift) {
  KernelMatrix kernel =
      (ComputeAdjacencyMatrix(graph) +
       ComputeConstraintPenaltyMatrix(must_links, cannot_links)).toDense() -
      DenseMatrix(ComputeVertexDegreeVector(graph).asDiagonal());
  if (diagonal_shift != 0.0) {
    kernel += diagonal_shift * DenseMatrix::Identity(
        boost::num_vertices(graph), boost::num_vertices(graph));
  }
  return kernel;
}

inline ComponentIndicatorMatrix ComputeCannotLinkedComponents(
    const ComponentIndices& component_indices,
    const ComponentIndicatorMatrix& components,
    const CannotLinks& cannot_links) {
  std::vector<bool> is_cannot_linked_component(components.cols());
  BOOST_FOREACH (
      const boost::graph_traits<CannotLinks>::vertex_descriptor& vertex,
      boost::vertices(cannot_links)) {
    ComponentIndex component = component_indices.find(vertex)->second;
    is_cannot_linked_component[component] = true;
  }
  std::vector<ComponentIndex> candidate_component_indices;
  for (ComponentIndex i = 0; i < components.cols(); ++i) {
    if (is_cannot_linked_component[i]) {
      candidate_component_indices.push_back(i);
    }
  }
  ComponentIndicatorMatrix candidate_components(
      components.rows(), candidate_component_indices.size());
  for (ComponentIndicatorMatrix::Index i = 0;
       i < candidate_components.cols();
       ++i) {
    candidate_components.col(i) =
        components.col(candidate_component_indices[i]);
  }
  return candidate_components;
}

inline ClusterIndicatorMatrix InitializeFarthestFirst(
    std::size_t k,
    const KernelMatrix& kernels,
    const MustLinks& must_links,
    const CannotLinks& cannot_links) {
  MustLinks inferred_must_links = ComputeTransitiveMustLinks(must_links);
  ComponentIndices component_indices;
  boost::connected_components(
      inferred_must_links,
      boost::associative_property_map<ComponentIndices>(component_indices));
  ComponentIndicatorMatrix components =
      ComputeComponentIndicatorMatrix(component_indices);
  // XXX: Maybe unnecessary
  CannotLinks inferred_cannot_links = ComputeTransitiveCannotLinks(
      cannot_links, component_indices, components);

  ComponentIndicatorMatrix candidate_components =
      ComputeCannotLinkedComponents(
          component_indices, components, inferred_cannot_links);
  if (candidate_components.cols() < k) {
    std::ostringstream message;
    message << "Failed to cluster initialize: Cannot set up " << k
            << " cluster(s) such that satisfies given constraints.";
    throw std::runtime_error(message.str());
  }

  // Finds the largest component from candidate components.
  std::vector<bool> is_chosen_component(candidate_components.cols());
  std::vector<ComponentIndex> chosen_component_indices;
  ComponentIndex largest_component = 0;
  for (ComponentIndicatorMatrix::Index i = 1;
       i < candidate_components.cols();
       ++i) {
    if (candidate_components.col(i).nonZeros() >
        candidate_components.col(largest_component).nonZeros()) {
      largest_component = i;
    }
  }
  is_chosen_component[largest_component] = true;
  chosen_component_indices.push_back(largest_component);

  DenseMatrix component_component_similarities =
      candidate_components.transpose() * kernels * candidate_components;
  while (chosen_component_indices.size() < k) {
    ComponentIndicatorMatrix::Index farthest_component = 0;
    double worst_similarity = std::numeric_limits<double>::infinity();
    for (ComponentIndicatorMatrix::Index i = 0;
         i < candidate_components.cols();
         ++i) {
      if (is_chosen_component[i]) { continue; }
      double similarity = 0.0;
      for (const auto j: chosen_component_indices) {
        similarity += component_component_similarities(i, j);
      }
      if (similarity < worst_similarity) {
        farthest_component = i;
        worst_similarity = similarity;
      }
    }
    is_chosen_component[farthest_component] = true;
    chosen_component_indices.push_back(farthest_component);
  }

  ComponentIndicatorMatrix chosen_components(components.rows(), k);
  for (ClusterIndicatorMatrix::Index i = 0; i < chosen_components.cols(); ++i) {
    chosen_components.col(i) =
        candidate_components.col(chosen_component_indices[i]);
  }
  DenseMatrix vertex_component_similarities = kernels * chosen_components;
  std::vector<SparseMatrixCoefficient> cluster_coeffs;
  cluster_coeffs.reserve(vertex_component_similarities.rows());
  for (ClusterIndicatorMatrix::Index i = 0;
       i < vertex_component_similarities.rows();
       ++i) {
    ClusterIndicatorMatrix::Index j;
    vertex_component_similarities.row(i).maxCoeff(&j);
    cluster_coeffs.push_back(SparseMatrixCoefficient(i, j, 1.0));
  }
  ClusterIndicatorMatrix clusters(components.rows(), k);
  clusters.setFromTriplets(cluster_coeffs.begin(), cluster_coeffs.end());
  return clusters;
}

}  // namespace internal

template <typename ConvergencePredicator>
inline ClusterIndicatorMatrix ExecuteSSKernelKMeans(
    int k,
    int k_min,
    ClusteringObjective objective,
    const UndirectedGraph& graph,
    const MustLinks& must_links,
    const CannotLinks& cannot_links,
    ConvergencePredicator& converged,
    double diagonal_shift = 0.0,
    bool less_memory = false) {
  switch (objective) {
    case kNormalizedCut: {
      WeightVector degrees = internal::ComputeVertexDegreeVector(graph);
      KernelMatrix kernels =
          internal::ComputeNormalizedCutKernelMatirx(
              graph, must_links, cannot_links, degrees, diagonal_shift);
      ClusterIndicatorMatrix clusters =
          internal::InitializeFarthestFirst(
              k, kernels, must_links, cannot_links);
      return ExecuteWeightedKernelKMeans(
          clusters, k_min, kernels, degrees, converged, less_memory);
    }
    case kRatioAssociation: {
      KernelMatrix kernels =
          internal::ComputeRatioAssociationKernelMatirx(
              graph, must_links, cannot_links, diagonal_shift);
      ClusterIndicatorMatrix clusters =
          internal::InitializeFarthestFirst(
              k, kernels, must_links, cannot_links);
      return ExecuteKernelKMeans(clusters, k_min, kernels, converged);
    }
    case kRatioCut: {
      KernelMatrix kernels =
          internal::ComputeRatioCutKernelMatrix(
              graph, must_links, cannot_links, diagonal_shift);
      ClusterIndicatorMatrix clusters = internal::InitializeFarthestFirst(
              k, kernels, must_links, cannot_links);
      return ExecuteKernelKMeans(clusters, k_min, kernels, converged);
    }
    default: {
      throw std::invalid_argument("Unknown clustering objective.");
    }
  }
}

}  // namespace sskkm

#endif  // SSKKM_SS_K_MEANS_H_
