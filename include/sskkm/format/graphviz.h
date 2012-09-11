#ifndef SSKKM_FORMAT_GRAPHVIZ_H_
#define SSKKM_FORMAT_GRAPHVIZ_H_

#include <boost/foreach.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/property_map/dynamic_property_map.hpp>

#include "sskkm/ss_kernel_k_means.h"

namespace sskkm {

namespace graphviz {

namespace internal {

inline void SetDefaultEdgeWeight(
    UndirectedGraph *graph, double default_weight) {
  BOOST_FOREACH (
      const boost::graph_traits<UndirectedGraph>::edge_descriptor edge,
      boost::edges(*graph)) {
    if (boost::get(boost::edge_weight, *graph, edge) == 0.0) {
      boost::put(boost::edge_weight, *graph, edge, default_weight);
    }
  }
}

}  // namespace internal

inline UndirectedGraph ReadUndirectedGraph(
    const std::string &source,
    double default_edge_weight) {
  UndirectedGraph graph;
  boost::dynamic_properties properties;
  properties.property("node_id", boost::get(boost::vertex_name, graph));
  properties.property("weight", boost::get(boost::edge_weight, graph));
  boost::read_graphviz(source, graph, properties);
  internal::SetDefaultEdgeWeight(&graph, default_edge_weight);
  return graph;
}

}  // namespace graphviz

}  // namespace sskkm

#endif  // SSKKM_FORMAT_GRAPHVIZ_H_
