#include "HeuristicGenerators.hpp"
#include "Chromosome.hpp"
#include "Graph.hpp"
#include "ListGraph.hpp"
#include "MatrixGraph.hpp"
#include <memory>
#include <vector>

Graph *HeuristicGenerators::copyGraph(const Graph &graph) {
  int n = graph.order();
  auto density = 0.0;
  if (n <= 1)
    density = 0.0;

  int m = graph.size();

  density = static_cast<double>(2 * m) / (n * (n - 1));

  const double DENSITY_THRESHOLD = 0.5;

  Graph *graph_copy = nullptr;
  if (density >= DENSITY_THRESHOLD) {
    graph_copy = new MatrixGraph(graph.order());
  } else {
    graph_copy = new ListGraph(graph.order());
  }

  graph.for_each_edge(
      [&graph_copy](int u, int v) { graph_copy->add_edge(u, v); });

  return graph_copy;
}

int HeuristicGenerators::getVertexWithHighestDegree(const Graph &graph) {
  int maxDegreeVertex = -1;
  int maxDegree = -1;

  graph.for_each_vertex([&](int v) {
    if (graph.contains(v)) {
      int degree = graph.degree(v);
      if (degree > maxDegree) {
        maxDegree = degree;
        maxDegreeVertex = v;
      }
    }
  });

  if (maxDegreeVertex == -1) {
    throw std::runtime_error("No vertices found in the graph");
  }

  return maxDegreeVertex;
}

bool HeuristicGenerators::isIsolatedVertex(const Graph &graph, int vertex) {
  return graph.contains(vertex) && graph.degree(vertex) == 0;
}

std::unordered_set<int>
HeuristicGenerators::getIsolatedVertices(const Graph &graph) {
  std::unordered_set<int> isolatedVertices;

  graph.for_each_vertex([&](int v) {
    if (graph.contains(v) && graph.degree(v) == 0) {
      isolatedVertices.insert(v);
    }
  });

  return isolatedVertices;
}

Chromosome HeuristicGenerators::h1(const Graph &graph) {
  Chromosome chromosome(graph.order());

  std::unique_ptr<Graph> h(copyGraph(graph));

  while (h->size() > 0 || !h->get_vertices().empty()) {
    int v = h->choose_rng();

    chromosome.set_value(v, LABEL_TWO);
    std::vector<int> neighbors;
    h->for_each_neighbor(v, [&neighbors](int w) { neighbors.push_back(w); });

    if (!neighbors.empty()) {
      int u = neighbors[0];
      chromosome.set_value(u, LABEL_ONE);

      for (size_t i = 1; i < neighbors.size(); i++) {
        chromosome.set_value(neighbors[i], LABEL_ZERO);
      }
    }

    h->remove_vertex(v);
    for (int neighbor : neighbors) {
      if (h->contains(neighbor)) {
        h->remove_vertex(neighbor);
      }
    }

    auto isolatedVertices = getIsolatedVertices(*h);
    while (!isolatedVertices.empty()) {
      auto z = *isolatedVertices.begin();
      chromosome.set_value(z, LABEL_ONE);

      int x = -1;
      graph.for_each_neighbor(z, [&x](int neighbor) {
        if (x == -1) {
          x = neighbor;
        }
      });

      if (x != -1) {
        chromosome.set_value(x, LABEL_ONE);
      }

      h->remove_vertex(z);
      isolatedVertices = getIsolatedVertices(*h);
    }
  }

  chromosome.calculate_fitness();
  return chromosome;
}
