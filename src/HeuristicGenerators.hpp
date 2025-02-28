#pragma once
#include "Chromosome.hpp"
#include "Graph.hpp"

class HeuristicGenerators {
public:
  static Chromosome h1(const Graph &graph);
  static Chromosome h2(const Graph &graph);
  static Chromosome h3(const Graph &graph);
  static Chromosome h4(const Graph &graph);
  static Chromosome h5(const Graph &graph);

private:
  static Graph *copyGraph(const Graph &graph);
  static bool isIsolatedVertex(const Graph &graph, int vertex);
  static std::unordered_set<int> getIsolatedVertices(const Graph &graph);
  static double calculateGraphDensity(const Graph &graph);
};
