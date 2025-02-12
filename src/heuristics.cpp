#include "../include/heuristics.hpp"
#include <algorithm>
#include <random>
#include <set>
#include <vector>

int get_random_int(int ini, int last) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dist(ini, last);
  return dist(gen);
}

Graph *create_graph_copy(const Graph &original) {
  Graph *copy = new BitListGraph(original.order());
  original.for_each_edge([&](int u, int v) { copy->add_edge(u, v); });

  return copy;
}

Chromosome h1(Graph &graph) {
  Chromosome solution(graph.order());

  Graph *h = create_graph_copy(graph);

  std::random_device rd;
  std::mt19937 gen(rd());

  while (!h->get_vertices().empty()) {
    auto vertices = h->get_vertices();
    std::vector<int> vertex_vec(vertices.begin(), vertices.end());
    std::uniform_int_distribution<> dis(0, vertex_vec.size() - 1);
    int v = vertex_vec[dis(gen)]; // Escolhe um vértice aleatório

    // Rotula o vértice escolhido com 2
    solution.set_value(v, 2);

    std::vector<int> neighbors;
    h->for_each_neighbor(v, [&neighbors](int u) { neighbors.push_back(u); });

    if (!neighbors.empty()) {
      // Escolhe o primeiro vizinho e rotula como 1
      int first_neighbor = neighbors[0];
      solution.set_value(first_neighbor, 1);

      // Os demais vizinhos são rotulados com 0
      for (size_t i = 1; i < neighbors.size(); i++) {
        solution.set_value(neighbors[i], 0);
      }
    }

    // Remove vértice e vizinhos
    h->remove_vertex(v);
    for (int neighbor : neighbors) {
      h->remove_vertex(neighbor);
    }

    // Enquanto houver vértices isolados no grafo, rotula cada um com 1
    std::unordered_set<int> isolated_vertices;
    h->for_each_vertex([&](int z) {
      if (h->degree(z) == 0) {
        isolated_vertices.insert(z);
      }
    });

    for (int z : isolated_vertices) {
      solution.set_value(z, 1);

      // Escolhe um primeiro vizinho no grafo original e rotula com 1
      bool set_neighbor = false;
      graph.for_each_neighbor(z, [&](int n) {
        if (!set_neighbor) {
          solution.set_value(n, 1);
          set_neighbor = true;
        }
      });

      h->remove_vertex(z);
    }
  }

  delete h;

  solution.calculate_fitness();
  return solution;
}

Chromosome h5(Graph &graph) { return Chromosome(graph.order(), 1); }
