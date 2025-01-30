
#include "../include/heuristics.hpp"
#include <random>

int get_random_int(int ini, int last) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dist(ini, last);
  return dist(gen);
}

Chromosome h1(DenseGraph &graph) {
  Chromosome chromosome(graph.getVertexCount());

  DenseGraph temp_graph = graph;

  while (temp_graph.getVertexCount() > 0) {
    std::unordered_set<int> unvisited = temp_graph.get_vertices();
    int random_index = get_random_int(0, unvisited.size() - 1);
    auto it = unvisited.begin();
    std::advance(it, random_index);
    chromosome.set_value(*it, LABEL_TWO);

    auto neighbors = temp_graph.get_neighbors_bitset(*it);
    bool first = true;
    for (size_t i = 0; i < neighbors.size(); ++i) {
      if (neighbors[i] == 1) {
        if (!first) {
          chromosome.set_value(i, LABEL_ZERO);
          temp_graph.delVertex(i);
        } else {
          first = false;
          chromosome.set_value(i, LABEL_ONE);
          temp_graph.delVertex(i);
        }
      }
    }
    temp_graph.delVertex(*it);

    // Look for isolated vertices
    unvisited = temp_graph.get_vertices();
    for (const auto &v : unvisited) {
      if (temp_graph.degree(v) == 0) {
        chromosome.set_value(v, LABEL_ONE);
        auto neighbors_of_v_in_G = graph.get_neighbors_bitset(v);
        for (size_t w = 0; w < neighbors_of_v_in_G.size(); ++w) {
          if (neighbors_of_v_in_G[w] == 1) {
            chromosome.set_value(w, LABEL_ONE);
            break;
          }
        }
        temp_graph.delVertex(v);
      }
    }
  }
  chromosome.calculate_fitness();
  return chromosome;
}

Chromosome h2(DenseGraph &graph) {
  Chromosome chromosome(graph.getVertexCount());
  // create a copy of the original graph
  DenseGraph temp_graph = graph;
  // repeat until there are vertices in the auxiliary graph
  while (temp_graph.getVertexCount() > 0) {
    // sort the vertices in descending order of their degrees
    std::list<std::pair<int, int>> pairs =
        temp_graph.get_all_pairs_vertex_degree();

    // select a vertex v with maximum degree
    auto max_degree_vertex_ptr = std::max_element(
        pairs.begin(), pairs.end(),
        [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
          return a.second < b.second;
        });

    int v = (*max_degree_vertex_ptr).first;
    chromosome.set_value(v, LABEL_TWO);

    // Pegue um vizinho u de v e faça f(u) = 1.
    // Para todos os demais vizinhos w de v faça f(w) = 0
    auto neighbors = temp_graph.get_neighbors_bitset(v);
    bool first = true;
    for (size_t i = 0; i < neighbors.size(); ++i) {
      if (neighbors[i] == 1) {
        if (!first) {
          chromosome.set_value(i, LABEL_ZERO);
          temp_graph.delVertex(i);
        } else {
          first = false;
          chromosome.set_value(i, LABEL_ONE);
          temp_graph.delVertex(i);
        }
      }
    }
    temp_graph.delVertex(v);

    // Look for isolated vertices
    auto unvisited = temp_graph.get_vertices();
    for (const auto &v : unvisited) {
      if (temp_graph.degree(v) == 0) {
        chromosome.set_value(v, LABEL_ONE);
        auto neighbors_of_v_in_G = graph.get_neighbors_bitset(v);
        bool has_neighbor_with_label_one = false;
        for (size_t w = 0; w < neighbors_of_v_in_G.size(); ++w) {
          if (neighbors_of_v_in_G[w] == 1 &&
              chromosome.get_value(w) == LABEL_ONE) {
            has_neighbor_with_label_one = true;
            break;
          }
        }
        if (!has_neighbor_with_label_one) {
          for (size_t w = 0; w < neighbors_of_v_in_G.size(); ++w) {
            if (neighbors_of_v_in_G[w] == 1) {
              chromosome.set_value(w, LABEL_ONE);
              break;
            }
          }
        }
        temp_graph.delVertex(v);
      }
    }
  }
  chromosome.calculate_fitness();
  return chromosome;
}

Chromosome h3(DenseGraph &graph) {
  // create an empty cromosome
  Chromosome chromosome(graph.getVertexCount());
  // create a copy of the original graph
  DenseGraph temp_graph = graph;
  // repeat until there are vertices in the auxiliary graph
  while (temp_graph.getVertexCount() > 0) {
    // sort the vertices in descending order of their degrees
    std::list<std::pair<int, int>> pairs =
        temp_graph.get_all_pairs_vertex_degree();

    // select a vertex v with maximum degree
    auto max_degree_vertex_ptr = std::max_element(
        pairs.begin(), pairs.end(),
        [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
          return a.second < b.second;
        });

    int v = (*max_degree_vertex_ptr).first;
    chromosome.set_value(v, LABEL_TWO);

    // Pegue um vizinho u de v e faça f(u) = 1.
    // Para todos os demais vizinhos w de v faça f(w) = 0
    auto neighbors = temp_graph.get_neighbors_bitset(v);
    std::list<std::pair<int, int>> list_pairs;
    for (size_t i = 0; i < neighbors.size(); ++i) {
      if (neighbors[i] == 1) {
        list_pairs.push_back({i, graph.degree(i)});
      }
    }
    // select the neighbor of v with maximum degree
    auto max_degree_vertex_ptr2 = std::max_element(
        list_pairs.begin(), list_pairs.end(),
        [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
          return a.second < b.second;
        });

    int u = (*max_degree_vertex_ptr2).first;
    chromosome.set_value(u, LABEL_ONE);
    for (auto &p : list_pairs) {
      if (p.first != u) {
        chromosome.set_value(p.first, LABEL_ZERO);
        temp_graph.delVertex(p.first);
      }
    }
    temp_graph.delVertex(u);
    temp_graph.delVertex(v);

    // Look for isolated vertices
    auto unvisited = temp_graph.get_vertices();
    for (const auto &v : unvisited) {
      if (temp_graph.degree(v) == 0) {
        chromosome.set_value(v, LABEL_ONE);
        auto neighbors_of_v_in_G = graph.get_neighbors_bitset(v);
        bool has_neighbor_with_label_one = false;
        for (size_t w = 0; w < neighbors_of_v_in_G.size(); ++w) {
          if (neighbors_of_v_in_G[w] == 1 &&
              chromosome.get_value(w) == LABEL_ONE) {
            has_neighbor_with_label_one = true;
            break;
          }
        }
        if (!has_neighbor_with_label_one) {
          int max_degree = 0;
          int degree = 0;
          int old_best_vertex = -1;
          for (size_t w = 0; w < neighbors_of_v_in_G.size(); ++w) {
            if (neighbors_of_v_in_G[w] == 1) {
              degree = graph.degree(w);
              if (degree > max_degree) {
                if (old_best_vertex != -1) {
                  chromosome.set_value(old_best_vertex, LABEL_ZERO);
                }
                max_degree = degree;
                old_best_vertex = w;
                chromosome.set_value(w, LABEL_ONE);
              } else {
                chromosome.set_value(w, LABEL_ZERO);
              }
            }
          }
        }
        temp_graph.delVertex(v);
      }
    }
  }
  chromosome.calculate_fitness();
  return chromosome;
}

Chromosome h5(DenseGraph &graph) {
  return Chromosome(graph.getVertexCount(), 1);
}
