#include "../include/heuristics.hpp"
#include <algorithm>
#include <random>

int get_random_int(int ini, int last) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dist(ini, last);
  return dist(gen);
}

Graph *create_graph_copy(const Graph &original) {
  Graph *copy = new BitListGraph(
      original.order()); // Using BitListGraph as default implementation

  original.for_each_edge([&](int u, int v) { copy->add_edge(u, v); });

  return copy;
}

Chromosome h1(Graph &graph) {
  Chromosome chromosome(graph.order());
  auto temp_graph =
      create_graph_copy(graph); // Assume we have this utility function

  while (temp_graph->size() > 0) {
    auto unvisited = temp_graph->get_vertices();
    if (unvisited.empty())
      break;

    int random_index = get_random_int(0, unvisited.size() - 1);
    auto it = unvisited.begin();
    std::advance(it, random_index);
    int selected_vertex = *it;
    chromosome.set_value(selected_vertex, LABEL_TWO);

    bool first = true;
    temp_graph->for_each_neighbor(selected_vertex, [&](int neighbor) {
      if (!first) {
        chromosome.set_value(neighbor, LABEL_ZERO);
        temp_graph->remove_vertex(neighbor);
      } else {
        first = false;
        chromosome.set_value(neighbor, LABEL_ONE);
        temp_graph->remove_vertex(neighbor);
      }
    });

    temp_graph->remove_vertex(selected_vertex);

    // Handle isolated vertices
    unvisited = temp_graph->get_vertices();
    for (int v : unvisited) {
      if (temp_graph->degree(v) == 0) {
        chromosome.set_value(v, LABEL_ONE);
        bool found_neighbor = false;
        graph.for_each_neighbor(v, [&](int w) {
          if (!found_neighbor) {
            chromosome.set_value(w, LABEL_ONE);
            found_neighbor = true;
          }
        });
        temp_graph->remove_vertex(v);
      }
    }
  }

  delete temp_graph;
  chromosome.calculate_fitness();
  return chromosome;
}

Chromosome h2(Graph &graph) {
  Chromosome chromosome(graph.order());
  Graph *temp_graph = create_graph_copy(graph);

  while (temp_graph->size() > 0) {
    // Find vertex with maximum degree
    int max_degree = -1;
    int max_degree_vertex = -1;
    temp_graph->for_each_vertex([&](int v) {
      int degree = temp_graph->degree(v);
      if (degree > max_degree) {
        max_degree = degree;
        max_degree_vertex = v;
      }
    });

    if (max_degree_vertex == -1)
      break;

    chromosome.set_value(max_degree_vertex, LABEL_TWO);
    bool first = true;

    temp_graph->for_each_neighbor(max_degree_vertex, [&](int neighbor) {
      if (!first) {
        chromosome.set_value(neighbor, LABEL_ZERO);
        temp_graph->remove_vertex(neighbor);
      } else {
        first = false;
        chromosome.set_value(neighbor, LABEL_ONE);
        temp_graph->remove_vertex(neighbor);
      }
    });

    temp_graph->remove_vertex(max_degree_vertex);

    // Handle isolated vertices
    auto unvisited = temp_graph->get_vertices();
    for (int v : unvisited) {
      if (temp_graph->degree(v) == 0) {
        chromosome.set_value(v, LABEL_ONE);
        bool has_neighbor_with_label_one = false;

        graph.for_each_neighbor(v, [&](int w) {
          if (chromosome.get_value(w) == LABEL_ONE) {
            has_neighbor_with_label_one = true;
          }
        });

        if (!has_neighbor_with_label_one) {
          bool found_neighbor = false;
          graph.for_each_neighbor(v, [&](int w) {
            if (!found_neighbor) {
              chromosome.set_value(w, LABEL_ONE);
              found_neighbor = true;
            }
          });
        }
        temp_graph->remove_vertex(v);
      }
    }
  }

  delete temp_graph;
  chromosome.calculate_fitness();
  return chromosome;
}

Chromosome h3(Graph &graph) {
  Chromosome chromosome(graph.order());
  Graph *temp_graph = create_graph_copy(graph);

  while (temp_graph->size() > 0) {
    // Find vertex with maximum degree
    int max_degree = -1;
    int max_degree_vertex = -1;
    temp_graph->for_each_vertex([&](int v) {
      int degree = temp_graph->degree(v);
      if (degree > max_degree) {
        max_degree = degree;
        max_degree_vertex = v;
      }
    });

    if (max_degree_vertex == -1)
      break;

    chromosome.set_value(max_degree_vertex, LABEL_TWO);

    // Find neighbor with maximum degree
    struct NeighborInfo {
      int vertex;
      int degree;
    };
    std::vector<NeighborInfo> neighbors;

    temp_graph->for_each_neighbor(max_degree_vertex, [&](int neighbor) {
      neighbors.push_back({neighbor, graph.degree(neighbor)});
    });

    if (!neighbors.empty()) {
      auto max_neighbor =
          std::max_element(neighbors.begin(), neighbors.end(),
                           [](const NeighborInfo &a, const NeighborInfo &b) {
                             return a.degree < b.degree;
                           });

      chromosome.set_value(max_neighbor->vertex, LABEL_ONE);

      for (const auto &info : neighbors) {
        if (info.vertex != max_neighbor->vertex) {
          chromosome.set_value(info.vertex, LABEL_ZERO);
          temp_graph->remove_vertex(info.vertex);
        }
      }
      temp_graph->remove_vertex(max_neighbor->vertex);
    }

    temp_graph->remove_vertex(max_degree_vertex);

    // Handle isolated vertices
    auto unvisited = temp_graph->get_vertices();
    for (int v : unvisited) {
      if (temp_graph->degree(v) == 0) {
        chromosome.set_value(v, LABEL_ONE);
        bool has_neighbor_with_label_one = false;

        graph.for_each_neighbor(v, [&](int w) {
          if (chromosome.get_value(w) == LABEL_ONE) {
            has_neighbor_with_label_one = true;
          }
        });

        if (!has_neighbor_with_label_one) {
          int max_degree = -1;
          int best_vertex = -1;

          graph.for_each_neighbor(v, [&](int w) {
            int degree = graph.degree(w);
            if (degree > max_degree) {
              if (best_vertex != -1) {
                chromosome.set_value(best_vertex, LABEL_ZERO);
              }
              max_degree = degree;
              best_vertex = w;
              chromosome.set_value(w, LABEL_ONE);
            } else {
              chromosome.set_value(w, LABEL_ZERO);
            }
          });
        }
        temp_graph->remove_vertex(v);
      }
    }
  }

  delete temp_graph;
  chromosome.calculate_fitness();
  return chromosome;
}
Chromosome h4(Graph &graph) {
  Chromosome chromosome(graph.order());
  Graph *temp_graph =
      create_graph_copy(graph); // Cria uma cópia do grafo original

  while (temp_graph->size() > 0) {
    // Passo 3: Encontre o vértice com o maior grau
    int max_degree = -1;
    int max_degree_vertex = -1;
    temp_graph->for_each_vertex([&](int v) {
      int degree = temp_graph->degree(v);
      if (degree > max_degree) {
        max_degree = degree;
        max_degree_vertex = v;
      }
    });

    if (max_degree_vertex == -1)
      break; // Se não há mais vértices, saia do loop

    // Passo 4: Atribua rótulo 2 ao vértice de maior grau
    chromosome.set_value(max_degree_vertex, LABEL_TWO);

    // Passo 5: Encontre o vizinho de maior grau e atribua rótulo 1
    int max_neighbor_degree = -1;
    int max_neighbor_vertex = -1;
    temp_graph->for_each_neighbor(max_degree_vertex, [&](int neighbor) {
      int degree = temp_graph->degree(neighbor);
      if (degree > max_neighbor_degree) {
        max_neighbor_degree = degree;
        max_neighbor_vertex = neighbor;
      }
    });

    if (max_neighbor_vertex != -1) {
      chromosome.set_value(max_neighbor_vertex, LABEL_ONE);
    }

    // Passo 6: Atribua rótulo 0 aos demais vizinhos
    temp_graph->for_each_neighbor(max_degree_vertex, [&](int neighbor) {
      if (neighbor != max_neighbor_vertex) {
        chromosome.set_value(neighbor, LABEL_ZERO);
      }
    });

    // Passo 7: Remova o vértice e seus vizinhos do grafo temporário
    temp_graph->remove_vertex(max_degree_vertex);
    temp_graph->for_each_neighbor(max_degree_vertex, [&](int neighbor) {
      temp_graph->remove_vertex(neighbor);
    });

    // Passo 8: Verifique vértices isolados
    auto unvisited = temp_graph->get_vertices();
    std::unordered_set<int> isolated_vertices;

    for (int v : unvisited) {
      if (temp_graph->degree(v) == 0) {
        isolated_vertices.insert(v);
      }
    }

    if (!isolated_vertices.empty()) {
      // Passo 9: Seja S o conjunto de vértices isolados
      std::unordered_set<int> S = isolated_vertices;

      // Passo 10: Seja N(S) o conjunto de vizinhos de S no grafo original
      std::unordered_set<int> N_S;
      for (int v : S) {
        graph.for_each_neighbor(v, [&](int neighbor) { N_S.insert(neighbor); });
      }

      // Passo 11: Para cada z em N(S) com pelo menos dois vizinhos em S
      for (int z : N_S) {
        int count = 0;
        graph.for_each_neighbor(z, [&](int neighbor) {
          if (S.count(neighbor)) {
            count++;
          }
        });

        if (count >= 2) {
          chromosome.set_value(z, LABEL_TWO);
          graph.for_each_neighbor(z, [&](int neighbor) {
            if (S.count(neighbor)) {
              chromosome.set_value(neighbor, LABEL_ZERO);
            }
          });
        }
      }

      // Passo 12: Pinte os demais vértices em N(S) com rótulo 2
      for (int z : N_S) {
        if (chromosome.get_value(z) == LABEL_UNDEFINED) {
          chromosome.set_value(z, LABEL_TWO);
        }
      }

      // Passo 13: Pinte com rótulo 0 os vértices em S que não foram pintados
      for (int v : S) {
        if (chromosome.get_value(v) == LABEL_UNDEFINED) {
          chromosome.set_value(v, LABEL_ZERO);
        }
      }

      // Passo 14: Remova todos os vértices de S do grafo temporário
      for (int v : S) {
        temp_graph->remove_vertex(v);
      }
    }
  }

  delete temp_graph;
  chromosome.calculate_fitness();
  return chromosome;
}

Chromosome h5(Graph &graph) { return Chromosome(graph.order(), 1); }
