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
  Graph *copy = new BitListGraph(original.order());

  original.for_each_edge([&](int u, int v) { copy->add_edge(u, v); });

  return copy;
}

Chromosome h1(Graph &graph) {
  Chromosome chromosome(graph.order());
  auto temp_graph = create_graph_copy(graph);

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
  Graph *temp_graph = create_graph_copy(graph);

  while (temp_graph->size() > 0) {
    int max_degree_vertex = -1;
    int max_degree = -1;
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

    int max_neighbor_vertex = -1;
    int max_neighbor_degree = -1;
    temp_graph->for_each_neighbor(max_degree_vertex, [&](int neighbor) {
      int degree = graph.degree(neighbor);
      if (degree > max_neighbor_degree) {
        max_neighbor_degree = degree;
        max_neighbor_vertex = neighbor;
      }
    });

    if (max_neighbor_vertex != -1) {
      chromosome.set_value(max_neighbor_vertex, LABEL_ONE);
    }

    std::vector<int> neighbors_to_remove;
    temp_graph->for_each_neighbor(max_degree_vertex, [&](int neighbor) {
      if (neighbor != max_neighbor_vertex) {
        chromosome.set_value(neighbor, LABEL_ZERO);
        neighbors_to_remove.push_back(neighbor);
      }
    });

    for (int neighbor : neighbors_to_remove) {
      temp_graph->remove_vertex(neighbor);
    }
    temp_graph->remove_vertex(max_degree_vertex);

    std::unordered_set<int> isolated_vertices;
    std::vector<int> temp_vertices;
    temp_graph->for_each_vertex([&](int v) {
      if (temp_graph->degree(v) == 0) {
        isolated_vertices.insert(v);
        temp_vertices.push_back(v);
      }
    });

    if (!isolated_vertices.empty()) {
      std::unordered_set<int> neighbors_of_isolated;
      for (int v : isolated_vertices) {
        graph.for_each_neighbor(
            v, [&](int neighbor) { neighbors_of_isolated.insert(neighbor); });
      }

      for (int z : neighbors_of_isolated) {
        int isolated_neighbor_count = 0;
        graph.for_each_neighbor(z, [&](int neighbor) {
          if (isolated_vertices.count(neighbor) > 0) {
            isolated_neighbor_count++;
          }
        });

        if (isolated_neighbor_count >= 2) {
          chromosome.set_value(z, LABEL_TWO);
          graph.for_each_neighbor(z, [&](int neighbor) {
            if (isolated_vertices.count(neighbor) > 0) {
              chromosome.set_value(neighbor, LABEL_ZERO);
            }
          });
        }
      }

      for (int z : neighbors_of_isolated) {
        if (chromosome.get_value(z) == LABEL_ZERO) {
          chromosome.set_value(z, LABEL_TWO);
        }
      }

      for (int v : isolated_vertices) {
        if (chromosome.get_value(v) == LABEL_ZERO) {
          chromosome.set_value(v, LABEL_ZERO);
        }
      }

      for (int v : temp_vertices) {
        temp_graph->remove_vertex(v);
      }
    }
  }

  delete temp_graph;
  chromosome.calculate_fitness();
  return chromosome;
}

Chromosome h5(Graph &graph) { return Chromosome(graph.order(), 1); }
