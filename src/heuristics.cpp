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

      temp_graph->for_each_neighbor(max_degree_vertex, [&](int neighbor) {
        if (neighbor != max_neighbor_vertex) {
          chromosome.set_value(neighbor, LABEL_ZERO);
          temp_graph->remove_vertex(neighbor);
        }
      });

      bool has_positive_neighbor = false;
      temp_graph->for_each_neighbor(max_neighbor_vertex, [&](int neighbor) {
        if (chromosome.get_value(neighbor) > 0) {
          has_positive_neighbor = true;
        }
      });

      if (!has_positive_neighbor) {
        temp_graph->for_each_neighbor(max_neighbor_vertex, [&](int neighbor) {
          if (!has_positive_neighbor && neighbor != max_degree_vertex &&
              chromosome.get_value(neighbor) == 0) {
            chromosome.set_value(neighbor, LABEL_ONE);
            has_positive_neighbor = true;
          }
        });
      }

      temp_graph->remove_vertex(max_neighbor_vertex);
    }

    temp_graph->remove_vertex(max_degree_vertex);

    auto vertices = temp_graph->get_vertices();
    std::vector<int> to_process;
    for (int v : vertices) {
      if (temp_graph->degree(v) <= 1) {
        to_process.push_back(v);
      }
    }

    for (int v : to_process) {
      if (temp_graph->contains(v)) {
        bool has_label_two_neighbor = false;
        graph.for_each_neighbor(v, [&](int neighbor) {
          if (chromosome.get_value(neighbor) == LABEL_TWO) {
            has_label_two_neighbor = true;
          }
        });

        if (!has_label_two_neighbor) {
          std::vector<int> uncovered_neighbors;
          graph.for_each_neighbor(v, [&](int neighbor) {
            if (chromosome.get_value(neighbor) == 0) {
              uncovered_neighbors.push_back(neighbor);
            }
          });

          if (!uncovered_neighbors.empty()) {
            int best_neighbor = uncovered_neighbors[0];
            int max_degree = graph.degree(best_neighbor);

            for (size_t i = 1; i < uncovered_neighbors.size(); i++) {
              int curr_degree = graph.degree(uncovered_neighbors[i]);
              if (curr_degree > max_degree) {
                max_degree = curr_degree;
                best_neighbor = uncovered_neighbors[i];
              }
            }

            chromosome.set_value(best_neighbor, LABEL_TWO);
            chromosome.set_value(v, LABEL_ZERO);
          } else {
            chromosome.set_value(v, LABEL_TWO);
          }
        } else {
          chromosome.set_value(v, LABEL_ZERO);
        }

        temp_graph->remove_vertex(v);
      }
    }
  }

  graph.for_each_vertex([&](int v) {
    if (chromosome.get_value(v) == LABEL_ONE) {
      bool has_positive_neighbor = false;
      graph.for_each_neighbor(v, [&](int neighbor) {
        if (chromosome.get_value(neighbor) > 0) {
          has_positive_neighbor = true;
        }
      });

      if (!has_positive_neighbor) {
        chromosome.set_value(v, LABEL_TWO);
      }
    }
  });

  delete temp_graph;
  chromosome.calculate_fitness();
  return chromosome;
}

Chromosome h5(Graph &graph) { return Chromosome(graph.order(), 1); }
