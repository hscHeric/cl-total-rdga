#pragma once
#include <functional>
#include <iostream>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class ListGraph {
private:
  std::unordered_map<int, std::vector<int>> _adjList;

public:
  explicit ListGraph(int n);
  ~ListGraph() = default;

  [[nodiscard]] int order() const noexcept;
  [[nodiscard]] int size() const noexcept;
  [[nodiscard]] int degree(int v) const;
  [[nodiscard]] std::pair<int, int> degree_range() const noexcept;
  [[nodiscard]] bool contains(int u, int v) const noexcept;
  [[nodiscard]] bool contains(int v) const noexcept;

  void add_edge(int u, int v);
  void remove_edge(int u, int v);
  void remove_vertex(int v);

  [[nodiscard]] const std::vector<int>& get_neighbors(int v) const;

  [[nodiscard]] std::unordered_set<int> get_vertices() const;
  [[nodiscard]] std::unordered_set<int> get_isolated_vertices() const;
  [[nodiscard]] bool is_isolated_vertex(int vertex) const;
  [[nodiscard]] float get_density() const;
  void print() const;

  [[nodiscard]] int choose_rng() const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::unordered_set<int> vertices = get_vertices();

    if (vertices.empty()) {
      throw std::runtime_error("Graph is empty, cannot choose a vertex.");
    }
    std::uniform_int_distribution<> distrib(0, vertices.size() - 1);
    auto it = vertices.begin();
    std::advance(it, distrib(gen));
    return *it;
  };
};
