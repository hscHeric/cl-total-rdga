#pragma once
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <functional>
#include <iostream>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using BitSet = boost::dynamic_bitset<>;

class MatrixGraph {
private:
  std::unordered_map<int, BitSet> _adjList;

public:
  explicit MatrixGraph(int n);
  MatrixGraph(const MatrixGraph &) = default;
  ~MatrixGraph() = default;

  [[nodiscard]] int order() const noexcept;
  [[nodiscard]] int size() const noexcept;
  [[nodiscard]] int degree(int v) const;
  [[nodiscard]] std::pair<int, int> degree_range() const noexcept;
  [[nodiscard]] bool contains(int u, int v) const noexcept;
  [[nodiscard]] bool contains(int v) const noexcept;

  void add_edge(int u, int v);
  void remove_edge(int u, int v);
  void remove_vertex(int v);

  [[nodiscard]] std::unordered_set<int> get_vertices() const;
  void print() const;

  [[nodiscard]] const BitSet& get_neighbors(int v) const;
  [[nodiscard]] std::unordered_set<int> get_isolated_vertices() const;
  [[nodiscard]] bool is_isolated_vertex(int vertex) const;
  [[nodiscard]] float get_density() const;

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
