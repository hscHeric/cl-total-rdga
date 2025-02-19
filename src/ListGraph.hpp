#pragma once
#include "Graph.hpp"

class ListGraph : public Graph {
private:
  std::unordered_map<int, std::vector<int>> _adjList;

public:
  explicit ListGraph(int n);
  ~ListGraph() override = default;

  [[nodiscard]] int order() const noexcept override;
  [[nodiscard]] int size() const noexcept override;
  [[nodiscard]] int degree(int v) const override;
  [[nodiscard]] std::pair<int, int> degree_range() const noexcept override;
  [[nodiscard]] bool contains(int u, int v) const noexcept override;
  [[nodiscard]] bool contains(int v) const noexcept override;

  void add_edge(int u, int v) override;
  void remove_edge(int u, int v) override;
  void remove_vertex(int v) override;

  void for_each_vertex(const VertexCallback &func) const override;
  void for_each_edge(const EdgeCallback &func) const override;
  void for_each_neighbor(int v, const VertexCallback &func) const override;
  [[nodiscard]] std::unordered_set<int> get_vertices() const override;
  void print() const override;
};
