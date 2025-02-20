#pragma once
#include "Graph.hpp"

// Estrutura auxiliar para hash de pares

class MatrixGraph : public Graph {
private:
  std::unordered_map<int, BitSet> _adjList;

public:
  explicit MatrixGraph(int n);
  MatrixGraph(const MatrixGraph &) = default;
  ~MatrixGraph() override = default;

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
