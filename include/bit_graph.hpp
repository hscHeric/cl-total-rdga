#pragma once
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <functional>
#include <iostream>
#include <unordered_set>
#include <vector>

class Graph;
class BitListGraph;
class BitMatrixGraph;

using BitSet = boost::dynamic_bitset<>;
using VertexCallback = std::function<void(int)>;
using EdgeCallback = std::function<void(int, int)>;

class Graph {
protected:
  int vertex_count;
  mutable size_t edge_count{0};

public:
  explicit Graph(int n) noexcept : vertex_count(n) {}
  virtual ~Graph() = default;
  Graph(const Graph &) = delete;
  Graph(Graph &&) = default;

  [[nodiscard]] size_t order() const noexcept { return vertex_count; }
  [[nodiscard]] virtual size_t size() const noexcept = 0;
  [[nodiscard]] virtual int degree(int v) const = 0;
  [[nodiscard]] virtual std::pair<int, int> degree_range() const noexcept = 0;
  [[nodiscard]] virtual bool contains(int u, int v) const noexcept = 0;
  [[nodiscard]] virtual bool contains(int v) const noexcept = 0;

  virtual void add_edge(int u, int v) = 0;
  virtual void remove_edge(int u, int v) = 0;
  virtual void remove_vertex(int v) = 0;

  [[nodiscard]] virtual const BitSet &neighbors(int v) const = 0;
  virtual void for_each_vertex(const VertexCallback &func) const = 0;
  virtual void for_each_edge(const EdgeCallback &func) const = 0;
  virtual void for_each_neighbor(int v, const VertexCallback &func) const = 0;
  [[nodiscard]] virtual std::unordered_set<int> get_vertices() const = 0;
  virtual void print() const = 0;
};

class BitMatrixGraph : public Graph {
private:
  std::vector<BitSet> adjacency_matrix;

public:
  explicit BitMatrixGraph(int n);
  ~BitMatrixGraph() override = default;

  [[nodiscard]] size_t size() const noexcept override;
  [[nodiscard]] int degree(int v) const override;
  [[nodiscard]] std::pair<int, int> degree_range() const noexcept override;
  [[nodiscard]] bool contains(int u, int v) const noexcept override;
  [[nodiscard]] bool contains(int v) const noexcept override;

  void add_edge(int u, int v) override;
  void remove_edge(int u, int v) override;
  void remove_vertex(int v) override;

  [[nodiscard]] const BitSet &neighbors(int v) const override;
  void for_each_vertex(const VertexCallback &func) const override;
  void for_each_edge(const EdgeCallback &func) const override;
  void for_each_neighbor(int v, const VertexCallback &func) const override;
  [[nodiscard]] std::unordered_set<int> get_vertices() const override;
  void print() const override;
};

class BitListGraph : public Graph {
private:
  std::vector<BitSet> adjacency_lists;
  BitSet active_vertices;

public:
  explicit BitListGraph(int n);
  ~BitListGraph() override = default;

  [[nodiscard]] size_t size() const noexcept override;
  [[nodiscard]] int degree(int v) const override;
  [[nodiscard]] std::pair<int, int> degree_range() const noexcept override;
  [[nodiscard]] bool contains(int u, int v) const noexcept override;
  [[nodiscard]] bool contains(int v) const noexcept override;

  void add_edge(int u, int v) override;
  void remove_edge(int u, int v) override;
  void remove_vertex(int v) override;

  [[nodiscard]] const BitSet &neighbors(int v) const override;
  void for_each_vertex(const VertexCallback &func) const override;
  void for_each_edge(const EdgeCallback &func) const override;
  void for_each_neighbor(int v, const VertexCallback &func) const override;
  [[nodiscard]] std::unordered_set<int> get_vertices() const override;
  void print() const override;
};
