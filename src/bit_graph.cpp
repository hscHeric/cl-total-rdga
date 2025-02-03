#include "../include/bit_graph.hpp"
#include <algorithm>
#include <iostream>
#include <stdexcept>

BitMatrixGraph::BitMatrixGraph(int n) : Graph(n) {
  adjacency_matrix.resize(n, BitSet(n));
}

size_t BitMatrixGraph::size() const noexcept { return edge_count; }

int BitMatrixGraph::degree(int v) const {
  if (!contains(v))
    throw std::out_of_range("Vertex not found");
  return adjacency_matrix[v].count();
}

std::pair<int, int> BitMatrixGraph::degree_range() const noexcept {
  int min_degree = vertex_count;
  int max_degree = 0;

  for (int v = 0; v < vertex_count; ++v) {
    if (contains(v)) {
      int d = adjacency_matrix[v].count();
      min_degree = std::min(min_degree, d);
      max_degree = std::max(max_degree, d);
    }
  }
  return {min_degree, max_degree};
}

bool BitMatrixGraph::contains(int u, int v) const noexcept {
  if (u < 0 || u >= vertex_count || v < 0 || v >= vertex_count)
    return false;
  return adjacency_matrix[u].test(v);
}

bool BitMatrixGraph::contains(int v) const noexcept {
  return v >= 0 && v < vertex_count && !adjacency_matrix[v].none();
}

void BitMatrixGraph::add_edge(int u, int v) {
  if (u < 0 || u >= vertex_count || v < 0 || v >= vertex_count)
    throw std::out_of_range("Vertex index out of range");

  if (u == v)
    throw std::invalid_argument("Self-loops are not allowed in simple graphs");

  if (!adjacency_matrix[u].test(v)) {
    adjacency_matrix[u].set(v);
    adjacency_matrix[v].set(u);
    edge_count++;
  }
}

void BitMatrixGraph::remove_edge(int u, int v) {
  if (u < 0 || u >= vertex_count || v < 0 || v >= vertex_count)
    throw std::out_of_range("Vertex index out of range");

  if (adjacency_matrix[u].test(v)) {
    adjacency_matrix[u].reset(v);
    adjacency_matrix[v].reset(u);
    edge_count--;
  }
}

void BitMatrixGraph::remove_vertex(int v) {
  if (!contains(v))
    return;

  edge_count -= adjacency_matrix[v].count();
  adjacency_matrix[v].reset();

  for (int u = 0; u < vertex_count; ++u) {
    if (adjacency_matrix[u].test(v)) {
      adjacency_matrix[u].reset(v);
    }
  }
}

const BitSet &BitMatrixGraph::neighbors(int v) const {
  if (!contains(v))
    throw std::out_of_range("Vertex not found");
  return adjacency_matrix[v];
}

void BitMatrixGraph::for_each_vertex(const VertexCallback &func) const {
  for (int v = 0; v < vertex_count; ++v) {
    if (contains(v))
      func(v);
  }
}

void BitMatrixGraph::for_each_edge(const EdgeCallback &func) const {
  for (int u = 0; u < vertex_count; ++u) {
    for (int v = u + 1; v < vertex_count; ++v) {
      if (contains(u, v))
        func(u, v);
    }
  }
}

void BitMatrixGraph::for_each_neighbor(int v,
                                       const VertexCallback &func) const {
  if (!contains(v))
    throw std::out_of_range("Vertex not found");

  for (size_t u = adjacency_matrix[v].find_first(); u != BitSet::npos;
       u = adjacency_matrix[v].find_next(u)) {
    func(u);
  }
}

std::unordered_set<int> BitMatrixGraph::get_vertices() const {
  std::unordered_set<int> vertices;
  for (int v = 0; v < vertex_count; ++v) {
    if (contains(v))
      vertices.insert(v);
  }
  return vertices;
}

void BitMatrixGraph::print() const {
  for (int v = 0; v < vertex_count; ++v) {
    if (contains(v)) {
      std::cout << v << ": ";
      for_each_neighbor(v, [](int u) { std::cout << u << " "; });
      std::cout << "\n";
    }
  }
}

BitListGraph::BitListGraph(int n) : Graph(n) {
  adjacency_lists.resize(n, BitSet(n));
  active_vertices = BitSet(n);
  active_vertices.set();
}

size_t BitListGraph::size() const noexcept { return edge_count; }

int BitListGraph::degree(int v) const {
  if (!contains(v))
    throw std::out_of_range("Vertex not found");
  return adjacency_lists[v].count();
}

std::pair<int, int> BitListGraph::degree_range() const noexcept {
  int min_degree = vertex_count;
  int max_degree = 0;

  for (size_t v = active_vertices.find_first(); v != BitSet::npos;
       v = active_vertices.find_next(v)) {
    int d = adjacency_lists[v].count();
    min_degree = std::min(min_degree, d);
    max_degree = std::max(max_degree, d);
  }
  return {min_degree, max_degree};
}

bool BitListGraph::contains(int u, int v) const noexcept {
  if (u < 0 || u >= vertex_count || v < 0 || v >= vertex_count)
    return false;
  return active_vertices.test(u) && adjacency_lists[u].test(v);
}

bool BitListGraph::contains(int v) const noexcept {
  return v >= 0 && v < vertex_count && active_vertices.test(v);
}

void BitListGraph::add_edge(int u, int v) {
  if (u < 0 || u >= vertex_count || v < 0 || v >= vertex_count)
    throw std::out_of_range("Vertex index out of range");

  if (u == v)
    throw std::invalid_argument("Self-loops are not allowed in simple graphs");

  if (!adjacency_lists[u].test(v)) {
    adjacency_lists[u].set(v);
    adjacency_lists[v].set(u);
    edge_count++;
  }
}

void BitListGraph::remove_edge(int u, int v) {
  if (u < 0 || u >= vertex_count || v < 0 || v >= vertex_count)
    throw std::out_of_range("Vertex index out of range");

  if (adjacency_lists[u].test(v)) {
    adjacency_lists[u].reset(v);
    adjacency_lists[v].reset(u);
    edge_count--;
  }
}

void BitListGraph::remove_vertex(int v) {
  if (!contains(v))
    return;

  for (size_t u = adjacency_lists[v].find_first(); u != BitSet::npos;
       u = adjacency_lists[v].find_next(u)) {
    adjacency_lists[u].reset(v);
    edge_count--;
  }

  adjacency_lists[v].reset();
  active_vertices.reset(v);
}

const BitSet &BitListGraph::neighbors(int v) const {
  if (!contains(v))
    throw std::out_of_range("Vertex not found");
  return adjacency_lists[v];
}

void BitListGraph::for_each_vertex(const VertexCallback &func) const {
  for (size_t v = active_vertices.find_first(); v != BitSet::npos;
       v = active_vertices.find_next(v)) {
    func(v);
  }
}

void BitListGraph::for_each_edge(const EdgeCallback &func) const {
  for (size_t u = active_vertices.find_first(); u != BitSet::npos;
       u = active_vertices.find_next(u)) {
    for (size_t v = adjacency_lists[u].find_first(); v != BitSet::npos;
         v = adjacency_lists[u].find_next(v)) {
      if (v > u && active_vertices.test(v))
        func(u, v);
    }
  }
}

void BitListGraph::for_each_neighbor(int v, const VertexCallback &func) const {
  if (!contains(v))
    throw std::out_of_range("Vertex not found");

  for (size_t u = adjacency_lists[v].find_first(); u != BitSet::npos;
       u = adjacency_lists[v].find_next(u)) {
    if (active_vertices.test(u))
      func(u);
  }
}

std::unordered_set<int> BitListGraph::get_vertices() const {
  std::unordered_set<int> vertices;
  for (size_t v = active_vertices.find_first(); v != BitSet::npos;
       v = active_vertices.find_next(v)) {
    vertices.insert(v);
  }
  return vertices;
}

void BitListGraph::print() const {
  for (size_t v = active_vertices.find_first(); v != BitSet::npos;
       v = active_vertices.find_next(v)) {
    std::cout << v << ": ";
    for_each_neighbor(v, [](int u) { std::cout << u << " "; });
    std::cout << "\n";
  }
}
