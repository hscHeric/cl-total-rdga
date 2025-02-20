#include "MatrixGraph.hpp"
#include <algorithm>
#include <stdexcept>

MatrixGraph::MatrixGraph(int n) : Graph(n) {
  // Inicializa o grafo com n vértices numerados de 0 a n-1
  for (int i = 0; i < n; i++) {
    _adjList[i] = BitSet(n);
  }
}

int MatrixGraph::order() const noexcept { return _adjList.size(); }

int MatrixGraph::size() const noexcept {
  int count = 0;
  for (const auto &[_, bitset] : _adjList) {
    count += bitset.count();
  }
  // Cada aresta é contada duas vezes (uma para cada extremidade)
  return count / 2;
}

int MatrixGraph::degree(int v) const {
  if (!contains(v)) {
    throw std::runtime_error("Vertex does not exist");
  }
  return _adjList.at(v).count();
}

std::pair<int, int> MatrixGraph::degree_range() const noexcept {
  if (_adjList.empty()) {
    return {0, 0};
  }

  int min_degree = std::numeric_limits<int>::max();
  int max_degree = 0;

  for (const auto &[_, bitset] : _adjList) {
    int deg = bitset.count();
    min_degree = std::min(min_degree, deg);
    max_degree = std::max(max_degree, deg);
  }

  return {min_degree, max_degree};
}

bool MatrixGraph::contains(int u, int v) const noexcept {
  if (!contains(u) || !contains(v)) {
    return false;
  }

  return _adjList.at(u)[v];
}

bool MatrixGraph::contains(int v) const noexcept {
  return _adjList.find(v) != _adjList.end();
}

void MatrixGraph::add_edge(int u, int v) {
  // Ignora self-loops
  if (u == v) {
    return;
  }

  // Verifica se os vértices existem
  if (!contains(u) || !contains(v)) {
    throw std::runtime_error("Vertex does not exist");
  }

  // Adiciona a aresta (nos dois sentidos, já que é um grafo não direcionado)
  _adjList[u][v] = true;
  _adjList[v][u] = true;
}

void MatrixGraph::remove_edge(int u, int v) {
  if (!contains(u) || !contains(v)) {
    throw std::runtime_error("Vertex does not exist");
  }

  // Remove a aresta (nos dois sentidos)
  _adjList[u][v] = false;
  _adjList[v][u] = false;
}

void MatrixGraph::remove_vertex(int v) {
  if (!contains(v)) {
    throw std::runtime_error("Vertex does not exist");
  }

  // Remove todas as arestas incidentes em v
  for (auto &[vertex, bitset] : _adjList) {
    if (vertex != v) {
      bitset[v] = false;
    }
  }

  // Remove o vértice v
  _adjList.erase(v);
}

void MatrixGraph::for_each_vertex(const VertexCallback &func) const {
  for (const auto &[vertex, _] : _adjList) {
    func(vertex);
  }
}

void MatrixGraph::for_each_edge(const EdgeCallback &func) const {
  std::unordered_set<std::pair<int, int>, PairHash> processed_edges;

  for (const auto &[u, bitset] : _adjList) {
    for (size_t v = 0; v < bitset.size(); ++v) {
      if (bitset[v]) {
        // Garante que cada aresta seja processada apenas uma vez
        auto edge = std::make_pair(std::min(u, static_cast<int>(v)),
                                   std::max(u, static_cast<int>(v)));

        if (processed_edges.find(edge) == processed_edges.end()) {
          func(u, static_cast<int>(v));
          processed_edges.insert(edge);
        }
      }
    }
  }
}

void MatrixGraph::for_each_neighbor(int v, const VertexCallback &func) const {
  if (!contains(v)) {
    throw std::runtime_error("Vertex does not exist");
  }

  const BitSet &bitset = _adjList.at(v);
  for (size_t i = 0; i < bitset.size(); ++i) {
    if (bitset[i]) {
      func(static_cast<int>(i));
    }
  }
}

std::unordered_set<int> MatrixGraph::get_vertices() const {
  std::unordered_set<int> vertices;
  for (const auto &[vertex, _] : _adjList) {
    vertices.insert(vertex);
  }
  return vertices;
}

void MatrixGraph::print() const {
  std::cout << "MatrixGraph: " << order() << " vertices, " << size()
            << " edges\n";

  for (const auto &[u, bitset] : _adjList) {
    std::cout << u << ": ";
    for (size_t v = 0; v < bitset.size(); ++v) {
      if (bitset[v]) {
        std::cout << v << " ";
      }
    }
    std::cout << "\n";
  }
}
