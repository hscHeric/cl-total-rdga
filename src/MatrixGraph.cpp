#include "MatrixGraph.hpp"

MatrixGraph::MatrixGraph(int n) : Graph(n) {
  for (int i = 0; i < n; ++i) {
    _adjList[i];
    _adjList[i].resize(n);
  }
}

int MatrixGraph::order() const noexcept {
  return static_cast<int>(_adjList.size());
}

int MatrixGraph::size() const noexcept {
  int edgeCount = 0;
  for (const auto &p : _adjList) {
    edgeCount += p.second.count();
  }
  return edgeCount / 2;
}

int MatrixGraph::degree(int v) const {
  // já lança uma exceção se não existir o vértice v no std::unordered_map
  return _adjList.at(v).count();
}

std::pair<int, int> MatrixGraph::degree_range() const noexcept {
  int min_degree = order();
  int max_degree = 0;

  for (const auto &p : _adjList) {
    int current_degree = degree(p.first);
    min_degree = std::min(min_degree, current_degree);
    max_degree = std::max(max_degree, current_degree);
  }
  return {min_degree, max_degree};
}

bool MatrixGraph::contains(int u, int v) const noexcept {
  if (contains(u) && contains(v)) {
    // já lança uma exceção se não tiver o vértice u
    return (_adjList.at(u)[v] == 1) ? true : false;
  } else {
    return false;
  }
}

bool MatrixGraph::contains(int v) const noexcept {
  return _adjList.count(v) != 0;
}

void MatrixGraph::add_edge(int u, int v) {
  if (!contains(u) || !contains(v)) {
    throw std::out_of_range("Endpoints of edge are not in the graph");
  }

  if (u == v) {
    throw std::invalid_argument("Self-loops are not allowed in simple graphs");
  }

  _adjList[u][v] = 1;
  _adjList[v][u] = 1;
}

void MatrixGraph::remove_edge(int u, int v) {
  // já lança exceção se não tiver u ou v
  _adjList.at(u)[v] = 0;
  _adjList.at(v)[u] = 0;
}

void MatrixGraph::remove_vertex(int v) {
  if (contains(v)) {
    for (unsigned i = 0; i < _adjList[v].size(); ++i) {
      if (_adjList[v][i] == 1) {
        _adjList.at(i)[v] = 0;
      }
    }
    _adjList.erase(v);
  }
}

void MatrixGraph::for_each_vertex(const VertexCallback &func) const {
  for (const auto &p : _adjList) {
    func(p.first);
  }
}

void MatrixGraph::for_each_edge(const EdgeCallback &func) const {
  for (const auto &par : _adjList) {
    for (size_t v = par.first + 1; v < par.second.size(); ++v) {
      if (contains(par.first, v)) {
        func(par.first, v);
      }
    }
  }
}

void MatrixGraph::for_each_neighbor(int v, const VertexCallback &func) const {
  if (!contains(v)) {
    throw std::out_of_range("Vertex not found");
  }

  for (size_t u = _adjList.at(v).find_first(); u != BitSet::npos;
       u = _adjList.at(v).find_next(u)) {
    func(u);
  }
}

std::unordered_set<int> MatrixGraph::get_vertices() const {
  std::unordered_set<int> vertices;
  for (const auto &p : _adjList) {
    vertices.insert(p.first);
  }
  return vertices;
}

void MatrixGraph::print() const {
  for (auto it = _adjList.cbegin(); it != _adjList.cend(); ++it) {
    std::cout << it->first << ": ";
    std::cout << it->second;
    std::cout << '\n';
  }
}
