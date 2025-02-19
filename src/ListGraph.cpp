#include "ListGraph.hpp"

ListGraph::ListGraph(int n) : Graph(n) {
  for (int i = 0; i < n; ++i) {
    _adjList[i];
  }
}

int ListGraph::order() const noexcept { return _adjList.size(); }

int ListGraph::size() const noexcept {
  size_t total = 0;
  for (const auto &p : _adjList) {
    total += p.second.size();
  }
  return static_cast<int>(total) / 2;
}

int ListGraph::degree(int v) const {
  // já lança uma exceção se não tiver o vértice v
  return _adjList.at(v).size();
}

std::pair<int, int> ListGraph::degree_range() const noexcept {
  int min_degree = order();
  int max_degree = 0;

  for (const auto &p : _adjList) {
    int current_degree = degree(p.first);
    min_degree = std::min(min_degree, current_degree);
    max_degree = std::max(max_degree, current_degree);
  }
  return {min_degree, max_degree};
}

bool ListGraph::contains(int u, int v) const noexcept {
  if (!contains(u) || !contains(v)) {
    return false;
  }
  for (const int &w : _adjList.at(u)) {
    if (w == v)
      return true;
  }
  return false;
}

bool ListGraph::contains(int v) const noexcept {
  return _adjList.count(v) != 0;
}

void ListGraph::add_edge(int u, int v) {
  if (!contains(u) || !contains(v)) {
    throw std::out_of_range("Endpoints of edge are not in the graph");
  }

  if (u == v) {
    throw std::invalid_argument("Self-loops are not allowed in simple graphs");
  }

  if (!contains(u, v)) { // não adiciona arestas repetidas
    _adjList[u].push_back(v);
    _adjList[v].push_back(u);
  }
}

void ListGraph::remove_edge(int u, int v) {
  if (contains(u, v)) {
    for (auto it = _adjList[u].begin(); it != _adjList[u].end(); ++it) {
      if (*it == v) {
        _adjList[u].erase(it);
        break;
      }
    }
    for (auto it = _adjList[v].begin(); it != _adjList[v].end(); ++it) {
      if (*it == u) {
        _adjList[v].erase(it);
        break;
      }
    }
  }
}

void ListGraph::remove_vertex(int v) {
  if (contains(v)) {
    for (const int &w : _adjList[v]) {
      for (auto it = _adjList[w].begin(); it != _adjList[w].end(); ++it) {
        if (*it == v) {
          _adjList[w].erase(it);
          break;
        }
      }
    }
    _adjList.erase(v);
  }
}

void ListGraph::for_each_vertex(const VertexCallback &func) const {
  for (const auto &par : _adjList) {
    func(par.first);
  }
}

void ListGraph::for_each_edge(const EdgeCallback &func) const {
  for (const auto &par : _adjList) {
    for (const auto &neighbor : par.second) {
      func(par.first, neighbor);
    }
  }
}

void ListGraph::for_each_neighbor(int v, const VertexCallback &func) const {
  if (!contains(v)) {
    throw std::out_of_range("Vertex not found");
  }
  for (auto it = _adjList.at(v).begin(); it != _adjList.at(v).end(); ++it) {
    func(*it);
  }
}

std::unordered_set<int> ListGraph::get_vertices() const {
  std::unordered_set<int> vertices;
  for (const auto &p : _adjList) {
    vertices.insert(p.first);
  }
  return vertices;
}

void ListGraph::print() const {
  for (auto it = _adjList.begin(); it != _adjList.end(); ++it) {
    std::cout << it->first << ": ";
    for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
      std::cout << *it2 << " ";
    }
    std::cout << '\n';
  }
}
