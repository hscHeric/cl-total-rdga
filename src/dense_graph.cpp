#include "../include/dense_graph.hpp"
#include <boost/dynamic_bitset.hpp>
#include <iostream>
#include <limits>
#include <list>
#include <string>
#include <unordered_map>
#include <unordered_set>

DenseGraph::DenseGraph(int n) {
  for (int i = 0; i < n; ++i) {
    _adjList[i];
    _adjList[i].resize(n);
  }
}

std::string DenseGraph::get_name() const { return name; }

void DenseGraph::set_name(const std::string &nm) { name = nm; }

int DenseGraph::getVertexCount() const {
  return static_cast<int>(_adjList.size());
}

int DenseGraph::getEdgeCount() const {
  int edgeCount = 0;
  for (const auto &p : _adjList) {
    edgeCount += p.second.count();
  }
  return edgeCount / 2;
}

int DenseGraph::degree(int v) const { return _adjList.at(v).count(); }

int DenseGraph::max_degree() const {
  int maxDegree = 0;
  for (const auto &p : _adjList) {
    int current_degree = degree(p.first);
    if (current_degree > maxDegree) {
      maxDegree = current_degree;
    }
  }
  return maxDegree;
}

unsigned DenseGraph::min_degree() const {
  int minDegree = std::numeric_limits<int>::max();
  for (const auto &p : _adjList) {
    int currentDegree = degree(p.first);
    if (currentDegree < minDegree) {
      minDegree = currentDegree;
    }
  }
  return minDegree;
}

bool DenseGraph::contains(short u, short v) const {
  return (_adjList.at(u)[v] == 1) ? true : false;
}

bool DenseGraph::contains(int v) const { return _adjList.count(v) != 0; }

void DenseGraph::addEdge(int u, int v) {
  if (contains(u, v) || u == v) {
    return;
  }
  _adjList[u][v] = 1;
  _adjList[v][u] = 1;
}

void DenseGraph::delEdge(int u, int v) {
  _adjList[u][v] = 0;
  _adjList[v][u] = 0;
}

void DenseGraph::delVertex(int v) {
  if (_adjList.count(v) > 0) {
    for (unsigned i = 0; i < _adjList[v].size(); ++i) {
      if (_adjList[v][i] == 1) {
        _adjList.at(i)[v] = 0;
      }
    }
    _adjList.erase(v);
  }
}

std::unordered_set<int> DenseGraph::get_vertices() const {
  std::unordered_set<int> V;
  for (const auto &p : _adjList) {
    V.insert(p.first);
  }
  return V;
}

boost::dynamic_bitset<> &DenseGraph::get_neighbors_bitset(int v) {
  return _adjList.at(v);
}

std::list<std::pair<int, int>> DenseGraph::get_all_pairs_vertex_degree() {
  std::list<std::pair<int, int>> list_pairs;

  for (const auto &p : _adjList) {
    list_pairs.push_back({p.first, degree(p.first)});
  }

  return list_pairs;
}

std::ostream &operator<<(std::ostream &out, const DenseGraph &g) {
  for (auto it = g._adjList.cbegin(); it != g._adjList.cend(); ++it) {
    out << it->first << ": ";
    out << it->second;
    out << '\n';
  }
  return out;
}

#if 0
void DenseGraph::dfs_visit(int u, std::vector<bool>& discovered, int& number_of_vertices, int& min_degree) {
    discovered[u] = true;
    if(degree(u) < (unsigned)min_degree) {
        min_degree = degree(u);
    }
    for(const int& v : get_neighbors(u)) {
        if(!discovered[v]) {
            number_of_vertices++;
            if(degree(v) < (unsigned)min_degree) {
                min_degree = degree(v);
            }
            dfs_visit(v, discovered, number_of_vertices, min_degree);
        }
    }
}

std::vector<std::pair<int,int>> DenseGraph::connected_components() {
    std::vector<bool> discovered(order(), false);
    std::vector<std::pair<int,int>> components;
    for(const int& vertex : get_vertices()) {
        if(!discovered[vertex]) {
            int number_of_vertices = 1;
            int min_degree = std::numeric_limits<int>::max();
            dfs_visit(vertex, discovered, number_of_vertices, min_degree);
            components.push_back({number_of_vertices, min_degree});
        }
    }
    return components;
}
#endif
