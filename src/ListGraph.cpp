#include "ListGraph.hpp"

#include <algorithm>
#include <stdexcept>

ListGraph::ListGraph( int n ) {
  if ( n <= 0 ) {
    throw std::runtime_error( "invalid graph order" );
  }
  // Inicializa o grafo com n vértices numerados de 0 a n-1
  for ( int i = 0; i < n; i++ ) {
    _adjList[i] = std::vector<int>();
  }
}

int ListGraph::order() const noexcept {
  return _adjList.size();
}

int ListGraph::size() const noexcept {
  int count = 0;
  for ( const auto & [vertex, neighbors] : _adjList ) {
    count += neighbors.size();
  }
  // Cada aresta é contada duas vezes (uma para cada extremidade)
  return count / 2;
}

int ListGraph::degree( int v ) const {
  if ( !contains( v ) ) {
    throw std::runtime_error( "Vertex does not exist (function degree)" );
  }
  return _adjList.at( v ).size();
}

std::pair<int, int> ListGraph::degree_range() const noexcept {
  if ( _adjList.empty() ) {
    return { 0, 0 };
  }

  int min_degree = std::numeric_limits<int>::max();
  int max_degree = 0;

  for ( const auto & [vertex, neighbors] : _adjList ) {
    int deg    = neighbors.size();
    min_degree = std::min( min_degree, deg );
    max_degree = std::max( max_degree, deg );
  }

  return { min_degree, max_degree };
}

bool ListGraph::contains( int u, int v ) const noexcept {
  if ( !contains( u ) || !contains( v ) ) {
    return false;
  }

  // Verifica se v está na lista de adjacência de u
  const auto & neighbors = _adjList.at( u );
  return std::find( neighbors.begin(), neighbors.end(), v ) != neighbors.end();
}

bool ListGraph::contains( int v ) const noexcept {
  return _adjList.find( v ) != _adjList.end();
}

void ListGraph::add_edge( int u, int v ) {
  // Ignora self-loops
  if ( u == v ) {
    return;
  }

  // Verifica se os vértices existem
  if ( !contains( u ) || !contains( v ) ) {
    throw std::runtime_error( "Vertex does not exist (function add_edge)" );
  }

  // Ignora arestas múltiplas (verifica se a aresta já existe)
  if ( contains( u, v ) ) {
    return;
  }

  // Adiciona a aresta (nos dois sentidos, já que é um grafo não direcionado)
  _adjList[u].push_back( v );
  _adjList[v].push_back( u );
}

void ListGraph::remove_edge( int u, int v ) {
  if ( !contains( u ) || !contains( v ) ) {
    throw std::runtime_error( "Vertex does not exist (function remove_edge)" );
  }

  // Remove v da lista de adjacência de u
  auto & u_neighbors = _adjList[u];
  u_neighbors.erase( std::remove( u_neighbors.begin(), u_neighbors.end(), v ), u_neighbors.end() );

  // Remove u da lista de adjacência de v
  auto & v_neighbors = _adjList[v];
  v_neighbors.erase( std::remove( v_neighbors.begin(), v_neighbors.end(), u ), v_neighbors.end() );
}

void ListGraph::remove_vertex( int v ) {
  if ( !contains( v ) ) {
    throw std::runtime_error( "Vertex does not exist (function remove_vertex)" );
  }

  // Remove todas as arestas incidentes em v
  for ( auto & [vertex, neighbors] : _adjList ) {
    if ( vertex != v ) {
      neighbors.erase( std::remove( neighbors.begin(), neighbors.end(), v ), neighbors.end() );
    }
  }

  // Remove o vértice v
  _adjList.erase( v );
}

/*
void ListGraph::for_each_vertex(const VertexCallback &func) const {
  for (const auto &[vertex, _] : _adjList) {
    func(vertex);
  }
}

void ListGraph::for_each_edge(const EdgeCallback &func) const {
  std::unordered_set<std::pair<int, int>, PairHash> processed_edges;

  for (const auto &[u, neighbors] : _adjList) {
    for (int v : neighbors) {
      // Garante que cada aresta seja processada apenas uma vez
      auto edge1 = std::make_pair(std::min(u, v), std::max(u, v));
      if (processed_edges.find(edge1) == processed_edges.end()) {
        func(u, v);
        processed_edges.insert(edge1);
      }
    }
  }
}

void ListGraph::for_each_neighbor(int v, const VertexCallback &func) const {
  if (!contains(v)) {
    throw std::runtime_error("Vertex does not exist");
  }

  for (int neighbor : _adjList.at(v)) {
    func(neighbor);
  }
}
  */

const std::vector<int> & ListGraph::get_neighbors( int v ) const {
  if ( !contains( v ) ) {
    throw std::runtime_error( "vertex does not exists (function get_neighbors)" );
  }
  return _adjList.at( v );
}

std::unordered_set<int> ListGraph::get_vertices() const {
  std::unordered_set<int> vertices;
  for ( const auto & [vertex, _] : _adjList ) {
    vertices.insert( vertex );
  }
  return vertices;
}

std::unordered_set<int> ListGraph::get_isolated_vertices() const {
  std::unordered_set<int> vertices;
  for ( const auto & [vertex, _] : _adjList ) {
    if ( degree( vertex ) == 0 ) {
      vertices.insert( vertex );
    }
  }
  return vertices;
}

bool ListGraph::is_isolated_vertex( int vertex ) const {
  return contains( vertex ) && degree( vertex ) == 0;
}

float ListGraph::get_density() const {
  return static_cast<float>( size() * 2 ) / ( order() * ( order() - 1 ) );
}

void ListGraph::print() const {
  std::cout << "ListGraph: " << order() << " vertices, " << size() << " edges\n";
  for ( const auto & [vertex, neighbors] : _adjList ) {
    std::cout << vertex << ": ";
    for ( int neighbor : neighbors ) {
      std::cout << neighbor << " ";
    }
    std::cout << "\n";
  }
}
