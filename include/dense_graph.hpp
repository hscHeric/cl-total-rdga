#pragma once

#include <boost/dynamic_bitset.hpp>
#include <iostream>
#include <list>
#include <string>
#include <unordered_map>
#include <unordered_set>

class DenseGraph {
private:
  // the adjacency list is implemented as a hash table whose keys
  // are the vertices and the associated values are the list of neighbors
  std::unordered_map<short, boost::dynamic_bitset<>> _adjList;

  // name of the graph
  std::string name;

public:
  /**
   * @brief Construct a new Graph object
   */
  DenseGraph() = default;

  /**
   * @brief Construct a new Graph object with n vertices numbered from 0 to n-1
   */
  explicit DenseGraph(int n);

  /**
   * @brief Get the name of the graph
   */
  std::string get_name() const;

  /**
   * @brief Set the name of the graph
   */
  void set_name(const std::string &nm);

  /**
   * @brief Return the number of vertices in the graph
   */
  int getVertexCount() const;

  /**
   * @brief Return the number of edges in the graph
   */
  int getEdgeCount() const;

  /**
   * @brief Returns the degree of a vertex v
   */
  int degree(int v) const;

  /**
   * @brief Return the maximum degree of the graph
   */
  int max_degree() const;

  /**
   * @brief Return the minimum degree of the graph
   */
  unsigned min_degree() const;

  /**
   * @brief Check if edge (u,v) is in the graph.
   * Attention: before invoquing this function, guarantee that u and v
   * are vertices of the graph; if they do not belong to the graph,
   * then an exception may be thrown.
   */
  bool contains(short u, short v) const;

  /**
   * @brief Check if vertex v is in the graph
   */
  bool contains(int v) const;

  /**
   * @brief Add edge (u,v) in the graph.
   * Atention: before adding the edge (u,v) guarantee that vertices u and v
   * belong to the graph; if they do not belong to the graph,
   * then an exception may be thrown.
   */
  void addEdge(int u, int v);

  /**
   * @brief Delete edge (u,v) from the graph
   */
  void delEdge(int u, int v);

  /**
   * @brief Delete vertex v and its incident edges
   */
  void delVertex(int v);

  /**
   * @brief Return the vertices as an std::unordered_set
   */
  std::unordered_set<int> get_vertices() const;

  /**
   * @brief Return a reference to the bitset array of each vertex v
   */
  boost::dynamic_bitset<> &get_neighbors_bitset(int v);

  /**
   * @brief Return a list os pairs (vertex, degree)
   */
  std::list<std::pair<int, int>> get_all_pairs_vertex_degree();

  /**
   * @brief Overload of operator<<
   */
  friend std::ostream &operator<<(std::ostream &out, const DenseGraph &g);

#if 0
    void dfs_visit(int u, std::vector<bool>& discovered, int& number_of_vertices, int& min_degree) {
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

    std::vector<std::pair<int,int>> connected_components() {
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
};
