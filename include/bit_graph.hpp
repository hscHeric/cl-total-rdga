#pragma once

#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <unordered_map>
#include <unordered_set>
#include <vector>

/**
 * @brief Abstract class that defines the interface for a bit-based graph.
 *
 * This class provides a base structure for representing a graph,
 * allowing different storage implementations such as adjacency lists
 * or adjacency matrices using bitsets for memory efficiency.
 */
class BitGraph {
protected:
  /// Number of vertices in the graph
  int vertex_count;

public:
  /**
   * @brief Constructs a BitGraph with a given number of vertices.
   * @param n Number of vertices in the graph.
   */
  explicit BitGraph(int n) : vertex_count(n) {}

  /// Virtual destructor to allow proper cleanup in derived classes.
  virtual ~BitGraph() = default;

  /**
   * @brief Returns the number of vertices (order) of the graph.
   * @return Number of vertices in the graph.
   */
  int order() const { return vertex_count; }

  /**
   * @brief Returns the number of edges in the graph.
   * @return Total number of edges in the graph.
   */
  virtual int size() const = 0;

  /**
   * @brief Returns the degree of a given vertex.
   * @param v The vertex whose degree is to be determined.
   * @return The number of edges incident to vertex v.
   */
  virtual int degree(int v) const = 0;

  /**
   * @brief Returns the maximum degree among all vertices.
   * @return The maximum vertex degree in the graph.
   */
  virtual int max_degree() const = 0;

  /**
   * @brief Returns the minimum degree among all vertices.
   * @return The minimum vertex degree in the graph.
   */
  virtual int min_degree() const = 0;

  /**
   * @brief Checks if an edge (u, v) exists in the graph.
   * @param u First vertex.
   * @param v Second vertex.
   * @return True if the edge (u, v) exists, otherwise false.
   */
  virtual bool contains(int u, int v) const = 0;

  /**
   * @brief Checks if a vertex exists in the graph.
   * @param v Vertex to check.
   * @return True if the vertex exists, otherwise false.
   */
  virtual bool contains(int v) const = 0;

  /**
   * @brief Adds an edge between two vertices.
   * @param u First vertex.
   * @param v Second vertex.
   */
  virtual void add_edge(int u, int v) = 0;

  /**
   * @brief Removes an edge between two vertices.
   * @param u First vertex.
   * @param v Second vertex.
   */
  virtual void remove_edge(int u, int v) = 0;

  /**
   * @brief Removes a vertex and all its incident edges.
   * @param v The vertex to be removed.
   */
  virtual void remove_vertex(int v) = 0;

  /**
   * @brief Returns a set of all vertices in the graph.
   * @return A set containing all vertices.
   */
  virtual std::unordered_set<int> get_vertices() const = 0;

  /**
   * @brief Prints the graph structure.
   */
  virtual void print() const = 0;
};

/**
 * @brief Implementation of BitGraph using an adjacency list with bitsets.
 *
 * This class represents a graph where each vertex maintains a bitset to
 * indicate its adjacent vertices, optimizing memory usage compared to
 * standard adjacency lists.
 */
class BitListGraph : public BitGraph {
private:
  /// Adjacency list representation using bitsets
  std::unordered_map<int, boost::dynamic_bitset<>> adj_list;

public:
  /**
   * @brief Constructs a BitListGraph with a given number of vertices.
   * @param n Number of vertices.
   */
  explicit BitListGraph(int n);

  int size() const override;
  int degree(int v) const override;
  int max_degree() const override;
  int min_degree() const override;

  bool contains(int u, int v) const override;
  bool contains(int v) const override;

  void add_edge(int u, int v) override;
  void remove_edge(int u, int v) override;
  void remove_vertex(int v) override;

  std::unordered_set<int> get_vertices() const override;
  void print() const override;
};

/**
 * @brief Implementation of BitGraph using an adjacency matrix with bitsets.
 *
 * This class represents a graph where edges are stored in a bitset matrix,
 * optimizing space compared to a traditional boolean matrix.
 */
class BitMatrixGraph : public BitGraph {
private:
  /// Adjacency matrix representation using bitsets
  std::vector<boost::dynamic_bitset<>> adj_matrix;

public:
  /**
   * @brief Constructs a BitMatrixGraph with a given number of vertices.
   * @param n Number of vertices.
   */
  explicit BitMatrixGraph(int n);

  int size() const override;
  int degree(int v) const override;
  int max_degree() const override;
  int min_degree() const override;

  bool contains(int u, int v) const override;
  bool contains(int v) const override;

  void add_edge(int u, int v) override;
  void remove_edge(int u, int v) override;
  void remove_vertex(int v) override;

  std::unordered_set<int> get_vertices() const override;
  void print() const override;
};
