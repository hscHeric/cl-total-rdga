#pragma once
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <functional>
#include <iostream>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using BitSet = boost::dynamic_bitset<>;
using VertexCallback = std::function<void(int)>;
using EdgeCallback = std::function<void(int, int)>;

struct PairHash {
  template <class T1, class T2>
  std::size_t operator()(const std::pair<T1, T2> &p) const {
    auto h1 = std::hash<T1>{}(p.first);
    auto h2 = std::hash<T2>{}(p.second);
    return h1 ^ (h2 << 1);
  }
};

/**
 * @brief Os grafos herdados desta classe possuem n vértices. Essa quantidade n
 * é informada no construtor da classe, quando o grafo é instanciado. Ademais, n
 * deve ser um inteiro maior que zero. Depois de instanciado o grafo, a
 * quantidade n de vértices não pode ser aumentada. Não é possível adicionar
 * vértices no grafo depois dele ser instanciado. Os vértices são numerados de 0
 * a n-1. Logo, as extremidades das arestas também deve estar dentro desse
 * intervalo de 0 a n-1. É possível remover vértices do grafo: a remoção de
 * vértices não altera a numeração inicial dada a cada vértice, ou seja, após a
 * remoção de um vértice v, os vértices que ficam continuam com a mesma
 * numeração que tinham desde a instanciação do grafo. Também é possível remover
 * e adicionar arestas no grafo depois dele instanciado.
 */
class Graph {
public:
  explicit Graph(int n) {
    if (n <= 0) {
      throw std::runtime_error("invalid graph order");
    }
  }
  virtual ~Graph() = default;
  // Graph(Graph &&) = default;

  // Ainda não entendi porquê você deletou o construtor de cópia.
  // Ainda vou estudar as implicações de remover ou de deixar ele.
  // Graph(const Graph &) = delete;
  // Graph& operator=(const Graph &) = default;

  [[nodiscard]] virtual int order() const noexcept = 0;
  [[nodiscard]] virtual int size() const noexcept = 0;
  [[nodiscard]] virtual int degree(int v) const = 0;
  [[nodiscard]] virtual std::pair<int, int> degree_range() const noexcept = 0;
  [[nodiscard]] virtual bool contains(int u, int v) const noexcept = 0;
  [[nodiscard]] virtual bool contains(int v) const noexcept = 0;

  virtual void add_edge(int u, int v) = 0;
  virtual void remove_edge(int u, int v) = 0;
  virtual void remove_vertex(int v) = 0;

  virtual void for_each_vertex(const VertexCallback &func) const = 0;
  virtual void for_each_edge(const EdgeCallback &func) const = 0;
  virtual void for_each_neighbor(int v, const VertexCallback &func) const = 0;
  [[nodiscard]] virtual std::unordered_set<int> get_vertices() const = 0;
  virtual void print() const = 0;
  [[nodiscard]] virtual int choose_rng() const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::unordered_set<int> vertices = get_vertices();

    if (vertices.empty()) {
      throw std::runtime_error("Graph is empty, cannot choose a vertex.");
    }
    std::uniform_int_distribution<> distrib(0, vertices.size() - 1);
    auto it = vertices.begin();
    std::advance(it, distrib(gen));
    return *it;
  };
};
