#pragma once
#include "ListGraph.hpp"
#include "MatrixGraph.hpp"
#include <boost/dynamic_bitset/dynamic_bitset.hpp>

#define LABEL_ZERO 0
#define LABEL_ONE 1
#define LABEL_TWO 2

class Chromosome {
private:
  int _fitness;
  int _size;
  boost::dynamic_bitset<> _genes0;
  boost::dynamic_bitset<> _genes1;
  boost::dynamic_bitset<> _genes2;

public:
  /**
   * @brief Construct a new Chromosome object
   */
  Chromosome() : _fitness(std::numeric_limits<int>::max()), _size(0) {}

  /**
   * @brief Construct a new Chromosome object
   */
  Chromosome(unsigned size);

  /**
   * @brief Construct a new Chromosome object
   */
  Chromosome(unsigned size, int value);

  /**
   * @brief Return the size of the chromosome
   */
  size_t size() const;

  /**
   * @brief Set the value object
   */
  void set_value(int index, int value);

  /**
   * @brief Get the value object
   */
  int get_value(int index) const;

  /**
   * @brief Get the fitness of the chromosome
   */
  int get_fitness() const;

  /**
   * @brief Calculate the fitness value and saves it in the attribute '_fitness'
   */
  int calculate_fitness();

  /**
   * @brief Overload of operator<< for being used with std::cout object
   */
  friend std::ostream &operator<<(std::ostream &out, const Chromosome &chr);

  /**
   * @brief Transform the chromosome into a valid solution for Total Roman
   Domination
   * @param graph The graph to fix the solution for
   */
  void fix_l(const ListGraph &graph);
  void fix_m(const MatrixGraph &graph);
};
