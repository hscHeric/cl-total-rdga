#pragma once

#include "../include/chromosome.hpp"
#include "../include/crossover.hpp"
#include "../include/selection.hpp"
#include "dense_graph.hpp"
#include <vector>

class Population {
public:
  Population(unsigned size, const std::vector<Chromosome> &initial_population);
  void evolve(const Selection &selector, const Crossover &crossover,
              DenseGraph &graph);
  Chromosome get_best_chromosome() const;
  size_t size() const;
  const std::vector<Chromosome> &get_chromosomes() const;

private:
  std::vector<Chromosome> _chromosomes;
  unsigned _size;
};
