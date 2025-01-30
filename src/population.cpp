#include "../include/population.hpp"
#include <algorithm>
#include <stdexcept>

Population::Population(unsigned size,
                       const std::vector<Chromosome> &initial_population)
    : _size(size), _chromosomes(initial_population) {
  if (initial_population.size() != size) {
    throw std::runtime_error(
        "Initial population size does not match the specified size.");
  }
}

void Population::evolve(const Selection &selector, const Crossover &crossover) {
  std::vector<Chromosome> new_population;
  while (new_population.size() < _size) {

    const Chromosome &parent1 = selector.select(*this);
    const Chromosome &parent2 = selector.select(*this);

    // Aplicação do crossover
    auto [child1, child2] = crossover.crossover(parent1, parent2);

    // Adição dos filhos à nova população
    new_population.push_back(child1);
    if (new_population.size() < _size) {
      new_population.push_back(child2);
    }
  }
  _chromosomes = std::move(new_population);
}

Chromosome Population::get_best_chromosome() const {
  return *std::max_element(_chromosomes.begin(), _chromosomes.end(),
                           [](const Chromosome &a, const Chromosome &b) {
                             return a.get_fitness() < b.get_fitness();
                           });
}

size_t Population::size() const { return _chromosomes.size(); }

const std::vector<Chromosome> &Population::get_chromosomes() const {
  return _chromosomes;
}
