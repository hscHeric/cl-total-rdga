#include "../include/crossover.hpp"
#include <random>

SinglePoint::SinglePoint(double crossover_rate)
    : _crossover_rate(crossover_rate) {}

std::pair<Chromosome, Chromosome>
SinglePoint::crossover(const Chromosome &parent1, const Chromosome &parent2,
                       DenseGraph &graph) const {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(0.0, 1.0);

  if (dist(gen) > _crossover_rate) {
    return {parent1, parent2};
  }

  int point = std::uniform_int_distribution<>(1, parent1.size() - 1)(gen);
  Chromosome child1 = parent1;
  Chromosome child2 = parent2;

  for (size_t i = point; i < parent1.size(); ++i) {
    int temp = child1.get_value(i);
    child1.set_value(i, child2.get_value(i));
    child2.set_value(i, temp);
  }

  child1.fix(graph);
  child2.fix(graph);

  return {child1, child2};
}
