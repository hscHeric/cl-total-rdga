#pragma once

#include "Chromosome.hpp"
#include <cstddef>
#include <vector>
class Population {
private:
  std::vector<Chromosome> chromosomes;
  size_t size;

public:
  Population(size_t size, const std::vector<Chromosome> &initialChromosomes);
  size_t getSize() const;
  const std::vector<Chromosome> &getChromosomes() const;
  void evolve_l(const ListGraph &graph, size_t tournamentSize,
                double crossoverRate, double elitismRate, double mutationRate);
  void evolve_m(const MatrixGraph &graph, size_t tournamentSize,
                double crossoverRate, double elitismRate, double mutationRate);

private:
  Chromosome tournamentSelection(size_t tournamentSize);
  std::pair<Chromosome, Chromosome> onePointCrossover(const Chromosome &p1,
                                                      const Chromosome &p2);
  void mutation_l(Chromosome &chr, float mutationRate, const ListGraph &graph);
  void mutation_m(Chromosome &chr, float mutationRate,
                  const MatrixGraph &graph);
};
