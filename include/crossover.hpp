#pragma once

#include "bit_graph.hpp"
#include "chromosome.hpp"
class Crossover {
public:
  virtual ~Crossover() = default;
  virtual std::pair<Chromosome, Chromosome> crossover(const Chromosome &parent1,
                                                      const Chromosome &parent2,
                                                      Graph &graph) const = 0;
};

class SinglePoint : public Crossover {
public:
  explicit SinglePoint(double crossover_rate);
  std::pair<Chromosome, Chromosome> crossover(const Chromosome &parent1,
                                              const Chromosome &parent2,
                                              Graph &graph) const override;

private:
  double _crossover_rate;
};
