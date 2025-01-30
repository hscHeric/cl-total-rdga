#pragma once

#include "../include/chromosome.hpp"

class Population;

class Selection {
public:
  virtual ~Selection() = default;
  virtual const Chromosome &select(Population &population) const = 0;
};

class KTournamentSelection : public Selection {
public:
  explicit KTournamentSelection(unsigned k);
  const Chromosome &select(Population &population) const override;

private:
  unsigned _k;
};
