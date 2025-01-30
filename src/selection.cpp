#include "../include/selection.hpp"
#include "../include/population.hpp"

#include <random>
KTournamentSelection::KTournamentSelection(unsigned k) : _k(k) {}

const Chromosome &KTournamentSelection::select(Population &population) const {
  if (population.size() == 0) {
    throw std::runtime_error("Population is empty.");
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<size_t> dist(0, population.size() - 1);

  std::vector<const Chromosome *> tournament;
  for (unsigned i = 0; i < _k; ++i) {
    size_t index = dist(gen);
    tournament.push_back(&population.get_chromosomes()[index]);
  }

  return **std::max_element(tournament.begin(), tournament.end(),
                            [](const Chromosome *a, const Chromosome *b) {
                              return a->get_fitness() < b->get_fitness();
                            });
}
