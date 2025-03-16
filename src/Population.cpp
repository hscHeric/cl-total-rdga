
#include "Population.hpp"

// Population.cpp
Population::Population(size_t size, const std::vector<Chromosome> &initialChromosomes)
    : chromosomes(initialChromosomes), size(size) {
  if (chromosomes.size() != size) {
    throw std::runtime_error(
        "Initial chromosome count must match population size");
  }
}

size_t Population::getSize() const { return size; }

const std::vector<Chromosome> &Population::getChromosomes() const {
  return chromosomes;
}

Chromosome Population::tournamentSelection(size_t tournamentSize) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dist(0, size - 1);

  Chromosome best = chromosomes[dist(gen)];
  for (size_t i = 1; i < tournamentSize; ++i) {
    Chromosome candidate = chromosomes[dist(gen)];
    if (candidate.get_fitness() < best.get_fitness()) {
      best = candidate;
    }
  }
  return best;
}

void Population::evolve_l(const ListGraph &graph, size_t tournamentSize, double crossoverRate) {

  std::vector<Chromosome> newChromosomes;
  newChromosomes.reserve(size);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(0.0, 1.0);

  while (newChromosomes.size() < size) {
    Chromosome parent1 = tournamentSelection(tournamentSize);
    Chromosome parent2 = tournamentSelection(tournamentSize);

    if (dist(gen) < crossoverRate) {
      auto [child1, child2] = onePointCrossover(parent1, parent2);

      // Aplicando a correção dos cromossomos antes de adicioná-los à população
      child1.fix_l(graph);
      child2.fix_l(graph);

      mutation_l(child1, 0.05, graph);
      mutation_l(child2, 0.05, graph);

      newChromosomes.push_back(child1);
      newChromosomes.push_back(child2);
    } else {
      newChromosomes.push_back(parent1);
      newChromosomes.push_back(parent2);
    }
  }

  if (newChromosomes.size() > size) {
    newChromosomes.pop_back();
  }

  chromosomes = newChromosomes;
}

void Population::evolve_m(const MatrixGraph &graph, size_t tournamentSize, double crossoverRate) {

  std::vector<Chromosome> newChromosomes;
  newChromosomes.reserve(size);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(0.0, 1.0);

  while (newChromosomes.size() < size) {
    Chromosome parent1 = tournamentSelection(tournamentSize);
    Chromosome parent2 = tournamentSelection(tournamentSize);

    if (dist(gen) < crossoverRate) {
      auto [child1, child2] = onePointCrossover(parent1, parent2);

      // Aplicando a correção dos cromossomos antes de adicioná-los à população
      child1.fix_m(graph);
      child2.fix_m(graph);

      newChromosomes.push_back(child1);
      newChromosomes.push_back(child2);
    } else {
      newChromosomes.push_back(parent1);
      newChromosomes.push_back(parent2);
    }

    
  }

  if (newChromosomes.size() > size) {
    newChromosomes.pop_back();
  }

  chromosomes.swap(newChromosomes); // troca o conteudo dos dois vectors
}

std::pair<Chromosome, Chromosome> 
Population::onePointCrossover(const Chromosome &p1, const Chromosome &p2) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dist(1, p1.size() - 1);

  size_t crossoverPoint = dist(gen);
  Chromosome child1 = p1;
  Chromosome child2 = p2;

  for (size_t i = crossoverPoint; i < p1.size(); ++i) {
    child1.set_value(i, p2.get_value(i));
    child2.set_value(i, p1.get_value(i));
  }

  child1.calculate_fitness();
  child2.calculate_fitness();

  return {child1, child2};
}


void Population::mutation_l(Chromosome &chr, float mutationRate, const ListGraph& graph) {
  // se o grafo for esparso, realiza mutacao
  if(graph.get_density() <= 0.5) {
    // Criar um motor de números aleatórios
    std::random_device rd;
    std::mt19937 gen(rd()); // Mersenne Twister 19937

    // Criar uma distribuição uniforme entre 0 e 1
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for(int i = 0; i < (int)chr.size(); ++i) {
      // Gerar um número aleatório
      double rn = dist(gen);

      if(rn < mutationRate) {
        chr.set_value(i, 0);
      }
    }

    chr.fix_l(graph);
  }
}
