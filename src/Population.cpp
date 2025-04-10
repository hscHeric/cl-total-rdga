#include "Population.hpp"

// Population.cpp
Population::Population( size_t size, const std::vector<Chromosome> & initialChromosomes ) : chromosomes( initialChromosomes ), size( size ) {
  if ( chromosomes.size() != size ) {
    std::string str = "Initial chromosome count of " + std::to_string( chromosomes.size() ) + " must match population size of " + std::to_string( size );
    throw std::runtime_error( str );
  }
}

size_t Population::getSize() const {
  return size;
}

const std::vector<Chromosome> & Population::getChromosomes() const {
  return chromosomes;
}

Chromosome Population::tournamentSelection( size_t tournamentSize ) {
  std::random_device              rd;
  std::mt19937                    gen( rd() );
  std::uniform_int_distribution<> dist( 0, size - 1 );

  Chromosome best = chromosomes[dist( gen )];
  for ( size_t i = 1; i < tournamentSize; ++i ) {
    Chromosome candidate = chromosomes[dist( gen )];
    if ( candidate.get_fitness() < best.get_fitness() ) {
      best = candidate;
    }
  }
  return best;
}

void Population::evolve_l( const ListGraph & graph, size_t tournamentSize, double crossoverRate, double elitismRate, double mutationRate ) {
  std::vector<Chromosome> newChromosomes;
  newChromosomes.reserve( size );

  // apply elitism -------------------------------------------
  // order the chromosomes in the population in
  // ascending order of their fitness values
  int size_elitism = size * elitismRate;
  if ( size_elitism > 0 ) {
    std::vector<Chromosome> temp_pop = chromosomes;

    std::sort( temp_pop.begin(), temp_pop.end(), []( const Chromosome & a, const Chromosome & b ) { return a.get_fitness() < b.get_fitness(); } );

    // leaves only the k first chromosomes with best fitness values,
    // where k = _size_population * _elitism_rate.
    for ( int i = 0; i < size_elitism; ++i ) {
      newChromosomes.push_back( temp_pop[i] );
    }
  }
  // end elitism ---------------------------------------------

  std::random_device               rd;
  std::mt19937                     gen( rd() );
  std::uniform_real_distribution<> dist( 0.0, 1.0 );

  while ( newChromosomes.size() < size ) {
    Chromosome parent1 = tournamentSelection( tournamentSize );
    Chromosome parent2 = tournamentSelection( tournamentSize );

    if ( dist( gen ) < crossoverRate ) {
      auto [child1, child2] = onePointCrossover( parent1, parent2 );

      // Aplicando a correção dos cromossomos antes de adicioná-los à população
      child1.fix_l( graph );
      child2.fix_l( graph );

      mutation_l( child1, mutationRate, graph );
      mutation_l( child2, mutationRate, graph );

      newChromosomes.push_back( child1 );
      newChromosomes.push_back( child2 );
    } else {
      mutation_l( parent1, mutationRate, graph );
      mutation_l( parent2, mutationRate, graph );

      newChromosomes.push_back( parent1 );
      newChromosomes.push_back( parent2 );
    }
  }

  if ( newChromosomes.size() > size ) {
    newChromosomes.pop_back();
  }

  chromosomes.swap( newChromosomes );  // troca o conteudo dos dois vectors
}

void Population::evolve_m( const MatrixGraph & graph, size_t tournamentSize, double crossoverRate, double elitismRate, double mutationRate ) {
  std::vector<Chromosome> newChromosomes;
  newChromosomes.reserve( size );

  // apply elitism -------------------------------------------
  // order the chromosomes in the population in
  // ascending order of their fitness values
  int size_elitism = size * elitismRate;
  if ( size_elitism > 0 ) {
    std::vector<Chromosome> temp_pop = chromosomes;

    std::sort( temp_pop.begin(), temp_pop.end(), []( const Chromosome & a, const Chromosome & b ) { return a.get_fitness() < b.get_fitness(); } );

    // leaves only the k first chromosomes with best fitness values,
    // where k = _size_population * _elitism_rate.
    for ( int i = 0; i < size_elitism; ++i ) {
      newChromosomes.push_back( temp_pop[i] );
    }
  }
  // end elitism ---------------------------------------------

  std::random_device               rd;
  std::mt19937                     gen( rd() );
  std::uniform_real_distribution<> dist( 0.0, 1.0 );

  while ( newChromosomes.size() < size ) {
    Chromosome parent1 = tournamentSelection( tournamentSize );
    Chromosome parent2 = tournamentSelection( tournamentSize );

    if ( dist( gen ) < crossoverRate ) {
      auto [child1, child2] = onePointCrossover( parent1, parent2 );

      // Aplicando a correção dos cromossomos antes de adicioná-los à população
      child1.fix_m( graph );
      child2.fix_m( graph );

      mutation_m( child1, mutationRate, graph );
      mutation_m( child2, mutationRate, graph );

      newChromosomes.push_back( child1 );
      newChromosomes.push_back( child2 );
    } else {
      mutation_m( parent1, mutationRate, graph );
      mutation_m( parent2, mutationRate, graph );

      newChromosomes.push_back( parent1 );
      newChromosomes.push_back( parent2 );
    }
  }

  if ( newChromosomes.size() > size ) {
    newChromosomes.pop_back();
  }

  chromosomes.swap( newChromosomes );  // troca o conteudo dos dois vectors
}

std::pair<Chromosome, Chromosome> Population::onePointCrossover( const Chromosome & p1, const Chromosome & p2 ) {
  std::random_device              rd;
  std::mt19937                    gen( rd() );
  std::uniform_int_distribution<> dist( 1, p1.size() - 1 );

  size_t     crossoverPoint = dist( gen );
  Chromosome child1         = p1;
  Chromosome child2         = p2;

  for ( size_t i = crossoverPoint; i < p1.size(); ++i ) {
    child1.set_value( i, p2.get_value( i ) );
    child2.set_value( i, p1.get_value( i ) );
  }

  child1.calculate_fitness();
  child2.calculate_fitness();

  return { child1, child2 };
}

void Population::mutation_l( Chromosome & chr, float mutationRate, const ListGraph & graph ) {
  // Criar um motor de números aleatórios
  std::random_device rd;
  std::mt19937       gen( rd() );  // Mersenne Twister 19937

  // Criar uma distribuição uniforme entre 0 e 1
  std::uniform_real_distribution<double> dist( 0.0, 1.0 );

  for ( int i = 0; i < (int)chr.size(); ++i ) {
    // Gerar um número aleatório
    double rn = dist( gen );

    // mutação no gene
    if ( rn < mutationRate ) {
      chr.set_value( i, 0 );
    }
  }

  // conserta o cromossomo
  chr.fix_l( graph );
}

void Population::mutation_m( Chromosome & chr, float mutationRate, const MatrixGraph & graph ) {
  // Criar um motor de números aleatórios
  std::random_device rd;
  std::mt19937       gen( rd() );  // Mersenne Twister 19937

  // Criar uma distribuição uniforme entre 0 e 1
  std::uniform_real_distribution<double> dist( 0.0, 1.0 );

  for ( int i = 0; i < (int)chr.size(); ++i ) {
    // Gerar um número aleatório
    double rn = dist( gen );

    // mutação no gene
    if ( rn < mutationRate ) {
      chr.set_value( i, 0 );
    }
  }

  // conserta o cromossomo
  chr.fix_m( graph );
}
