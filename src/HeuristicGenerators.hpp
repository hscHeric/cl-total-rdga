#pragma once
#include "Chromosome.hpp"
#include "ListGraph.hpp"
#include "MatrixGraph.hpp"

class HeuristicGenerators {
public:
  static Chromosome h1_l( const ListGraph & graph );
  static Chromosome h2_l( const ListGraph & graph );
  static Chromosome h3_l( const ListGraph & graph );
  static Chromosome h4_l( const ListGraph & graph );
  static Chromosome h5_l( const ListGraph & graph );

  static Chromosome h1_m( const MatrixGraph & graph );
  static Chromosome h2_m( const MatrixGraph & graph );
  static Chromosome h3_m( const MatrixGraph & graph );
  static Chromosome h4_m( const MatrixGraph & graph );
  static Chromosome h5_m( const MatrixGraph & graph );
};
