#include "Chromosome.hpp"
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <iostream>

Chromosome::Chromosome(unsigned size, int value) {
  if (value < 0 || value > 2) {
    throw std::runtime_error("wrong value");
  }
  _size = size;
  _genes0.resize(size, 0);
  _genes1.resize(size, 0);
  _genes2.resize(size, 0);

  if (value == 0) {
    _genes0.set();
  } else if (value == 1) {
    _genes1.set();
  } else if (value == 2) {
    _genes2.set();
  }
  calculate_fitness();
}

Chromosome::Chromosome(unsigned size) : Chromosome(size, 0) {}

size_t Chromosome::size() const { return static_cast<size_t>(_size); }

void Chromosome::set_value(int index, int value) {
  if (static_cast<int>(_genes0.size()) < index) {
    std::cerr << "error: size of chromossome: " << _genes0.size() << std::endl;
    throw std::runtime_error("error: chromossome size < index 11");
  }
  if (value == 0) {
    _genes0[index] = 1;
    _genes1[index] = 0;
    _genes2[index] = 0;
  } else if (value == 1) {
    _genes0[index] = 0;
    _genes1[index] = 1;
    _genes2[index] = 0;
  } else if (value == 2) {
    _genes0[index] = 0;
    _genes1[index] = 0;
    _genes2[index] = 1;
  }
}

int Chromosome::get_value(int index) const {
  if (static_cast<int>(_genes0.size()) < index) {
    throw std::runtime_error("error: chromossome size < index 22");
  }
  return _genes1[index] * 1 + _genes2[index] * 2;
}

int Chromosome::get_fitness() const { return _fitness; }

int Chromosome::calculate_fitness() {
  _fitness = 0;
  for (int i = 0; i < _size; ++i) {
    _fitness += get_value(i);
  }
  return _fitness;
}

std::ostream &operator<<(std::ostream &out, const Chromosome &chr) {
  out << "[";
  for (int i = 0; i < chr._size; ++i) {
    out << chr.get_value(i) << " ";
  }
  out << "|| fitness: " << chr._fitness << "]";
  return out;
}

/*
void Chromosome::fix(const Graph &graph) {
  boost::dynamic_bitset<> already_dominated;
  already_dominated.resize(this->size(), 0);

  for (size_t u = 0; u < this->size(); ++u) {
    if (this->get_value(u) == LABEL_ZERO && !already_dominated[u]) {
      unsigned short num_active = 0;
      bool has_neighbor_with_label_two = false;
      int vertex_with_label_one = -1;
      int last_neighbor = -1;

      graph.for_each_neighbor(u, [&](int w) {
        last_neighbor = w;
        if (this->get_value(w) == LABEL_ONE) {
          num_active++;
          vertex_with_label_one = w;
        } else if (this->get_value(w) == LABEL_TWO) {
          num_active++;
          has_neighbor_with_label_two = true;
        }
      });

      if (num_active == 0 && last_neighbor != -1) {
        this->set_value(last_neighbor, LABEL_TWO);
      } else if (!has_neighbor_with_label_two && vertex_with_label_one != -1) {
        this->set_value(vertex_with_label_one, LABEL_TWO);
      }
      already_dominated[u] = 1;
    } else if (this->get_value(u) == LABEL_TWO) {
      unsigned short num_active = 0;
      int last_neighbor = -1;

      graph.for_each_neighbor(u, [&](int w) {
        last_neighbor = w;
        already_dominated[w] = 1;
        if (this->get_value(w) >= LABEL_ONE) {
          num_active++;
        }
      });

      if (num_active == 0 && last_neighbor != -1) {
        this->set_value(last_neighbor, LABEL_ONE);
      }
    } else if (this->get_value(u) == LABEL_ONE) {
      unsigned short num_active = 0;
      int last_neighbor = -1;

      graph.for_each_neighbor(u, [&](int w) {
        last_neighbor = w;
        if (this->get_value(w) >= LABEL_ONE) {
          num_active++;
        }
      });

      if (num_active == 0 && last_neighbor != -1) {
        this->set_value(last_neighbor, LABEL_ONE);
      }
    }
  }

  this->calculate_fitness();
}
*/

void Chromosome::fix_l(const ListGraph &graph) {
  boost::dynamic_bitset<> already_dominated;
  already_dominated.resize(this->size(), 0);

  std::vector<int> vertices(this->size());

  // Preenche o vetor com os valores 0, 1, ..., size()-1
  std::iota(vertices.begin(), vertices.end(), 0);

  // Gera um motor de números aleatórios
  std::random_device rd;
  std::mt19937 gen(rd());

  // Embaralha o vetor
  std::shuffle(vertices.begin(), vertices.end(), gen);

  for (const int &u : vertices) {
    if (this->get_value(u) == LABEL_ZERO && !already_dominated[u]) {
      unsigned short num_active = 0;
      bool has_neighbor_with_label_two = false;
      int vertex_with_label_one = -1;
      int last_neighbor = -1;

      for (const int &w : graph.get_neighbors(u)) {
        last_neighbor = w;
        if (this->get_value(w) == LABEL_ONE) {
          num_active++;
          vertex_with_label_one = w;
        } else if (this->get_value(w) == LABEL_TWO) {
          num_active++;
          has_neighbor_with_label_two = true;
          break; // pode parar o laço
        }
      }

      if (num_active == 0 && last_neighbor != -1) {
        this->set_value(last_neighbor, LABEL_TWO);
      } else if (!has_neighbor_with_label_two && vertex_with_label_one != -1) {
        this->set_value(vertex_with_label_one, LABEL_TWO);
      }
      already_dominated[u] = 1;
    } else if (this->get_value(u) == LABEL_TWO) {
      unsigned short num_active = 0;
      int last_neighbor = -1;

      for (const int &w : graph.get_neighbors(u)) {
        last_neighbor = w;
        if (w >= 0 && w < static_cast<int>(already_dominated.size())) {
          already_dominated[w] = 1;
        }
        if (this->get_value(w) >= LABEL_ONE) {
          num_active++;
        }
      }

      if (num_active == 0 && last_neighbor != -1) {
        this->set_value(last_neighbor, LABEL_ONE);
      }
    } else if (this->get_value(u) == LABEL_ONE && !already_dominated[u]) {
      unsigned short num_active = 0;
      int last_neighbor = -1;

      for (const int &w : graph.get_neighbors(u)) {
        last_neighbor = w;
        if (this->get_value(w) >= LABEL_ONE) {
          num_active++;
          break; // pode parar o laço
        }
      }

      if (num_active == 0 && last_neighbor != -1) {
        this->set_value(last_neighbor, LABEL_ONE);
        if (last_neighbor >= 0 &&
            last_neighbor < static_cast<int>(already_dominated.size())) {
          already_dominated[last_neighbor] = 1;
        }
      }
      already_dominated[u] = 1;
    }
  }

  this->calculate_fitness();
}

void Chromosome::fix_m(const MatrixGraph &graph) {
  boost::dynamic_bitset<> already_dominated;
  already_dominated.resize(this->size(), 0);

  std::vector<int> vertices(this->size());

  // Preenche o vetor com os valores 0, 1, ..., size()-1
  std::iota(vertices.begin(), vertices.end(), 0);

  // Gera um motor de números aleatórios
  std::random_device rd;
  std::mt19937 gen(rd());

  // Embaralha o vetor
  std::shuffle(vertices.begin(), vertices.end(), gen);

  for (const int &u : vertices) {
    if (this->get_value(u) == LABEL_ZERO && !already_dominated[u]) {
      unsigned short num_active = 0;
      bool has_neighbor_with_label_two = false;
      int vertex_with_label_one = -1;
      int last_neighbor = -1;

      auto &neighbors = graph.get_neighbors(u);
      for (int w = 0; w < static_cast<int>(neighbors.size()); ++w) {
        if (neighbors[w] == 1) {
          last_neighbor = w;
          if (this->get_value(w) == LABEL_ONE) {
            num_active++;
            vertex_with_label_one = w;
          } else if (this->get_value(w) == LABEL_TWO) {
            num_active++;
            has_neighbor_with_label_two = true;
            break; // pode parar de percorrer
          }
        }
      }

      if (num_active == 0 && last_neighbor != -1) {
        this->set_value(last_neighbor, LABEL_TWO);
      } else if (!has_neighbor_with_label_two && vertex_with_label_one != -1) {
        this->set_value(vertex_with_label_one, LABEL_TWO);
      }
      already_dominated[u] = 1;
    } else if (this->get_value(u) == LABEL_TWO) {
      unsigned short num_active = 0;
      int last_neighbor = -1;

      auto &neighbors = graph.get_neighbors(u);
      for (int w = 0; w < static_cast<int>(neighbors.size()); ++w) {
        if (neighbors[w] == 1) {
          last_neighbor = w;
          if (w >= 0 && w < static_cast<int>(already_dominated.size())) {
            already_dominated[w] = 1;
          }
          if (this->get_value(w) >= LABEL_ONE) {
            num_active++;
          }
        }
      }

      if (num_active == 0 && last_neighbor != -1) {
        this->set_value(last_neighbor, LABEL_ONE);
      }
    } else if (this->get_value(u) == LABEL_ONE && !already_dominated[u]) {
      unsigned short num_active = 0;
      int last_neighbor = -1;

      auto &neighbors = graph.get_neighbors(u);
      for (int w = 0; w < static_cast<int>(neighbors.size()); ++w) {
        if (neighbors[w] == 1) {
          last_neighbor = w;
          if (this->get_value(w) >= LABEL_ONE) {
            num_active++;
            break; // pode parar o laço
          }
        }
      }

      if (num_active == 0 && last_neighbor != -1) {
        this->set_value(last_neighbor, LABEL_ONE);
        if (last_neighbor >= 0 &&
            last_neighbor < static_cast<int>(already_dominated.size())) {
          already_dominated[last_neighbor] = 1;
        }
      }
      already_dominated[u] = 1;
    }
  }

  this->calculate_fitness();
}

// void Chromosome::fix_l(const ListGraph &graph) {
//   boost::dynamic_bitset<> already_dominated;
//   already_dominated.resize(this->size(), 0);
//
//   std::vector<int> vertices(this->size());
//
//   // Preenche o vetor com os valores 0, 1, ..., size()-1
//   std::iota(vertices.begin(), vertices.end(), 0);
//
//   // Gera um motor de números aleatórios
//   std::random_device rd;
//   std::mt19937 gen(rd());
//
//   // Embaralha o vetor
//   std::shuffle(vertices.begin(), vertices.end(), gen);
//
//   //for (size_t u = 0; u < this->size(); ++u) {
//   for(const int& u : vertices) {
//     if (this->get_value(u) == LABEL_ZERO && !already_dominated[u]) {
//       unsigned short num_active = 0;
//       bool has_neighbor_with_label_two = false;
//       int vertex_with_label_one = -1;
//       int last_neighbor = -1;
//
//       for(const int& w : graph.get_neighbors(u)) {
//         last_neighbor = w;
//         if (this->get_value(w) == LABEL_ONE) {
//           num_active++;
//           vertex_with_label_one = w;
//         } else if (this->get_value(w) == LABEL_TWO) {
//           num_active++;
//           has_neighbor_with_label_two = true;
//           break; // pode parar o laço
//         }
//       }
//
//       if (num_active == 0 && last_neighbor != -1) {
//         this->set_value(last_neighbor, LABEL_TWO);
//       } else if (!has_neighbor_with_label_two && vertex_with_label_one != -1)
//       {
//         this->set_value(vertex_with_label_one, LABEL_TWO);
//       }
//       already_dominated[u] = 1;
//     }
//     else if (this->get_value(u) == LABEL_TWO) {
//       unsigned short num_active = 0;
//       int last_neighbor = -1;
//
//       for(const int& w : graph.get_neighbors(u)) {
//         last_neighbor = w;
//         already_dominated[w] = 1;
//         if (this->get_value(w) >= LABEL_ONE) {
//           num_active++;
//         }
//       }
//
//       if (num_active == 0 && last_neighbor != -1) {
//         this->set_value(last_neighbor, LABEL_ONE);
//       }
//     }
//     else if (this->get_value(u) == LABEL_ONE && !already_dominated[u]) {
//       unsigned short num_active = 0;
//       int last_neighbor = -1;
//
//       for(const int& w : graph.get_neighbors(u)) {
//         last_neighbor = w;
//         if (this->get_value(w) >= LABEL_ONE) {
//           num_active++;
//           break; // pode parar o laço
//         }
//       }
//
//       if (num_active == 0 && last_neighbor != -1) {
//         this->set_value(last_neighbor, LABEL_ONE);
//       }
//       already_dominated[u] = 1;
//       already_dominated[last_neighbor] = 1;
//     }
//   }
//
//   this->calculate_fitness();
// }
//
//
// void Chromosome::fix_m(const MatrixGraph &graph) {
//   boost::dynamic_bitset<> already_dominated;
//   already_dominated.resize(this->size(), 0);
//
//   std::vector<int> vertices(this->size());
//
//   // Preenche o vetor com os valores 0, 1, ..., size()-1
//   std::iota(vertices.begin(), vertices.end(), 0);
//
//   // Gera um motor de números aleatórios
//   std::random_device rd;
//   std::mt19937 gen(rd());
//
//   // Embaralha o vetor
//   std::shuffle(vertices.begin(), vertices.end(), gen);
//
//   //for (size_t u = 0; u < this->size(); ++u) {
//   for(const int& u : vertices) {
//     if (this->get_value(u) == LABEL_ZERO && !already_dominated[u]) {
//       unsigned short num_active = 0;
//       bool has_neighbor_with_label_two = false;
//       int vertex_with_label_one = -1;
//       int last_neighbor = -1;
//
//       auto& neighbors = graph.get_neighbors(u);
//       for(int w = 0; w < (int)neighbors.size(); ++w) {
//         if(neighbors[w] == 1) {
//           last_neighbor = w;
//           if (this->get_value(w) == LABEL_ONE) {
//             num_active++;
//             vertex_with_label_one = w;
//           } else if (this->get_value(w) == LABEL_TWO) {
//             num_active++;
//             has_neighbor_with_label_two = true;
//             break; // pode parar de percorrer
//           }
//         }
//       }
//
//       if (num_active == 0 && last_neighbor != -1) {
//         this->set_value(last_neighbor, LABEL_TWO);
//       } else if (!has_neighbor_with_label_two && vertex_with_label_one != -1)
//       {
//         this->set_value(vertex_with_label_one, LABEL_TWO);
//       }
//       already_dominated[u] = 1;
//     }
//     else if (this->get_value(u) == LABEL_TWO) {
//       unsigned short num_active = 0;
//       int last_neighbor = -1;
//
//       auto& neighbors = graph.get_neighbors(u);
//       for(int w = 0; w < (int)neighbors.size(); ++w) {
//         if(neighbors[w] == 1) {
//           last_neighbor = w;
//           already_dominated[w] = 1;
//           if (this->get_value(w) >= LABEL_ONE) {
//             num_active++;
//           }
//         }
//       }
//
//       if (num_active == 0 && last_neighbor != -1) {
//         this->set_value(last_neighbor, LABEL_ONE);
//       }
//     }
//     else if (this->get_value(u) == LABEL_ONE && !already_dominated[u]) {
//       unsigned short num_active = 0;
//       int last_neighbor = -1;
//
//       auto& neighbors = graph.get_neighbors(u);
//       for(int w = 0; w < (int)neighbors.size(); ++w) {
//         if(neighbors[w] == 1) {
//           last_neighbor = w;
//           if (this->get_value(w) >= LABEL_ONE) {
//             num_active++;
//             break; // pode parar o laço
//           }
//         }
//       }
//
//       if (num_active == 0 && last_neighbor != -1) {
//         this->set_value(last_neighbor, LABEL_ONE);
//       }
//       already_dominated[u] = 1;
//       already_dominated[last_neighbor] = 1;
//     }
//   }
//
//   this->calculate_fitness();
// }
