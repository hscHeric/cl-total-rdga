#include "HeuristicGenerators.hpp"
#include "Chromosome.hpp"
#include "ListGraph.hpp"
#include "MatrixGraph.hpp"
#include <vector>

Chromosome HeuristicGenerators::h1_l(const ListGraph &graph) {
  // 1. Criar uma cópia do grafo G original
  ListGraph H(graph);

  // Inicializar o cromossomo com o tamanho igual à ordem
  // do grafo e valor padrão 0
  Chromosome chromosome(graph.order(), 0);

  // 2. Enquanto tiver vértices em H faça
  // No inicio de cada iteracao desse laço garantimos
  // que H não tem vértices isolados
  while (H.order() > 0) {
    // 3. Escolha aleatoriamente um vértice v qualquer de H
    int v = H.choose_rng();

    // 4. Faça f(v) = 2
    chromosome.set_value(v, LABEL_TWO);

    // 5. Pegue um vizinho u de v e faça f(u) = 1 (um vizinho aleatório)
    auto &neighbors = H.get_neighbors(v);

    // Gera um motor de números aleatórios
    std::random_device rd;
    std::mt19937 gen(rd());

    // Criar uma distribuição uniforme no intervalo dos índices do vetor
    std::uniform_int_distribution<size_t> dist(0, neighbors.size() - 1);

    // Escolher um vizinho aleatório
    int random_index = dist(gen);
    int u = neighbors[random_index];

    chromosome.set_value(u, LABEL_ONE);

    // 6. Para todos os demais vizinhos w de v faça f(w) = 0
    for (size_t i = 0; i < neighbors.size(); ++i) {
      if (static_cast<int>(i) != random_index) {
        chromosome.set_value(neighbors[i], LABEL_ZERO);
      }
    }

    // 7. Remova do grafo H o vértice v e todos os seus vizinhos
    for (int j = 0; j < (int)neighbors.size(); ++j) {
      H.remove_vertex(neighbors[j]);
    }
    H.remove_vertex(v);

    // 8. Se restarem vértices isolados em H atribuímos o
    // rótulo 1 a todos eles.
    auto isolatedVertices = H.get_isolated_vertices();

    for (const int &z : isolatedVertices) {
      // 9. Faça f(z) = 1
      chromosome.set_value(z, LABEL_ONE);

      // 10. Seja x um vizinho de z no grafo G
      // Escolhemos o primeiro na lista de adjacencia
      auto &g_neighbors = graph.get_neighbors(z);
      int x = g_neighbors[0];

      // 11. Mude a cor do vértice x para f(x) = 1
      chromosome.set_value(x, LABEL_ONE);

      // 12. Remova z do grafo H
      H.remove_vertex(z);
    }
  }

  // Calcula o fitness do cromossomo antes de retorná-lo
  chromosome.calculate_fitness();

  return chromosome;
}

Chromosome HeuristicGenerators::h1_m(const MatrixGraph &graph) {
  // 1. Criar uma cópia do grafo G original
  MatrixGraph H(graph);

  // Inicializar o cromossomo com o tamanho igual à ordem
  // do grafo e valor padrão 0
  Chromosome chromosome(graph.order(), 0);

  // 2. Enquanto tiver vértices em H faça
  // No inicio de cada iteracao desse laço garantimos
  // que H não tem vértices isolados
  while (H.order() > 0) {
    // 3. Escolha aleatoriamente um vértice v qualquer de H
    int v = H.choose_rng();

    // 4. Faça f(v) = 2
    chromosome.set_value(v, LABEL_TWO);

    // 5. Pegue um vizinho u de v e faça f(u) = 1 (o primeiro vizinho que
    // aparecer)
    auto &neighbors = H.get_neighbors(v);
    int u = -1;
    for (u = 0; u < static_cast<int>(neighbors.size()); ++u) {
      if (neighbors[u] == 1) {
        chromosome.set_value(u, LABEL_ONE);
        break;
      }
    }

    // 6. Para todos os demais vizinhos w de v faça f(w) = 0
    for (u = u + 1; u < static_cast<int>(neighbors.size()); ++u) {
      if (neighbors[u] == 1) {
        chromosome.set_value(u, LABEL_ZERO);
      }
    }

    // 7. Remova do grafo H o vértice v e todos os seus vizinhos
    for (int w = 0; w < static_cast<int>(neighbors.size()); ++w) {
      if (neighbors[w] == 1) {
        H.remove_vertex(w);
      }
    }
    H.remove_vertex(v);

    // 8. Se restarem vértices isolados em H atribuímos o
    // rótulo 1 a todos eles.
    auto isolatedVertices = H.get_isolated_vertices();

    for (const int &z : isolatedVertices) {
      // 9. Faça f(z) = 1
      chromosome.set_value(z, LABEL_ONE);

      // 10. Seja x um vizinho de z no grafo G.
      // Escolhemos o primeiro vizinho que aparecer.
      // 11. Mude a cor do vértice x para f(x) = 1
      auto &g_neighbors = graph.get_neighbors(z);
      for (int x = 0; x < static_cast<int>(g_neighbors.size()); ++x) {
        if (g_neighbors[x] == 1) {
          chromosome.set_value(x, LABEL_ONE);
          break;
        }
      }

      // 12. Remova z do grafo H
      H.remove_vertex(z);
    }
  }

  // Calcula o fitness do cromossomo antes de retorná-lo
  chromosome.calculate_fitness();

  return chromosome;
}

/*
Chromosome HeuristicGenerators::h2(const Graph &graph) {
  std::unique_ptr<Graph> H(copyGraph(graph));

  Chromosome chromosome(graph.order(), 0);

  while (H->order() > 0) {
    int v = -1;
    int maxDegree = -1;
    H->for_each_vertex([&H, &v, &maxDegree](int vertex) {
      int currentDegree = H->degree(vertex);
      if (currentDegree > maxDegree) {
        maxDegree = currentDegree;
        v = vertex;
      }
    });

    chromosome.set_value(v, LABEL_TWO);

    std::vector<int> neighbors;
    H->for_each_neighbor(
        v, [&neighbors](int neighbor) { neighbors.push_back(neighbor); });
    // 5. Pegue um vizinho u de v e faça f(u) = 1 (o primeiro da lista de
    // adjacência)
    if (!neighbors.empty()) {
      int u = neighbors[0];
      chromosome.set_value(u, LABEL_ONE);
      // 6. Para todos os demais vizinhos w de v faça f(w) = 0
      for (size_t i = 1; i < neighbors.size(); ++i) {
        chromosome.set_value(neighbors[i], LABEL_ZERO);
      }
    }
    // 7. Remova do grafo H o vértice v e todos os seus vizinhos
    for (int neighbor : neighbors) {
      if (H->contains(neighbor)) {
        H->remove_vertex(neighbor);
      }
    }
    if (H->contains(v)) {
      H->remove_vertex(v);
    }
    // 8. Enquanto houver vértice isolado z em H faça
    auto isolatedVertices = getIsolatedVertices(*H);
    while (!isolatedVertices.empty()) {
      // Pega um vértice isolado qualquer
      int z = *isolatedVertices.begin();
      // 9. Faça f(z) = 1
      chromosome.set_value(z, LABEL_ONE);
      // 10. Se no grafo G, z não tiver um vizinho com rótulo 1 então:
      bool hasLabel1Neighbor = false;
      graph.for_each_neighbor(
          z, [&hasLabel1Neighbor, &chromosome](int neighbor) {
            if (chromosome.get_value(neighbor) == LABEL_ONE) {
              hasLabel1Neighbor = true;
            }
          });

      if (!hasLabel1Neighbor) {
        // 11. Seja x um vizinho de z no grafo G
        std::vector<int> g_neighbors;
        graph.for_each_neighbor(z, [&g_neighbors](int neighbor) {
          g_neighbors.push_back(neighbor);
        });
        if (!g_neighbors.empty()) {
          // Escolhe um vizinho qualquer (o primeiro)
          int x = g_neighbors[0];
          // 12. Mude a cor do vértice x para f(x) = 1
          chromosome.set_value(x, LABEL_ONE);
        }
      }
      // 13. Remova z do grafo H
      H->remove_vertex(z);
      // Atualize a lista de vértices isolados
      isolatedVertices = getIsolatedVertices(*H);
    }
  }
  // Calcula o fitness do cromossomo antes de retorná-lo
  chromosome.calculate_fitness();
  return chromosome;
}

Chromosome HeuristicGenerators::h3(const Graph &graph) {
  std::unique_ptr<Graph> H(copyGraph(graph));
  Chromosome chromosome(graph.order(), 0);
  while (H->order() > 0) {
    // Encontrar o vértice v com maior grau
    int v = -1;
    int maxDegree = -1;
    H->for_each_vertex([&H, &v, &maxDegree](int vertex) {
      int currentDegree = H->degree(vertex);
      if (currentDegree > maxDegree) {
        maxDegree = currentDegree;
        v = vertex;
      }
    });
    chromosome.set_value(v, LABEL_TWO);

    // Pegar todos os vizinhos de v
    std::vector<int> neighbors;
    H->for_each_neighbor(
        v, [&neighbors](int neighbor) { neighbors.push_back(neighbor); });

    if (!neighbors.empty()) {
      // A diferença: escolher o vizinho u com o maior grau para receber label 1
      int u = -1;
      int maxNeighborDegree = -1;
      for (int neighbor : neighbors) {
        int neighborDegree = H->degree(neighbor);
        if (neighborDegree > maxNeighborDegree) {
          maxNeighborDegree = neighborDegree;
          u = neighbor;
        }
      }

      // Atribuir label 1 ao vizinho com maior grau
      chromosome.set_value(u, LABEL_ONE);

      // Restante dos vizinhos recebem label 0
      for (int neighbor : neighbors) {
        if (neighbor != u) {
          chromosome.set_value(neighbor, LABEL_ZERO);
        }
      }
    }

    // Remover do grafo H o vértice v e todos os seus vizinhos
    for (int neighbor : neighbors) {
      if (H->contains(neighbor)) {
        H->remove_vertex(neighbor);
      }
    }
    if (H->contains(v)) {
      H->remove_vertex(v);
    }

    // Enquanto houver vértice isolado z em H
    auto isolatedVertices = getIsolatedVertices(*H);
    while (!isolatedVertices.empty()) {
      int z = *isolatedVertices.begin();
      chromosome.set_value(z, LABEL_ONE);

      // Verificar se z não tem um vizinho com rótulo 1 no grafo original
      bool hasLabel1Neighbor = false;
      graph.for_each_neighbor(
          z, [&hasLabel1Neighbor, &chromosome](int neighbor) {
            if (chromosome.get_value(neighbor) == LABEL_ONE) {
              hasLabel1Neighbor = true;
            }
          });

      if (!hasLabel1Neighbor) {
        // Encontrar um vizinho de z no grafo original G e atribuir label 1
        std::vector<int> g_neighbors;
        graph.for_each_neighbor(z, [&g_neighbors](int neighbor) {
          g_neighbors.push_back(neighbor);
        });

        if (!g_neighbors.empty()) {
          int x = g_neighbors[0];
          chromosome.set_value(x, LABEL_ONE);
        }
      }

      // Remover z do grafo H
      H->remove_vertex(z);
      isolatedVertices = getIsolatedVertices(*H);
    }
  }

  chromosome.calculate_fitness();
  return chromosome;
}

// Chromosome HeuristicGenerators::h4(const Graph &graph) {
//   std::unique_ptr<Graph> H(copyGraph(graph));
//   Chromosome chromosome(graph.order(), 0);
//
//   while (H->order() > 0) {
//     // Encontrar o vértice v com maior grau
//     int v = -1;
//     int maxDegree = -1;
//     H->for_each_vertex([&H, &v, &maxDegree](int vertex) {
//       int currentDegree = H->degree(vertex);
//       if (currentDegree > maxDegree) {
//         maxDegree = currentDegree;
//         v = vertex;
//       }
//     });
//     chromosome.set_value(v, LABEL_TWO);
//
//     // Pegar todos os vizinhos de v
//     std::vector<int> neighbors;
//     H->for_each_neighbor(
//         v, [&neighbors](int neighbor) { neighbors.push_back(neighbor); });
//
//     if (!neighbors.empty()) {
//       // A diferença: escolher o vizinho u com o maior grau para receber
//       label 1 int u = -1; int maxNeighborDegree = -1; for (int neighbor :
//       neighbors) {
//         int neighborDegree = H->degree(neighbor);
//         if (neighborDegree > maxNeighborDegree) {
//           maxNeighborDegree = neighborDegree;
//           u = neighbor;
//         }
//       }
//
//       // Atribuir label 1 ao vizinho com maior grau
//       chromosome.set_value(u, LABEL_ONE);
//
//       // Restante dos vizinhos recebem label 0
//       for (int neighbor : neighbors) {
//         if (neighbor != u) {
//           chromosome.set_value(neighbor, LABEL_ZERO);
//         }
//       }
//     }
//
//     // Remover do grafo H o vértice v e todos os seus vizinhos
//     for (int neighbor : neighbors) {
//       if (H->contains(neighbor)) {
//         H->remove_vertex(neighbor);
//       }
//     }
//     if (H->contains(v)) {
//       H->remove_vertex(v);
//     }
//
//     auto isolated_vertices = getIsolatedVertices(*H);
//     while (!isolated_vertices.empty()) {
//       std::set<int> S(isolated_vertices.begin(), isolated_vertices.end());
//       std::set<int> N_S;
//
//       for (auto s : S) {
//         graph.for_each_neighbor(s,
//                                 [&N_S](int neighbor) { N_S.insert(neighbor);
//                                 });
//       }
//
//       for (int z : N_S) {
//         int count = 0;
//         graph.for_each_neighbor(z, [&S, &count](int neighbor) {
//           if (S.count(neighbor))
//             count++;
//         });
//
//         if (count >= 2) {
//           chromosome.set_value(z, LABEL_TWO);
//           graph.for_each_neighbor(z, [&chromosome, &S](int neighbor) {
//             if (S.count(neighbor)) {
//               chromosome.set_value(neighbor, LABEL_ZERO);
//             }
//           });
//         }
//       }
//
//       // Atribuir rótulo 2 aos vértices restantes de N(S)
//       for (int z : N_S) {
//         if (chromosome.get_value(z) == LABEL_ZERO) {
//           chromosome.set_value(z, LABEL_TWO);
//         }
//       }
//
//       // Atribuir rótulo 0 aos vértices restantes de S
//       for (int s : S) {
//         if (chromosome.get_value(s) == LABEL_ZERO) {
//           chromosome.set_value(s, LABEL_ZERO);
//         }
//       }
//
//       // Remover S de H
//       for (int s : S) {
//         H->remove_vertex(s);
//       }
//
//       isolated_vertices = getIsolatedVertices(*H);
//     }
//   }
//
//   chromosome.calculate_fitness();
//   return chromosome;
// }

Chromosome HeuristicGenerators::h4(const Graph &graph) {
  std::unique_ptr<Graph> H(copyGraph(graph));
  Chromosome chromosome(graph.order(), 0);

  while (H->order() > 0) {
    int v = -1;
    int maxDegree = -1;
    H->for_each_vertex([&H, &v, &maxDegree](int vertex) {
      int currentDegree = H->degree(vertex);
      if (currentDegree > maxDegree) {
        maxDegree = currentDegree;
        v = vertex;
      }
    });

    chromosome.set_value(v, LABEL_TWO);

    std::vector<int> neighbors;
    H->for_each_neighbor(
        v, [&neighbors](int neighbor) { neighbors.push_back(neighbor); });

    std::sort(neighbors.begin(), neighbors.end(),
              [&H](int a, int b) { return H->degree(b) < H->degree(a); });

    if (!neighbors.empty()) {
      chromosome.set_value(neighbors[0], LABEL_ONE);

      for (size_t i = 1; i < neighbors.size(); ++i) {
        chromosome.set_value(neighbors[i], LABEL_ZERO);
      }
    }

    for (int neighbor : neighbors) {
      if (H->contains(neighbor)) {
        H->remove_vertex(neighbor);
      }
    }
    if (H->contains(v)) {
      H->remove_vertex(v);
    }

    while (true) {
      auto isolatedVertices = getIsolatedVertices(*H);
      if (isolatedVertices.empty()) {
        break;
      }

      std::set<int> neighborSet;
      for (int s : isolatedVertices) {
        graph.for_each_neighbor(
            s, [&neighborSet](int neighbor) { neighborSet.insert(neighbor); });
      }
      std::vector<int> neighbors(neighborSet.begin(), neighborSet.end());

      for (int z : neighbors) {
        int isolatedNeighborCount = 0;
        for (int s : isolatedVertices) {
          if (graph.contains(s, z)) {
            isolatedNeighborCount++;
          }
        }

        if (isolatedNeighborCount >= 2) {
          chromosome.set_value(z, LABEL_TWO);

          for (int s : isolatedVertices) {
            if (graph.contains(s, z)) {
              chromosome.set_value(s, LABEL_ZERO);
            }
          }
        } else {
          chromosome.set_value(z, LABEL_TWO);
        }
      }

      for (int s : isolatedVertices) {
        if (chromosome.get_value(s) == LABEL_ZERO) {
          chromosome.set_value(s, LABEL_ZERO);
        }
      }

      for (int s : isolatedVertices) {
        H->remove_vertex(s);
      }
    }
  }

  chromosome.calculate_fitness();
  return chromosome;
}
*/

Chromosome HeuristicGenerators::h5_l(const ListGraph &graph) {
  Chromosome chromosome(graph.order(), 1);

  chromosome.calculate_fitness();
  return chromosome;
}

Chromosome HeuristicGenerators::h5_m(const MatrixGraph &graph) {
  Chromosome chromosome(graph.order(), 1);

  chromosome.calculate_fitness();
  return chromosome;
}

Chromosome HeuristicGenerators::h2_l(const ListGraph &graph) {
  ListGraph H(graph);

  Chromosome chromosome(graph.order(), 0);

  while (H.order() > 0) {
    int v = -1;
    int maxDegree = -1;

    auto vertices = H.get_vertices();
    for (int vertex : vertices) {
      int currentDegree = H.degree(vertex);
      if (currentDegree > maxDegree) {
        maxDegree = currentDegree;
        v = vertex;
      }
    }

    chromosome.set_value(v, LABEL_TWO);

    const std::vector<int> &neighbors = H.get_neighbors(v);

    if (!neighbors.empty()) {
      int u = neighbors[0];
      chromosome.set_value(u, LABEL_ONE);

      for (size_t i = 1; i < neighbors.size(); ++i) {
        chromosome.set_value(neighbors[i], LABEL_ZERO);
      }
    }

    std::vector<int> neighbors_copy = neighbors;
    for (int neighbor : neighbors_copy) {
      if (H.contains(neighbor)) {
        H.remove_vertex(neighbor);
      }
    }

    if (H.contains(v)) {
      H.remove_vertex(v);
    }

    auto isolatedVertices = H.get_isolated_vertices();

    while (!isolatedVertices.empty()) {
      int z = *isolatedVertices.begin();

      chromosome.set_value(z, LABEL_ONE);

      bool hasLabel1Neighbor = false;

      for (int neighbor : graph.get_neighbors(z)) {
        if (chromosome.get_value(neighbor) == LABEL_ONE) {
          hasLabel1Neighbor = true;
          break;
        }
      }

      if (!hasLabel1Neighbor) {
        const std::vector<int> &g_neighbors = graph.get_neighbors(z);
        if (!g_neighbors.empty()) {
          int x = g_neighbors[0];
          chromosome.set_value(x, LABEL_ONE);
        }
      }

      H.remove_vertex(z);

      isolatedVertices = H.get_isolated_vertices();
    }
  }

  chromosome.calculate_fitness();
  return chromosome;
}

Chromosome HeuristicGenerators::h2_m(const MatrixGraph &graph) {
  MatrixGraph H(graph);

  Chromosome chromosome(graph.order(), 0);

  while (H.order() > 0) {
    int v = -1;
    int maxDegree = -1;

    auto vertices = H.get_vertices();
    for (int vertex : vertices) {
      int currentDegree = H.degree(vertex);
      if (currentDegree > maxDegree) {
        maxDegree = currentDegree;
        v = vertex;
      }
    }

    chromosome.set_value(v, LABEL_TWO);

    const BitSet &neighborBits = H.get_neighbors(v);
    std::vector<int> neighbors;

    for (size_t i = 0; i < neighborBits.size(); ++i) {
      if (neighborBits[i]) {
        neighbors.push_back(static_cast<int>(i));
      }
    }

    if (!neighbors.empty()) {
      int u = neighbors[0];
      chromosome.set_value(u, LABEL_ONE);

      for (size_t i = 1; i < neighbors.size(); ++i) {
        chromosome.set_value(neighbors[i], LABEL_ZERO);
      }
    }

    std::vector<int> neighbors_copy = neighbors;
    for (int neighbor : neighbors_copy) {
      if (H.contains(neighbor)) {
        H.remove_vertex(neighbor);
      }
    }

    if (H.contains(v)) {
      H.remove_vertex(v);
    }

    auto isolatedVertices = H.get_isolated_vertices();

    while (!isolatedVertices.empty()) {
      int z = *isolatedVertices.begin();

      chromosome.set_value(z, LABEL_ONE);

      bool hasLabel1Neighbor = false;
      const BitSet &zNeighbors = graph.get_neighbors(z);

      for (size_t i = 0; i < zNeighbors.size(); ++i) {
        if (zNeighbors[i] &&
            chromosome.get_value(static_cast<int>(i)) == LABEL_ONE) {
          hasLabel1Neighbor = true;
          break;
        }
      }

      if (!hasLabel1Neighbor) {
        const BitSet &gNeighbors = graph.get_neighbors(z);
        for (size_t i = 0; i < gNeighbors.size(); ++i) {
          if (gNeighbors[i]) {
            int x = static_cast<int>(i);
            chromosome.set_value(x, LABEL_ONE);
            break;
          }
        }
      }

      H.remove_vertex(z);

      isolatedVertices = H.get_isolated_vertices();
    }
  }

  chromosome.calculate_fitness();
  return chromosome;
}
