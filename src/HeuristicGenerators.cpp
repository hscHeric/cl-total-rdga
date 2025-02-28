#include "HeuristicGenerators.hpp"
#include "Chromosome.hpp"
#include "Graph.hpp"
#include "ListGraph.hpp"
#include "MatrixGraph.hpp"
#include <memory>
#include <vector>

Graph *HeuristicGenerators::copyGraph(const Graph &graph) {
  int n = graph.order();
  auto density = 0.0;
  if (n <= 1)
    density = 0.0;

  int m = graph.size();

  density = static_cast<double>(2 * m) / (n * (n - 1));

  const double DENSITY_THRESHOLD = 0.5;

  Graph *graph_copy = nullptr;
  if (density >= DENSITY_THRESHOLD) {
    graph_copy = new MatrixGraph(graph.order());
  } else {
    graph_copy = new ListGraph(graph.order());
  }

  graph.for_each_edge(
      [&graph_copy](int u, int v) { graph_copy->add_edge(u, v); });

  return graph_copy;
}

bool HeuristicGenerators::isIsolatedVertex(const Graph &graph, int vertex) {
  return graph.contains(vertex) && graph.degree(vertex) == 0;
}

std::unordered_set<int>
HeuristicGenerators::getIsolatedVertices(const Graph &graph) {
  std::unordered_set<int> isolatedVertices;

  graph.for_each_vertex([&](int v) {
    if (graph.contains(v) && graph.degree(v) == 0) {
      isolatedVertices.insert(v);
    }
  });

  return isolatedVertices;
}

Chromosome HeuristicGenerators::h1(const Graph &graph) {
  // 1. Criar uma cópia do grafo G original
  std::unique_ptr<Graph> H(copyGraph(graph));

  // Inicializar o cromossomo com o tamanho igual à ordem do grafo e valor
  // padrão 0
  Chromosome chromosome(graph.order(), 0);

  // 2. Enquanto tiver vértices em H faça
  while (H->order() > 0) {
    // 3. Escolha aleatoriamente um vértice v qualquer de H
    int v = H->choose_rng();

    // 4. Faça f(v) = 2
    chromosome.set_value(v, LABEL_TWO);

    // Armazenar vizinhos de v para processamento e remoção posterior
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

      // 10. Seja x um vizinho qualquer de z no grafo G
      std::vector<int> g_neighbors;
      graph.for_each_neighbor(
          z, [&g_neighbors](int neighbor) { g_neighbors.push_back(neighbor); });

      if (!g_neighbors.empty()) {
        // Escolhe um vizinho qualquer (o primeiro)
        int x = g_neighbors[0];

        // 11. Mude a cor do vértice x para f(x) = 1
        chromosome.set_value(x, LABEL_ONE);
      }

      // 12. Remova z do grafo H
      H->remove_vertex(z);

      // Atualize a lista de vértices isolados
      isolatedVertices = getIsolatedVertices(*H);
    }
  }

  // Calcula o fitness do cromossomo antes de retorná-lo
  chromosome.calculate_fitness();

  return chromosome;
}

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
