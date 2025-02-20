#include "Chromosome.hpp"
#include "Graph.hpp"
#include "HeuristicGenerators.hpp"
#include "ListGraph.hpp"
#include "MatrixGraph.hpp"
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <memory>
#include <set>
#include <string>

struct AlgorithmParams {
  size_t max_stagnant = 100;
  size_t generations = 1000;
  size_t tournament_size = 5;
  double crossover_rate = 0.9;
  double population_factor = 1.5;
  std::string file_path;
  size_t trials = 1;
  std::string output_file = "results.csv";
};

struct TrialResult {
  std::string graph_name;
  size_t node_count;
  size_t edge_count;
  int fitness;
  uint64_t elapsed_micros;
};

std::unique_ptr<Graph> createOptimalGraph(int numVertices, double density) {
  // eficiente Para grafos esparsos, ListGraph é melhor
  if (density > 0.3) {
    return std::make_unique<MatrixGraph>(numVertices);
  } else {
    return std::make_unique<ListGraph>(numVertices);
  }
}

std::unique_ptr<Graph> load_graph(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Não foi possível abrir o arquivo: " + filename);
  }

  std::set<int> vertices;
  std::vector<std::pair<int, int>> validEdges;
  std::string line;

  while (std::getline(file, line)) {
    std::istringstream iss(line);
    int u, v;
    if (!(iss >> u >> v)) {
      continue;
    }

    vertices.insert(u);
    vertices.insert(v);

    if (u != v) {
      validEdges.emplace_back(u, v);
    }
  }

  if (vertices.empty()) {
    throw std::runtime_error("Nenhum vértice encontrado na entrada");
  }

  int maxVertex = *vertices.rbegin();
  int numVertices = maxVertex + 1;

  std::set<std::pair<int, int>> uniqueEdges;
  for (const auto &edge : validEdges) {
    int u = edge.first;
    int v = edge.second;
    if (u > v)
      std::swap(u, v); // Normaliza a aresta
    uniqueEdges.insert({u, v});
  }

  double maxPossibleEdges = numVertices * (numVertices - 1) / 2.0;
  double density = uniqueEdges.size() / maxPossibleEdges;

  std::cout << "Estatísticas do grafo:" << std::endl;
  std::cout << "  Vértices: " << numVertices << std::endl;
  std::cout << "  Arestas únicas: " << uniqueEdges.size() << std::endl;
  std::cout << "  Densidade: " << density << std::endl;
  std::cout << "  Implementação selecionada: "
            << (density > 0.3 ? "MatrixGraph" : "ListGraph") << std::endl;

  std::unique_ptr<Graph> graph = createOptimalGraph(numVertices, density);

  file.clear();
  file.seekg(0, std::ios::beg);

  for (const auto &edge : uniqueEdges) {
    graph->add_edge(edge.first, edge.second);
  }

  return graph;
}

AlgorithmParams parse_args(int argc, char *argv[]) {
  AlgorithmParams params;

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <graph_file> [options]\n"
              << "Options:\n"
              << "  --crossover VALUE\n"
              << "  --stagnation VALUE\n"
              << "  --generations VALUE\n"
              << "  --population VALUE\n"
              << "  --tournament VALUE\n"

              << "  --trials VALUE\n"
              << "  --output FILE\n";
    exit(1);
  }

  params.file_path = argv[1];

  for (int i = 2; i < argc; i++) {
    std::string arg = argv[i];
    if (arg == "--crossover" && i + 1 < argc) {
      params.crossover_rate = std::stod(argv[++i]);
    } else if (arg == "--stagnation" && i + 1 < argc) {
      params.max_stagnant = std::stoul(argv[++i]);
    } else if (arg == "--generations" && i + 1 < argc) {
      params.generations = std::stoul(argv[++i]);
    } else if (arg == "--population" && i + 1 < argc) {
      params.population_factor = std::stod(argv[++i]);
    } else if (arg == "--tournament" && i + 1 < argc) {
      params.tournament_size = std::stoul(argv[++i]);
    } else if (arg == "--trials" && i + 1 < argc) {
      params.trials = std::stoul(argv[++i]);
    } else if (arg == "--output" && i + 1 < argc) {
      params.output_file = argv[++i];
    } else {
      std::cerr << "Unknown argument: " << arg << std::endl;
      exit(1);
    }
  }

  return params;
}

bool valid_totalrd(const Graph *graph, const Chromosome &chromosome) {
  for (int v : graph->get_vertices()) {
    int label = chromosome.get_value(v);
    bool hasValidNeighbor = false;

    graph->for_each_neighbor(v, [&](int neighbor) {
      int neighbor_label = chromosome.get_value(neighbor);
      if ((label == 0 && neighbor_label == 2) ||
          (label > 0 && neighbor_label > 0)) {
        hasValidNeighbor = true;
      }
    });

    if (!hasValidNeighbor) {
      std::cout << "Solução inválida: vértice " << v << " com f(v) = " << label
                << " não possui vizinho válido." << std::endl;
      return false;
    }
  }
  return true;
}

int main(int argc, char *argv[]) {
  auto params = parse_args(argc, argv);
  auto graph = load_graph(params.file_path);
  HeuristicGenerators heuristicHandle;

  if (graph->order() == 0) {
    std::cerr << "Error: The graph has no nodes." << std::endl;
    return 1;
  }

  auto solution = heuristicHandle.h1(*graph);
  auto is_valid = valid_totalrd(graph.get(), solution);

  if (is_valid) {
    std::cout << "A solução é válida" << "\n";
  } else {
    std::cout << "A solução não é válida" << "\n";
  }
}
