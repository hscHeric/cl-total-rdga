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

std::unique_ptr<Graph> load_and_normalize_graph(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Não foi possível abrir o arquivo: " + filename);
  }

  // Primeira passagem: coletar todos os vértices únicos
  std::set<int> vertices;
  std::vector<std::pair<int, int>> edges;
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
      edges.emplace_back(u, v);
    }
  }

  if (vertices.empty()) {
    throw std::runtime_error("Nenhum vértice encontrado na entrada");
  }

  // Identificar lacunas nos IDs de vértices
  std::vector<int> missing_vertices;
  if (!vertices.empty()) {
    int min_vertex = *vertices.begin();
    int max_vertex = *vertices.rbegin();

    for (int i = min_vertex; i <= max_vertex; i++) {
      if (vertices.find(i) == vertices.end()) {
        missing_vertices.push_back(i);
      }
    }
  }

  std::unordered_map<int, int> id_map;
  int new_id = 0;
  for (int original_id : vertices) {
    id_map[original_id] = new_id++;
  }

  int num_vertices = vertices.size();
  double max_possible_edges = num_vertices * (num_vertices - 1) / 2.0;
  std::set<std::pair<int, int>> unique_edges;

  for (const auto &edge : edges) {
    int norm_u = id_map[edge.first];
    int norm_v = id_map[edge.second];
    if (norm_u > norm_v) {
      std::swap(norm_u, norm_v);
    }
    unique_edges.insert({norm_u, norm_v});
  }

  double density = unique_edges.size() / max_possible_edges;

  // Imprimir estatísticas
  std::cout << "Estatísticas do grafo normalizado:" << std::endl;
  std::cout << "  Vértices originais: " << vertices.size() << std::endl;
  std::cout << "  Intervalo original: [" << *vertices.begin() << ", "
            << *vertices.rbegin() << "]" << std::endl;

  if (!missing_vertices.empty()) {
    std::cout << "  Vértices faltantes: ";
    for (size_t i = 0; i < missing_vertices.size(); i++) {
      std::cout << missing_vertices[i];
      if (i < missing_vertices.size() - 1) {
        std::cout << ", ";
      }
    }
    std::cout << std::endl;
    std::cout << "  Total de lacunas: " << missing_vertices.size() << std::endl;
  } else {
    std::cout << "  Não há lacunas na numeração dos vértices" << std::endl;
  }

  std::cout << "  Mapeamento aplicado: ";
  for (const auto &[orig, norm] : id_map) {
    std::cout << orig << "->" << norm << " ";
  }
  std::cout << std::endl;
  std::cout << "  Vértices normalizados: " << num_vertices << std::endl;
  std::cout << "  Arestas únicas: " << unique_edges.size() << std::endl;
  std::cout << "  Densidade: " << density << std::endl;
  std::cout << "  Implementação selecionada: "
            << (density > 0.3 ? "MatrixGraph" : "ListGraph") << std::endl;

  // Criar grafo com implementação adequada
  std::unique_ptr<Graph> graph;
  if (density > 0.3) {
    graph = std::make_unique<MatrixGraph>(num_vertices);
  } else {
    graph = std::make_unique<ListGraph>(num_vertices);
  }

  // Adicionar arestas normalizadas
  for (const auto &edge : unique_edges) {
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
  auto graph = load_and_normalize_graph(params.file_path);
  HeuristicGenerators heuristicHandle;

  if (graph->order() == 0) {
    std::cerr << "Error: The graph has no nodes." << std::endl;
    return 1;
  }

  // Test h1 heuristic
  std::cout << "Testing h1 heuristic:\n";
  auto solution_h1 = heuristicHandle.h1(*graph);
  auto is_valid_h1 = valid_totalrd(graph.get(), solution_h1);
  std::cout << "Solution h1 fitness: " << solution_h1.get_fitness()
            << std::endl;
  if (is_valid_h1) {
    std::cout << "h1 solution is valid\n";
  } else {
    std::cout << "h1 solution is not valid\n";
  }

  // Test h2 heuristic
  std::cout << "\nTesting h2 heuristic:\n";
  auto solution_h2 = heuristicHandle.h2(*graph);
  auto is_valid_h2 = valid_totalrd(graph.get(), solution_h2);
  std::cout << "Solution h2 fitness: " << solution_h2.get_fitness()
            << std::endl;
  if (is_valid_h2) {
    std::cout << "h2 solution is valid\n";
  } else {
    std::cout << "h2 solution is not valid\n";
  }

  std::cout << "\nTesting h3 heuristic:\n";
  auto solution_h3 = heuristicHandle.h2(*graph);
  auto is_valid_h3 = valid_totalrd(graph.get(), solution_h2);
  std::cout << "Solution h3 fitness: " << solution_h2.get_fitness()
            << std::endl;
  if (is_valid_h3) {
    std::cout << "h3 solution is valid\n";
  } else {
    std::cout << "h3 solution is not valid\n";
  }

  // Compare both solutions
  std::cout << "\nComparison:\n";
  if (solution_h1.get_fitness() < solution_h2.get_fitness()) {
    std::cout << "h1 found a better solution\n";
  } else if (solution_h1.get_fitness() > solution_h2.get_fitness()) {
    std::cout << "h2 found a better solution\n";
  } else if (solution_h3.get_fitness() < solution_h2.get_fitness() &&
             solution_h3.get_fitness() < solution_h1.get_fitness()) {
    std::cout << "h3 found a better solution\n";

  } else {
    std::cout << "Both heuristics found solutions with equal fitness\n";
  }

  return 0;
}
