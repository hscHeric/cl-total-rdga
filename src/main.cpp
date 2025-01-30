#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "../include/chromosome.hpp"
#include "../include/crossover.hpp"
#include "../include/dense_graph.hpp"
#include "../include/heuristics.hpp"
#include "../include/population.hpp"
#include "../include/selection.hpp"

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

void write_results_to_csv(const std::vector<TrialResult> &results,
                          const std::string &output_file) {
  std::ofstream file(output_file, std::ios::app);
  if (!file) {
    std::cerr << "Failed to open output file: " << output_file << std::endl;
    return;
  }

  if (file.tellp() == 0) {
    file << "graph_name,graph_order,graph_size,fitness_value,elapsed_time("
            "microsecond)\n";
  }

  for (const auto &result : results) {
    file << result.graph_name << "," << result.node_count << ","
         << result.edge_count << "," << result.fitness << ","
         << result.elapsed_micros << "\n";
  }
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

DenseGraph load_graph_from_file(const std::string &file_path) {
  std::ifstream file(file_path);
  if (!file) {
    std::cerr << "Error: Cannot open file " << file_path << std::endl;
    exit(1);
  }

  int num_vertices;
  if (!(file >> num_vertices) || num_vertices <= 0) {
    std::cerr << "Error: Invalid number of vertices in file " << file_path
              << std::endl;
    exit(1);
  }

  DenseGraph graph(num_vertices);
  std::unordered_map<int, bool>
      valid_vertices; // Para garantir que os vértices são válidos

  int u, v;
  while (file >> u >> v) {
    if (u < 0 || u >= num_vertices || v < 0 || v >= num_vertices) {
      std::cerr << "Warning: Ignoring invalid edge (" << u << ", " << v << ")"
                << std::endl;
      continue;
    }

    if (u != v) {
      graph.addEdge(u, v);
      valid_vertices[u] = true;
      valid_vertices[v] = true;
    }
  }

  // Verificar se o grafo contém pelo menos um nó válido
  if (valid_vertices.empty()) {
    std::cerr << "Error: The graph has no valid edges or nodes." << std::endl;
    exit(1);
  }

  return graph;
}

void execute_trial(size_t trial, DenseGraph &graph,
                   const AlgorithmParams &params,
                   std::vector<TrialResult> &results) {
  auto start_time = std::chrono::high_resolution_clock::now();

  std::vector<Chromosome> initial_population;
  for (size_t i = 0; i < static_cast<size_t>(graph.getVertexCount() /
                                             params.population_factor);
       ++i) {
    initial_population.push_back(h1(graph));
  }

  Population population(initial_population.size(), initial_population);
  SinglePoint crossover(params.crossover_rate);
  KTournamentSelection selector(params.tournament_size);

  Chromosome best_solution = population.get_best_chromosome();
  size_t stagnant_generations = 0;

  for (size_t generation = 0; generation < params.generations; ++generation) {
    population.evolve(selector, crossover);
    Chromosome new_best_solution = population.get_best_chromosome();

    if (new_best_solution.get_fitness() < best_solution.get_fitness()) {
      best_solution = new_best_solution;
      stagnant_generations = 0;
    } else {
      stagnant_generations++;
    }

    if (stagnant_generations >= params.max_stagnant) {
      break;
    }
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  uint64_t elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(
                              end_time - start_time)
                              .count();

  // Extrair apenas o nome do arquivo (removendo o caminho)
  std::string file_name =
      params.file_path.substr(params.file_path.find_last_of("/\\") + 1);

  // Remover a extensão ".txt" se existir
  size_t dot_pos = file_name.find_last_of(".");
  if (dot_pos != std::string::npos && file_name.substr(dot_pos) == ".txt") {
    file_name = file_name.substr(0, dot_pos);
  }

  results.push_back({file_name, static_cast<size_t>(graph.getVertexCount()),
                     static_cast<size_t>(graph.getEdgeCount()),
                     best_solution.get_fitness(), elapsed_time});
}

int main(int argc, char *argv[]) {
  AlgorithmParams params = parse_args(argc, argv);

  std::cout << "Loading graph from file: " << params.file_path << std::endl;
  DenseGraph graph = load_graph_from_file(params.file_path);

  if (graph.getVertexCount() == 0) {
    std::cerr << "Error: The graph has no nodes." << std::endl;
    return 1;
  }

  std::cout << "Graph loaded - Nodes: " << graph.getVertexCount()
            << ", Edges: " << graph.getEdgeCount() << std::endl;

  std::vector<TrialResult> results;
  results.reserve(params.trials);

  for (size_t trial = 0; trial < params.trials; ++trial) {
    std::cout << "Executing trial " << (trial + 1) << "..." << std::endl;
    execute_trial(trial, graph, params, results);
  }

  write_results_to_csv(results, params.output_file);

  std::cout << "Execution completed. Results saved to " << params.output_file
            << std::endl;
  return 0;
}
