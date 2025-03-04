#include "Chromosome.hpp"
#include "Graph.hpp"
#include "HeuristicGenerators.hpp"
#include "ListGraph.hpp"
#include "MatrixGraph.hpp"
#include "Population.hpp"
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <filesystem>
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
  bool matches_heuristic;
  std::string matched_heuristic;
  bool is_valid;
  bool density;
};

std::unique_ptr<Graph> createOptimalGraph(int numVertices, double density) {
  // eficiente Para grafos esparsos, ListGraph é melhor
  if (density > 0.5) {
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
  if (density > 0.5) {
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

bool chromosomes_equal(const Chromosome &c1, const Chromosome &c2) {
  if (c1.size() != c2.size()) {
    return false;
  }

  for (size_t i = 0; i < c1.size(); ++i) {
    if (c1.get_value(i) != c2.get_value(i)) {
      return false;
    }
  }

  return true;
}

void ensure_csv_header(const std::string &filename) {
  bool file_exists = std::filesystem::exists(filename);

  std::ofstream file;
  if (!file_exists) {
    file.open(filename);
    file << "graph_name,graph_order,graph_size,fitness_value,elapsed_time("
            "microsecond),matches_heuristic,heuristic_matched,is_valid,"
            "density\n";
  } else {
    file.open(filename, std::ios::app);
  }
  file.close();
}

void write_result_to_csv(const std::string &filename,
                         const TrialResult &result) {
  std::ofstream file(filename, std::ios::app);
  file << result.graph_name << "," << result.node_count << ","
       << result.edge_count << "," << result.fitness << ","
       << result.elapsed_micros << ","
       << (result.matches_heuristic ? "yes" : "no") << ","
       << result.matched_heuristic << "," << (result.is_valid ? "yes" : "no")
       << ","                     // Nova coluna: is_valid
       << result.density << "\n"; // Nova coluna: density
  file.close();
}

TrialResult run_genetic_algorithm(const std::unique_ptr<Graph> &graph,
                                  const AlgorithmParams &params) {
  HeuristicGenerators heuristicHandle;

  size_t pop_size =
      static_cast<size_t>(graph->order() * params.population_factor);

  auto heuristic_h1 = heuristicHandle.h1(*graph);
  auto heuristic_h2 = heuristicHandle.h2(*graph);
  auto heuristic_h3 = heuristicHandle.h3(*graph);
  // auto heuristic_h4 = heuristicHandle.h4(*graph);
  auto heuristic_h5 = heuristicHandle.h5(*graph);

  std::vector<Chromosome> initialChromosomes;
  initialChromosomes.push_back(heuristic_h1);
  initialChromosomes.push_back(heuristic_h2);
  initialChromosomes.push_back(heuristic_h3);
  // initialChromosomes.push_back(heuristic_h4);
  initialChromosomes.push_back(heuristic_h5);

  while (initialChromosomes.size() < pop_size) {
    initialChromosomes.push_back(heuristicHandle.h1(*graph));
  }

  auto start_time = std::chrono::high_resolution_clock::now();

  Population population(pop_size, initialChromosomes);

  size_t generation = 0;
  size_t stagnant_generations = 0;
  int best_fitness = std::numeric_limits<int>::max();
  Chromosome best_chromosome = initialChromosomes[0];

  for (const auto &chromosome : initialChromosomes) {
    if (chromosome.get_fitness() < best_fitness) {
      best_fitness = chromosome.get_fitness();
      best_chromosome = chromosome;
    }
  }

  while (generation < params.generations &&
         stagnant_generations < params.max_stagnant) {
    population.evolve(*graph, params.tournament_size, params.crossover_rate);

    const auto &chromosomes = population.getChromosomes();
    auto best_it =
        std::min_element(chromosomes.begin(), chromosomes.end(),
                         [](const Chromosome &a, const Chromosome &b) {
                           return a.get_fitness() < b.get_fitness();
                         });

    if (best_it->get_fitness() < best_fitness) {
      best_fitness = best_it->get_fitness();
      best_chromosome = *best_it;
      stagnant_generations = 0;

      std::cout << "Nova melhor solução na geração " << generation
                << ": fitness = " << best_fitness << std::endl;
    } else {
      stagnant_generations++;
    }

    generation++;
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(
      end_time - start_time);

  std::cout << "Algoritmo genético concluído após " << generation << " gerações"
            << std::endl;
  std::cout << "Melhor fitness: " << best_fitness << std::endl;

  bool matches_heuristic = false;
  std::string matched_heuristic = "none";

  if (chromosomes_equal(best_chromosome, heuristic_h1)) {
    matches_heuristic = true;
    matched_heuristic = "h1";
  } else if (chromosomes_equal(best_chromosome, heuristic_h2)) {
    matches_heuristic = true;
    matched_heuristic = "h2";
  } else if (chromosomes_equal(best_chromosome, heuristic_h3)) {
    matches_heuristic = true;
    matched_heuristic = "h3";
    // } else if (chromosomes_equal(best_chromosome, heuristic_h4)) {
    //   matches_heuristic = true;
    //   matched_heuristic = "h4";
  } else if (chromosomes_equal(best_chromosome, heuristic_h5)) {
    matches_heuristic = true;
    matched_heuristic = "h5";
  }

  // Verificar se a solução é válida
  bool is_valid = valid_totalrd(graph.get(), best_chromosome);
  if (!is_valid) {
    std::cout << "ATENÇÃO: A melhor solução encontrada não é válida!"
              << std::endl;
  } else {
    std::cout << "A solução encontrada é válida." << std::endl;
  }

  double density = static_cast<double>(graph->size()) /
                   (graph->order() * (graph->order() - 1) / 2.0);

  std::filesystem::path input_path(params.file_path);
  TrialResult result;
  result.graph_name = input_path.filename().string();
  result.node_count = graph->order();
  result.edge_count = graph->size();
  result.fitness = best_chromosome.get_fitness();
  result.elapsed_micros = elapsed.count();
  result.matches_heuristic = matches_heuristic;
  result.matched_heuristic = matched_heuristic;
  result.is_valid = is_valid;
  result.density = density;
  return result;
}

int main(int argc, char *argv[]) {
  auto params = parse_args(argc, argv);

  ensure_csv_header(params.output_file);

  auto graph = load_and_normalize_graph(params.file_path);

  for (size_t trial = 0; trial < params.trials; ++trial) {
    std::cout << "\n===== Executando trial " << (trial + 1) << " de "
              << params.trials << " =====\n"
              << std::endl;

    auto result = run_genetic_algorithm(graph, params);

    std::cout << "\nResultado do trial " << (trial + 1) << ":" << std::endl;
    std::cout << "  Nome do grafo: " << result.graph_name << std::endl;
    std::cout << "  Ordem (vértices): " << result.node_count << std::endl;
    std::cout << "  Tamanho (arestas): " << result.edge_count << std::endl;
    std::cout << "  Fitness: " << result.fitness << std::endl;
    std::cout << "  Tempo (microssegundos): " << result.elapsed_micros
              << std::endl;
    std::cout << "  Corresponde a heurística: "
              << (result.matches_heuristic ? "Sim" : "Não");
    if (result.matches_heuristic) {
      std::cout << " (" << result.matched_heuristic << ")";
    }
    std::cout << std::endl;

    // Escrever no CSV
    write_result_to_csv(params.output_file, result);
  }

  std::cout << "\nTodos os resultados foram salvos em: " << params.output_file
            << std::endl;

  return 0;
  // if (graph->order() == 0) {
  //   std::cerr << "Error: The graph has no nodes." << std::endl;
  //   return 1;
  // }
  //
  // // Test h1 heuristic
  // std::cout << "Testing h1 heuristic:\n";
  // auto solution_h1 = heuristicHandle.h1(*graph);
  // auto is_valid_h1 = valid_totalrd(graph.get(), solution_h1);
  // std::cout << "Solution h1 fitness: " << solution_h1.get_fitness()
  //           << std::endl;
  // std::cout << (is_valid_h1 ? "h1 solution is valid\n"
  //                           : "h1 solution is not valid\n");
  //
  // // Test h2 heuristic
  // std::cout << "\nTesting h2 heuristic:\n";
  // auto solution_h2 = heuristicHandle.h2(*graph);
  // auto is_valid_h2 = valid_totalrd(graph.get(), solution_h2);
  // std::cout << "Solution h2 fitness: " << solution_h2.get_fitness()
  //           << std::endl;
  // std::cout << (is_valid_h2 ? "h2 solution is valid\n"
  //                           : "h2 solution is not valid\n");
  //
  // // Test h3 heuristic
  // std::cout << "\nTesting h3 heuristic:\n";
  // auto solution_h3 = heuristicHandle.h3(*graph);
  // auto is_valid_h3 = valid_totalrd(graph.get(), solution_h3);
  // std::cout << "Solution h3 fitness: " << solution_h3.get_fitness()
  //           << std::endl;
  // std::cout << (is_valid_h3 ? "h3 solution is valid\n"
  //                           : "h3 solution is not valid\n");
  //
  // // Test h4 heuristic
  // std::cout << "\nTesting h4 heuristic:\n";
  // auto solution_h4 = heuristicHandle.h4(*graph);
  // auto is_valid_h4 = valid_totalrd(graph.get(), solution_h4);
  // std::cout << "Solution h4 fitness: " << solution_h4.get_fitness()
  //           << std::endl;
  // std::cout << (is_valid_h4 ? "h4 solution is valid\n"
  //                           : "h4 solution is not valid\n");
  //
  // // Test h5 heuristic
  // std::cout << "\nTesting h5 heuristic:\n";
  // auto solution_h5 = heuristicHandle.h5(*graph);
  // auto is_valid_h5 = valid_totalrd(graph.get(), solution_h5);
  // std::cout << "Solution h5 fitness: " << solution_h5.get_fitness()
  //           << std::endl;
  // std::cout << (is_valid_h5 ? "h5 solution is valid\n"
  //                           : "h5 solution is not valid\n");
  //
  // // Compare solutions
  // std::cout << "\nComparison:\n";
  // std::vector<std::pair<std::string, double>> heuristics = {
  //     {"h1", solution_h1.get_fitness()},
  //     {"h2", solution_h2.get_fitness()},
  //     {"h3", solution_h3.get_fitness()},
  //     {"h4", solution_h4.get_fitness()},
  //     {"h5", solution_h5.get_fitness()}};
  //
  // auto best_solution = *std::min_element(
  //     heuristics.begin(), heuristics.end(),
  //     [](const auto &a, const auto &b) { return a.second < b.second; });
  //
  // std::cout << best_solution.first << " found the best solution with fitness
  // "
  //           << best_solution.second << "\n";
  //
}
