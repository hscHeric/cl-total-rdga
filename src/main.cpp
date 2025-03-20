#include "Chromosome.hpp"
#include "HeuristicGenerators.hpp"
#include "ListGraph.hpp"
#include "MatrixGraph.hpp"
#include "Population.hpp"
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <thread>

#define IRACE 0
#define DEBUG ((~IRACE) & 0)
#define PARALLEL 1
#define NUM_THREADS 8

struct AlgorithmParams {
  size_t max_stagnant = 100;
  size_t generations = 341;
  size_t tournament_size = 3;
  double crossover_rate = 0.8933;
  double population_factor = 4;
  double elitism_rate = 0.1;
  double mutation_rate = 0.05;
  std::string file_path;
  size_t trials = 1;
  std::string output_file = "results.csv";
};

struct TrialResult {
  std::string graph_name;
  size_t node_count;
  size_t edge_count;
  float graph_density;
  int fitness;
  uint64_t elapsed_micros;
  bool matches_heuristic;
  std::string matched_heuristic;
  bool is_valid;
  bool is_dense;
};

struct HeuristicBuffer {
  Chromosome h2;
  Chromosome h3;
  Chromosome h4;
  Chromosome h5;
  bool initialized = false;
};

HeuristicBuffer global_heuristic_buffer_l;
HeuristicBuffer global_heuristic_buffer_m;

void init_heuristic_buffer_l(const ListGraph &graph) {
  if (!global_heuristic_buffer_l.initialized) {
    HeuristicGenerators heuristicHandle;
    global_heuristic_buffer_l.h2 = heuristicHandle.h2_l(graph);
    global_heuristic_buffer_l.h3 = heuristicHandle.h3_l(graph);
    global_heuristic_buffer_l.h4 = heuristicHandle.h4_l(graph);
    global_heuristic_buffer_l.h5 = heuristicHandle.h5_l(graph);
    global_heuristic_buffer_l.initialized = true;

#if (DEBUG)
    std::cout << "[DEBUG] Inicializando o HeuristicBuffer - ListGraph\n";
#endif
  }
}

void init_heuristic_buffer_m(const MatrixGraph &graph) {
  if (!global_heuristic_buffer_m.initialized) {
    HeuristicGenerators heuristicHandle;
    global_heuristic_buffer_m.h2 = heuristicHandle.h2_m(graph);
    global_heuristic_buffer_m.h3 = heuristicHandle.h3_m(graph);
    global_heuristic_buffer_m.h4 = heuristicHandle.h4_m(graph);
    global_heuristic_buffer_m.h5 = heuristicHandle.h5_m(graph);
    global_heuristic_buffer_m.initialized = true;

#if (DEBUG)
    std::cout << "[DEBUG] Inicializando o HeuristicBuffer - MatrixGraph\n";
#endif
  }
}

void load_and_normalize_graph(const std::string &filename,
                              bool &graph_is_matrix, ListGraph &graph_l,
                              MatrixGraph &graph_m);
AlgorithmParams parse_args(int argc, char *argv[]);
bool valid_totalrd_l(const ListGraph &graph, const Chromosome &chromosome);
bool valid_totalrd_m(const MatrixGraph &graph, const Chromosome &chromosome);
bool chromosomes_equal(const Chromosome &c1, const Chromosome &c2);
void ensure_csv_header(const std::string &filename);
void write_result_to_csv(const std::string &filename,
                         const TrialResult &result);
TrialResult run_genetic_algorithm_l(const ListGraph &graph,
                                    const AlgorithmParams &params);
TrialResult run_genetic_algorithm_m(const MatrixGraph &graph,
                                    const AlgorithmParams &params);

/**
 * @brief Main
 */
int main(int argc, char *argv[]) {
  auto params = parse_args(argc, argv);
#if !IRACE
  ensure_csv_header(params.output_file);
#endif
  // Criar os grafos e determinar qual será utilizado
  ListGraph graph_l(1);
  MatrixGraph graph_m(1);
  bool graph_is_matrix = false;
  load_and_normalize_graph(params.file_path, graph_is_matrix, graph_l, graph_m);
  if (graph_is_matrix) {
    if (!graph_m.get_isolated_vertices().empty()) {
#if DEBUG
      std::cout << "A execução foi interrompida pois o grafo contém vértices "
                   "isolados."
                << std::endl;
#endif
      return 1;
    }
    init_heuristic_buffer_m(graph_m);
  } else {
    if (!graph_l.get_isolated_vertices().empty()) {
#if DEBUG
      std::cout << "A execução foi interrompida pois o grafo contém vértices "
                   "isolados."
                << std::endl;
#endif
      return 1;
    }
    init_heuristic_buffer_l(graph_l);
  }
  // Modo IRACE: Executa apenas uma tentativa e imprime apenas o fitness
#if IRACE
  TrialResult result;
  if (graph_is_matrix) {
    result = run_genetic_algorithm_m(graph_m, params);
  } else {
    result = run_genetic_algorithm_l(graph_l, params);
  }
  std::cout << result.fitness << std::endl;
  return 0;
#else
#if PARALLEL
  std::vector<TrialResult> results(params.trials);
  std::vector<std::thread> threads;
  std::mutex csv_mutex; // Mutex para proteger a escrita no arquivo CSV

  auto thread_func = [&](size_t trial_id) {
    TrialResult result;
    if (graph_is_matrix) {
      result = run_genetic_algorithm_m(graph_m, params);
    } else {
      result = run_genetic_algorithm_l(graph_l, params);
    }

#if DEBUG
    std::cout << "\nResultado do trial " << (trial_id + 1) << ":" << std::endl;
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
#endif

    std::lock_guard<std::mutex> lock(csv_mutex);
    write_result_to_csv(params.output_file, result);
  };

  // Create and start threads
  for (size_t trial = 0; trial < params.trials; ++trial) {
    if (threads.size() < NUM_THREADS) {
      threads.emplace_back(thread_func, trial);
    } else {
      threads[trial % NUM_THREADS].join();
      threads[trial % NUM_THREADS] = std::thread(thread_func, trial);
    }
  }

  // Wait for all remaining threads to finish
  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }

#else
  for (size_t trial = 0; trial < params.trials; ++trial) {
#if DEBUG
    std::cout << "\n===== Executando trial " << (trial + 1) << " de "
              << params.trials << " =====\n"
              << std::endl;
#endif
    TrialResult result;
    if (graph_is_matrix) {
      result = run_genetic_algorithm_m(graph_m, params);
    } else {
      result = run_genetic_algorithm_l(graph_l, params);
    }
#if DEBUG
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
#endif
    write_result_to_csv(params.output_file, result);
  }
#endif // PARALLEL

#if DEBUG
  std::cout << "\nTodos os resultados foram salvos em: " << params.output_file
            << std::endl;
  return 0;
#endif // DEBUG
#endif // !IRACE
}

// ----------------------------------------------------------------------------------

void load_and_normalize_graph(const std::string &filename,
                              bool &graph_is_matrix, ListGraph &graph_l,
                              MatrixGraph &graph_m) {

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
  /*
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
  */

  // std::cout << "  Mapeamento aplicado: ";
  // for (const auto &[orig, norm] : id_map) {
  //   std::cout << orig << "->" << norm << " ";
  // }
  /*
  std::cout << std::endl;
  std::cout << "  Vértices normalizados: " << num_vertices << std::endl;
  std::cout << "  Arestas únicas: " << unique_edges.size() << std::endl;
  std::cout << "  Densidade: " << density << std::endl;
  std::cout << "  Implementação selecionada: "
            << (density > 0.3 ? "MatrixGraph" : "ListGraph") << std::endl;
            */

  // Criar grafo com implementação adequada
  // std::unique_ptr<Graph> graph;

  if (density > 0.5) {
    graph_is_matrix = true;
    graph_m = MatrixGraph(num_vertices);

    // Adicionar arestas normalizadas
    for (const auto &edge : unique_edges) {
      graph_m.add_edge(edge.first, edge.second);
    }
  } else {
    graph_is_matrix = false;
    graph_l = ListGraph(num_vertices);

    // Adicionar arestas normalizadas
    for (const auto &edge : unique_edges) {
      graph_l.add_edge(edge.first, edge.second);
    }
  }
}

AlgorithmParams parse_args(int argc, char *argv[]) {
  AlgorithmParams params;

#if IRACE
  if (argc < 5) {
    std::cerr << "Usage: " << argv[0] << " [options] <graph_file>\n";
    exit(1);
  }
  params.file_path = argv[4];

#else
  if (argc < 2) { // Modo normal: arquivo de entrada em argv[1]
    std::cerr << "Usage: " << argv[0] << " <graph_file> [options]\n"
              << "Options:\n"
              << "  --crossover VALUE\n"
              << "  --stagnation VALUE\n"
              << "  --generations VALUE\n"
              << "  --population VALUE\n"
              << "  --tournament VALUE\n"
              << "  --elitism VALUE\n"
              << "  --mutation VALUE\n"
              << "  --trials VALUE\n"
              << "  --output FILE\n";
    exit(1);
  }
  params.file_path =
      argv[1]; // No modo normal, o caminho do arquivo está em argv[1]
#endif

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
    } else if (arg == "--elitism" && i + 1 < argc) {
      params.elitism_rate = std::stod(argv[++i]);
    } else if (arg == "--mutation" && i + 1 < argc) {
      params.mutation_rate = std::stod(argv[++i]);
#if !IRACE
    } else if (arg == "--trials" && i + 1 < argc) {
      params.trials = std::stoul(argv[++i]);
    } else if (arg == "--output" && i + 1 < argc) {
      params.output_file = argv[++i];
#endif
#if IRACE
    } else {
      continue;
#endif
    }
  }

#if IRACE
  params.trials = 1;
#endif

  return params;
}

void ensure_csv_header(const std::string &filename) {
  bool file_exists = std::filesystem::exists(filename);

  std::ofstream file;
  if (!file_exists) {
    file.open(filename);
    file << "graph_name,graph_order,graph_size,fitness_value,elapsed_time("
            "microsecond),matches_heuristic,heuristic_matched,is_valid,"
            "density,is_dense\n";
  } else {
    file.open(filename, std::ios::app);
  }
  file.close();
}

bool valid_totalrd_m(const MatrixGraph &graph, const Chromosome &chromosome) {
  for (int v : graph.get_vertices()) {
    int label = chromosome.get_value(v);
    bool hasValidNeighbor = false;

    auto &neighbors = graph.get_neighbors(v);
    for (int i = 0; i < static_cast<int>(neighbors.size()); ++i) {
      if (neighbors[i] == 1) {
        int neighbor_label = chromosome.get_value(i);
        if ((label == 0 && neighbor_label == 2) ||
            (label > 0 && neighbor_label > 0)) {
          hasValidNeighbor = true;
        }
      }
    }

    if (!hasValidNeighbor) {
#if DEBUG
      std::cout << "Solução inválida: vértice " << v << " com f(v) = " << label
                << " não possui vizinho válido." << std::endl;
#endif
      return false;
    }
  }
  return true;
}

bool valid_totalrd_l(const ListGraph &graph, const Chromosome &chromosome) {
  for (int v : graph.get_vertices()) {
    int label = chromosome.get_value(v);
    bool hasValidNeighbor = false;

    for (const int &neighbor : graph.get_neighbors(v)) {
      int neighbor_label = chromosome.get_value(neighbor);
      if ((label == 0 && neighbor_label == 2) ||
          (label > 0 && neighbor_label > 0)) {
        hasValidNeighbor = true;
      }
    }

    if (!hasValidNeighbor) {
#if DEBUG
      std::cout << "Solução inválida: vértice " << v << " com f(v) = " << label
                << " não possui vizinho válido." << std::endl;
#endif
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

void write_result_to_csv(const std::string &filename,
                         const TrialResult &result) {
  std::ofstream file(filename, std::ios::app);
  file << result.graph_name << "," << result.node_count << ","
       << result.edge_count << "," << result.fitness << ","
       << result.elapsed_micros << ","
       << (result.matches_heuristic ? "yes" : "no") << ","
       << result.matched_heuristic << "," << (result.is_valid ? "yes" : "no")
       << "," // Nova coluna: is_valid
       << result.graph_density << "," << std::boolalpha << result.is_dense
       << "\n"; // Nova coluna: is_dense
  file.close();
}

TrialResult run_genetic_algorithm_m(const MatrixGraph &graph,
                                    const AlgorithmParams &params) {
  HeuristicGenerators heuristicHandle;

  // o tamanho da população é uma fração do tamanho do grafo
  size_t pop_size =
      static_cast<size_t>(graph.order() / params.population_factor);

  auto start_time = std::chrono::high_resolution_clock::now();
  auto heuristic_h1 = heuristicHandle.h1_m(graph);
#if DEBUG
  auto heuristic_h2 = heuristicHandle.h2_m(graph);
  auto heuristic_h3 = heuristicHandle.h3_m(graph);
  auto heuristic_h4 = heuristicHandle.h4_m(graph);
  auto heuristic_h5 = heuristicHandle.h5_m(graph);

  std::cout << "========================================\n";
  std::cout << "          IMPLEMENTAÇÃO: MatrixGraph    \n";
  std::cout << "========================================\n";

  // Lista de heurísticas e seus resultados
  std::vector<std::pair<std::string, decltype(heuristic_h1)>> heuristics = {
      {"h1", heuristic_h1},
      {"h2", heuristic_h2},
      {"h3", heuristic_h3},
      {"h4", heuristic_h4},
      {"h5", heuristic_h5}};

  // Itera sobre as heurísticas e valida cada uma
  for (const auto &[name, heuristic] : heuristics) {
    std::cout << "[INFO] Testando solução com " << name << "...\n";

    if (valid_totalrd_m(graph, heuristic)) {
      std::cout << "[OK] Solução utilizando " << name << " é VÁLIDA\n";
    } else {
      std::cout << "[ERRO - MatrixGraph] Solução utilizando " << name
                << " é INVÁLIDA\n";
    }

    std::cout << "----------------------------------------\n";
  }
#endif

  std::vector<Chromosome> initialChromosomes;

  int s_h1 = 0.6 * pop_size;
  int s_h2 = 0.1 * pop_size;
  int s_h3 = 0.1 * pop_size;
  int s_h4 = 0.1 * pop_size;

  int count = 0;
  while (initialChromosomes.size() < pop_size) {
    if (count < s_h1) {
      initialChromosomes.push_back(heuristicHandle.h1_m(graph));
    } else if (count < s_h1 + s_h2) {
      initialChromosomes.push_back(global_heuristic_buffer_m.h2);
    } else if (count < s_h1 + s_h2 + s_h3) {
      initialChromosomes.push_back(global_heuristic_buffer_m.h3);
    } else if (count < s_h1 + s_h2 + s_h3 + s_h4) {
      initialChromosomes.push_back(global_heuristic_buffer_m.h4);
    } else {
      initialChromosomes.push_back(global_heuristic_buffer_m.h5);
    }
    count++;
  }

  while (initialChromosomes.size() < pop_size) {
    initialChromosomes.push_back(heuristicHandle.h1_m(graph));
  }

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

    population.evolve_m(graph, params.tournament_size, params.crossover_rate,
                        params.elitism_rate, params.mutation_rate);

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
#if DEBUG
      std::cout << "Nova melhor solução na geração " << generation
                << ": fitness = " << best_fitness << std::endl;
#endif
    } else {
      stagnant_generations++;
    }

    generation++;
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(
      end_time - start_time);

#if DEBUG
  std::cout << "Algoritmo genético concluído após " << generation << " gerações"
            << std::endl;
  std::cout << "Melhor fitness: " << best_fitness << std::endl;
#endif

  bool matches_heuristic = false;
  std::string matched_heuristic = "none";

  if (chromosomes_equal(best_chromosome, heuristic_h1)) {
    matches_heuristic = true;
    matched_heuristic = "h1";
  } else if (chromosomes_equal(best_chromosome, global_heuristic_buffer_m.h2)) {
    matches_heuristic = true;
    matched_heuristic = "h2";
  } else if (chromosomes_equal(best_chromosome, global_heuristic_buffer_m.h3)) {
    matches_heuristic = true;
    matched_heuristic = "h3";
  } else if (chromosomes_equal(best_chromosome, global_heuristic_buffer_m.h4)) {
    matches_heuristic = true;
    matched_heuristic = "h4";
  } else if (chromosomes_equal(best_chromosome, global_heuristic_buffer_m.h5)) {
    matches_heuristic = true;
    matched_heuristic = "h5";
  }

  // Verificar se a solução é válida
  bool is_valid = valid_totalrd_m(graph, best_chromosome);
  if (!is_valid) {
#if DEBUG
    std::cout << "ATENÇÃO: A melhor solução encontrada não é válida!"
              << std::endl;
#endif
  } else {
#if DEBUG
    std::cout << "A solução encontrada é válida." << std::endl;
#endif
  }

  double density = static_cast<double>(graph.size()) /
                   (graph.order() * (graph.order() - 1) / 2.0);

  std::filesystem::path input_path(params.file_path);
  TrialResult result;
  result.graph_name = input_path.filename().string();
  result.node_count = graph.order();
  result.edge_count = graph.size();
  result.fitness = best_chromosome.get_fitness();
  result.elapsed_micros = elapsed.count();
  result.matches_heuristic = matches_heuristic;
  result.matched_heuristic = matched_heuristic;
  result.is_valid = is_valid;
  result.graph_density = density;
  result.is_dense = (density > 0.5) ? true : false;
  return result;
}

TrialResult run_genetic_algorithm_l(const ListGraph &graph,
                                    const AlgorithmParams &params) {
  HeuristicGenerators heuristicHandle;

  // o tamanho da população é uma fração do tamanho do grafo
  size_t pop_size =
      static_cast<size_t>(graph.order() / params.population_factor);

  auto start_time = std::chrono::high_resolution_clock::now();
  auto heuristic_h1 = heuristicHandle.h1_l(graph);

#if DEBUG
  auto heuristic_h2 = heuristicHandle.h2_l(graph);
  auto heuristic_h3 = heuristicHandle.h3_l(graph);
  auto heuristic_h4 = heuristicHandle.h4_l(graph);
  auto heuristic_h5 = heuristicHandle.h5_l(graph);

  std::cout << "========================================\n";
  std::cout << "          IMPLEMENTAÇÃO: ListGraph       \n";
  std::cout << "========================================\n";

  std::vector<std::pair<std::string, decltype(heuristic_h1)>> heuristics = {
      {"h1", heuristic_h1},
      {"h2", heuristic_h2},
      {"h3", heuristic_h3},
      {"h4", heuristic_h4},
      {"h5", heuristic_h5}};

  for (const auto &[name, heuristic] : heuristics) {
    std::cout << "[INFO] Testando solução com " << name << "...\n";

    if (valid_totalrd_l(graph, heuristic)) {
      std::cout << "[OK] Solução utilizando " << name << " é VÁLIDA\n";
    } else {
      std::cout << "[ERRO - ListGraph] Solução utilizando " << name
                << " é INVÁLIDA\n";
    }

    std::cout << "----------------------------------------\n";
  }
#endif

  std::vector<Chromosome> initialChromosomes;

  int s_h1 = 0.6 * pop_size;
  int s_h2 = 0.1 * pop_size;
  int s_h3 = 0.1 * pop_size;
  int s_h4 = 0.1 * pop_size;

  int count = 0;
  while (initialChromosomes.size() < pop_size) {
    if (count < s_h1) {
      initialChromosomes.push_back(heuristicHandle.h1_l(graph));
    } else if (count < s_h1 + s_h2) {
      initialChromosomes.push_back(global_heuristic_buffer_l.h2);
    } else if (count < s_h1 + s_h2 + s_h3) {
      initialChromosomes.push_back(global_heuristic_buffer_l.h3);
    } else if (count < s_h1 + s_h2 + s_h3 + s_h4) {
      initialChromosomes.push_back(global_heuristic_buffer_l.h4);
    } else {
      initialChromosomes.push_back(global_heuristic_buffer_l.h5);
    }
    count++;
  }

  // initialChromosomes.push_back(heuristic_h2);
  // initialChromosomes.push_back(heuristic_h3);
  // initialChromosomes.push_back(heuristic_h4);
  // initialChromosomes.push_back(heuristic_h5);

  // while (initialChromosomes.size() < pop_size) {
  //   initialChromosomes.push_back(heuristicHandle.h1_l(graph));
  // }

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

    population.evolve_l(graph, params.tournament_size, params.crossover_rate,
                        params.elitism_rate, params.mutation_rate);

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

#if DEBUG
      std::cout << "Nova melhor solução na geração " << generation
                << ": fitness = " << best_fitness << std::endl;
#endif
    } else {
      stagnant_generations++;
    }

    generation++;
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(
      end_time - start_time);

#if DEBUG
  std::cout << "Algoritmo genético concluído após " << generation << " gerações"
            << std::endl;
  std::cout << "Melhor fitness: " << best_fitness << std::endl;
#endif

  bool matches_heuristic = false;
  std::string matched_heuristic = "none";

  if (chromosomes_equal(best_chromosome, heuristic_h1)) {
    matches_heuristic = true;
    matched_heuristic = "h1";
  } else if (chromosomes_equal(best_chromosome, global_heuristic_buffer_l.h2)) {
    matches_heuristic = true;
    matched_heuristic = "h2";
  } else if (chromosomes_equal(best_chromosome, global_heuristic_buffer_l.h3)) {
    matches_heuristic = true;
    matched_heuristic = "h3";
  } else if (chromosomes_equal(best_chromosome, global_heuristic_buffer_l.h4)) {
    matches_heuristic = true;
    matched_heuristic = "h4";
  } else if (chromosomes_equal(best_chromosome, global_heuristic_buffer_l.h5)) {
    matches_heuristic = true;
    matched_heuristic = "h5";
  }

  // Verificar se a solução é válida
  bool is_valid = valid_totalrd_l(graph, best_chromosome);
  if (!is_valid) {
#if DEBUG
    std::cout << "ATENÇÃO: A melhor solução encontrada não é válida!"
              << std::endl;
#endif
  } else {
#if DEBUG
    std::cout << "A solução encontrada é válida." << std::endl;
#endif
  }

  double density = static_cast<double>(graph.size()) /
                   (graph.order() * (graph.order() - 1) / 2.0);

  std::filesystem::path input_path(params.file_path);
  TrialResult result;
  result.graph_name = input_path.filename().string();
  result.node_count = graph.order();
  result.edge_count = graph.size();
  result.fitness = best_chromosome.get_fitness();
  result.elapsed_micros = elapsed.count();
  result.matches_heuristic = matches_heuristic;
  result.matched_heuristic = matched_heuristic;
  result.is_valid = is_valid;
  result.graph_density = density;
  result.is_dense = (density > 0.5) ? true : false;
  return result;
}
