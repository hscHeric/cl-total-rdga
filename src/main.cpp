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
#include <set>
#include <string>

#define DEBUG 1

struct AlgorithmParams {
  size_t max_stagnant = 200;
  size_t generations = 800;
  size_t tournament_size = 5;
  double crossover_rate = 0.8;
  double population_factor = 3;
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
void init_heuristic_buffer_l(const ListGraph &graph);
void init_heuristic_buffer_m(const MatrixGraph &graph);

/**
 * @brief Main
 */
int main(int argc, char *argv[]) {
  auto params = parse_args(argc, argv);

  ensure_csv_header(params.output_file);

  // Cria dois grafos, e depois decide qual deles será inicializado,
  // dependendo da densidade do grafo do arquivo lido
  ListGraph graph_l(1);
  MatrixGraph graph_m(1);

  // Flag que indica se o grafo lido é uma matriz ou uma lista de adjacencias
  bool graph_is_matrix = false;

  load_and_normalize_graph(params.file_path, graph_is_matrix, graph_l, graph_m);

  if (graph_is_matrix) {
    init_heuristic_buffer_m(graph_m);
  } else {

    init_heuristic_buffer_l(graph_l);
  }

  for (size_t trial = 0; trial < params.trials; ++trial) {
    // std::cout << "\n===== Executando trial " << (trial + 1) << " de "
    //         << params.trials << " =====\n"
    //       << std::endl;

    TrialResult result;

    if (graph_is_matrix) {
      result = run_genetic_algorithm_m(graph_m, params);
    } else {
      result = run_genetic_algorithm_l(graph_l, params);
    }

    // std::cout << "\nResultado do trial " << (trial + 1) << ":" << std::endl;
    // std::cout << "  Nome do grafo: " << result.graph_name << std::endl;
    // std::cout << "  Ordem (vértices): " << result.node_count << std::endl;
    // std::cout << "  Tamanho (arestas): " << result.edge_count << std::endl;
    // std::cout << "  Fitness: " << result.fitness << std::endl;
    // std::cout << "  Tempo (microssegundos): " << result.elapsed_micros
    //          << std::endl;
    // std::cout << "  Corresponde a heurística: "
    //          << (result.matches_heuristic ? "Sim" : "Não");
    // if (result.matches_heuristic) {
    // std::cout << " (" << result.matched_heuristic << ")";
    // }
    // std::cout << std::endl;

    // Escrever no CSV
    write_result_to_csv(params.output_file, result);
  }

  std::cout << "\nTodos os resultados foram salvos em: " << params.output_file
            << std::endl;

  return 0;
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
      // std::cout << "Solução inválida: vértice " << v << " com f(v) = " <<
      // label
      //           << " não possui vizinho válido." << std::endl;
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
      // std::cout << "Solução inválida: vértice " << v << " com f(v) = " <<
      // label
      //           << " não possui vizinho válido." << std::endl;
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

  size_t pop_size =
      std::max(static_cast<size_t>(graph.order() / params.population_factor),
               static_cast<size_t>(5)); // Tamanho minimo da população como 5

  auto start_time = std::chrono::high_resolution_clock::now();
#if DEBUG
  auto heuristic_h1 = heuristicHandle.h1_m(graph);
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

  // Verificando se as heuristicas são deterministicas

  for (const auto &[name, heuristic] : heuristics) {
    std::vector<Chromosome> solutions;

    for (int i = 0; i < 10; ++i) {
      solutions.push_back(heuristic);
    }

    bool all_equal = true;
    for (size_t i = 1; i < solutions.size(); ++i) {
      if (!chromosomes_equal(solutions[0], solutions[i])) {
        all_equal = false;
        break;
      }
    }

    if (all_equal) {
      std::cout << "[DEBUG] Todas as 10 soluções geradas por " << name
                << " são idênticas. Pode haver falta de variação.\n";
    }
  }
#endif

  std::vector<Chromosome> initialChromosomes;
  size_t h1_count = static_cast<size_t>(pop_size * 0.6);
  size_t h2_count = static_cast<size_t>(pop_size * 0.1);
  size_t h3_count = static_cast<size_t>(pop_size * 0.1);
  size_t h4_count = static_cast<size_t>(pop_size * 0.1);
  size_t h5_count = static_cast<size_t>(pop_size * 0.1);

  if (h1_count + h2_count + h3_count + h4_count + h5_count < pop_size) {
    h1_count +=
        pop_size - (h1_count + h2_count + h3_count + h4_count + h5_count);
  }

  for (size_t i = 0; i < h1_count; ++i) {
    initialChromosomes.push_back(heuristicHandle.h1_m(graph));
  }

  for (size_t i = 0; i < h2_count; ++i) {
    initialChromosomes.push_back(global_heuristic_buffer_m.h2);
  }

  for (size_t i = 0; i < h3_count; ++i) {
    initialChromosomes.push_back(global_heuristic_buffer_m.h3);
  }

  for (size_t i = 0; i < h4_count; ++i) {
    initialChromosomes.push_back(global_heuristic_buffer_m.h4);
  }

  for (size_t i = 0; i < h5_count; ++i) {
    initialChromosomes.push_back(global_heuristic_buffer_m.h5);
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

    // std::cout << "-- geração " << generation << std::endl;
    population.evolve_m(graph, params.tournament_size, params.crossover_rate);

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

      // std::cout << "Nova melhor solução na geração " << generation
      //           << ": fitness = " << best_fitness << std::endl;
    } else {
      stagnant_generations++;
    }

    generation++;
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(
      end_time - start_time);

  // std::cout << "Algoritmo genético concluído após " << generation << "
  // gerações"
  //           << std::endl;
  // std::cout << "Melhor fitness: " << best_fitness << std::endl;
  //
  bool matches_heuristic = false;
  std::string matched_heuristic = "none";

  if (chromosomes_equal(best_chromosome, heuristic_h1)) {
    matches_heuristic = true;
    matched_heuristic = "h1";
  }
  /*else if (chromosomes_equal(best_chromosome, heuristic_h2)) {
    matches_heuristic = true;
    matched_heuristic = "h2";
  } else if (chromosomes_equal(best_chromosome, heuristic_h3)) {
    matches_heuristic = true;
    matched_heuristic = "h3";
  } else if (chromosomes_equal(best_chromosome, heuristic_h4)) {
    matches_heuristic = true;
    matched_heuristic = "h4";
  }*/
  else if (chromosomes_equal(best_chromosome, heuristic_h5)) {
    matches_heuristic = true;
    matched_heuristic = "h5";
  }

  // Verificar se a solução é válida
  bool is_valid = valid_totalrd_m(graph, best_chromosome);
  if (!is_valid) {
    // std::cout << "ATENÇÃO: A melhor solução encontrada não é válida!"
    //           << std::endl;
  } else {
    // std::cout << "A solução encontrada é válida." << std::endl;
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

  size_t pop_size =
      std::max(static_cast<size_t>(graph.order() / params.population_factor),
               static_cast<size_t>(5)); // Tamanho minimo da população como 5

  auto start_time = std::chrono::high_resolution_clock::now();

#if DEBUG
  auto heuristic_h1 = heuristicHandle.h1_l(graph);
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

  for (const auto &[name, heuristic] : heuristics) {
    std::vector<Chromosome> solutions;

    for (int i = 0; i < 10; ++i) {
      solutions.push_back(heuristic);
    }

    bool all_equal = true;
    for (size_t i = 1; i < solutions.size(); ++i) {
      if (!chromosomes_equal(solutions[0], solutions[i])) {
        all_equal = false;
        break;
      }
    }

    if (all_equal) {
      std::cout << "[DEBUG] Todas as 10 soluções geradas por " << name
                << " são idênticas. Pode haver falta de variação.\n";
    } else {
      std::cout << "[DEBUG] Existe soluções diferentes geradas por " << name
                << "\n";
    }
  }
#endif

  std::vector<Chromosome> initialChromosomes;

  size_t h1_count = static_cast<size_t>(pop_size * 0.6);
  size_t h2_count = static_cast<size_t>(pop_size * 0.1);
  size_t h3_count = static_cast<size_t>(pop_size * 0.1);
  size_t h4_count = static_cast<size_t>(pop_size * 0.1);
  size_t h5_count = static_cast<size_t>(pop_size * 0.1);

  if (h1_count + h2_count + h3_count + h4_count + h5_count < pop_size) {
    h1_count +=
        pop_size - (h1_count + h2_count + h3_count + h4_count + h5_count);
  }

  for (size_t i = 0; i < h1_count; ++i) {
    initialChromosomes.push_back(heuristicHandle.h1_l(graph));
  }

  for (size_t i = 0; i < h2_count; ++i) {
    initialChromosomes.push_back(global_heuristic_buffer_l.h2);
  }

  for (size_t i = 0; i < h3_count; ++i) {
    initialChromosomes.push_back(global_heuristic_buffer_l.h3);
  }

  for (size_t i = 0; i < h4_count; ++i) {
    initialChromosomes.push_back(global_heuristic_buffer_l.h4);
  }

  for (size_t i = 0; i < h5_count; ++i) {
    initialChromosomes.push_back(global_heuristic_buffer_l.h5);
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
    // std::cout << "-- geração " << generation << std::endl;
    population.evolve_l(graph, params.tournament_size, params.crossover_rate);

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

      // std::cout << "Nova melhor solução na geração " << generation
      //         << ": fitness = " << best_fitness << std::endl;
    } else {
      stagnant_generations++;
    }

    generation++;
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(
      end_time - start_time);

  // std::cout << "Algoritmo genético concluído após " << generation << "
  // gerações"
  //         << std::endl;
  // std::cout << "Melhor fitness: " << best_fitness << std::endl;

  bool matches_heuristic = false;
  std::string matched_heuristic = "none";

  if (chromosomes_equal(best_chromosome, heuristic_h1)) {
    matches_heuristic = true;
    matched_heuristic = "h1";
  }
  /*else if (chromosomes_equal(best_chromosome, heuristic_h2)) {
    matches_heuristic = true;
    matched_heuristic = "h2";
  } else if (chromosomes_equal(best_chromosome, heuristic_h3)) {
    matches_heuristic = true;
    matched_heuristic = "h3";
  } else if (chromosomes_equal(best_chromosome, heuristic_h4)) {
    matches_heuristic = true;
    matched_heuristic = "h4";
  }*/
  else if (chromosomes_equal(best_chromosome, heuristic_h5)) {
    matches_heuristic = true;
    matched_heuristic = "h5";
  }

  // Verificar se a solução é válida
  bool is_valid = valid_totalrd_l(graph, best_chromosome);
  if (!is_valid) {
    // std::cout << "ATENÇÃO: A melhor solução encontrada não é válida!"
    //         << std::endl;
  } else {
    // std::cout << "A solução encontrada é válida." << std::endl;
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
