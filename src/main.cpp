#include "Chromosome.hpp"
#include "Graph.hpp"
#include "HeuristicGenerators.hpp"
#include "ListGraph.hpp"
#include "MatrixGraph.hpp"
#include <cstddef>
#include <cstdint>
#include <fstream>
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

Graph *load_graph(const std::string &filename, double density_threshold = 0.5) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Não foi possível abrir o arquivo: " + filename);
  }

  int num_vertices, num_edges;
  file >> num_vertices >> num_edges;

  // Primeira passagem para encontrar o maior índice de vértice
  int max_vertex_index = -1;
  std::vector<std::pair<int, int>> edges;

  int u, v;
  file.clear();
  file.seekg(0);                     // Volta ao início do arquivo
  file >> num_vertices >> num_edges; // Pula o cabeçalho novamente

  while (file >> u >> v) {
    max_vertex_index = std::max(max_vertex_index, std::max(u, v));
    edges.push_back({u, v});
  }

  int actual_vertices = max_vertex_index + 1;

  double max_edges = (actual_vertices * (actual_vertices - 1)) / 2.0;
  double density =
      (max_edges > 0) ? (static_cast<double>(num_edges) / max_edges) : 0.0;

  Graph *graph;
  if (density < density_threshold) {
    graph = new ListGraph(actual_vertices);
  } else {
    graph = new MatrixGraph(actual_vertices);
  }

  // Adiciona as arestas
  std::unordered_set<std::string> edge_set;
  for (const auto &edge : edges) {
    u = edge.first;
    v = edge.second;

    if (u == v)
      continue;

    std::string edgeKey = (u < v) ? std::to_string(u) + "_" + std::to_string(v)
                                  : std::to_string(v) + "_" + std::to_string(u);

    if (edge_set.find(edgeKey) == edge_set.end()) {
      try {
        graph->add_edge(u, v);
        edge_set.insert(edgeKey);
      } catch (const std::exception &e) {
        std::cerr << "Aviso: " << e.what() << " para aresta (" << u << ", " << v
                  << ")\n";
      }
    }
  }

  file.close();
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
  bool isValid = true;
  graph->for_each_vertex([&](int v) {
    int label = chromosome.get_value(v);
    bool hasValidNeighbor = false;

    graph->for_each_neighbor(v, [&](int neighbor) {
      int neighborLabel = chromosome.get_value(neighbor);

      if (label == 0 && neighborLabel == 2) {
        hasValidNeighbor = true;
      }
      if ((label == 1 || label == 2) &&
          (neighborLabel == 1 || neighborLabel == 2)) {
        hasValidNeighbor = true;
      }
    });

    if ((label == 0 && !hasValidNeighbor) ||
        ((label == 1 || label == 2) && !hasValidNeighbor)) {
      // std::cerr << "Erro: vértice " << v << " (rótulo " << label
      //           << ") não atende aos critérios de dominância.\n";
      isValid = false;
    }
  });

  return isValid;
}

int main(int argc, char *argv[]) {
  auto params = parse_args(argc, argv);
  auto graph = load_graph(params.file_path);
  HeuristicGenerators heuristicHandle;

  if (graph->order() == 0) {
    std::cerr << "Error: The graph has no nodes." << std::endl;
    delete graph;
    return 1;
  }

  // Gerar uma solução usando h1
  Chromosome solution = heuristicHandle.h1(*graph);

  // Verificar se a solução é válida
  auto a = valid_totalrd(graph, solution);
  std::cout << a << "\n";
  delete graph;
  return 0;
}
