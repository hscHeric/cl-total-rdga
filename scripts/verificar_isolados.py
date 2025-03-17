import os
import sys
import networkx as nx


def carregar_grafos(base_dir):
    """Percorre a pasta e verifica grafos com vértices isolados."""
    grafos_com_isolados = []

    # Percorre todas as pastas dentro do diretório fornecido
    for root, _, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".txt"):  # Considera apenas arquivos .txt
                file_path = os.path.join(root, file)

                # Criando o grafo não-direcionado
                G = nx.Graph()

                # Lendo as arestas do arquivo
                with open(file_path, "r") as f:
                    for line in f:
                        u, v = map(int, line.split())
                        G.add_edge(u, v)

                # Verifica se há vértices isolados
                isolated_nodes = list(nx.isolates(G))
                if isolated_nodes:
                    grafos_com_isolados.append((file, isolated_nodes))

    return grafos_com_isolados


if __name__ == "__main__":
    # Verifica se o usuário forneceu o caminho da pasta
    if len(sys.argv) != 2:
        print("Uso: python verificar_isolados.py <caminho_dos_grafos>")
        sys.exit(1)

    # Obtém o diretório da linha de comando
    caminho_dos_grafos = sys.argv[1]

    # Verifica se o caminho fornecido existe
    if not os.path.isdir(caminho_dos_grafos):
        print(f"Erro: O diretório '{caminho_dos_grafos}' não existe.")
        sys.exit(1)

    # Processa os grafos e exibe os que têm vértices isolados
    resultado = carregar_grafos(caminho_dos_grafos)

    if resultado:
        print("Grafos que possuem vértices isolados:")
        for nome, isolados in resultado:
            print(f"{nome} -> Vértices isolados: {isolados}")
    else:
        print("Nenhum grafo com vértices isolados encontrado.")
