import networkx as nx
import matplotlib.pyplot as plt
import sys
import matplotlib

matplotlib.use("Agg")  # Configura o backend para não precisar de interface gráfica


def ler_arquivo(nome_arquivo):
    grafo = nx.Graph()  # Criar um grafo não direcionado
    try:
        with open(nome_arquivo, "r") as file:
            for linha in file:
                # Remover espaços e quebras de linha
                aresta = linha.strip().split()
                if len(aresta) == 2:
                    # Adicionar aresta ao grafo (deve ser no formato "nó1 nó2")
                    grafo.add_edge(int(aresta[0]), int(aresta[1]))
    except FileNotFoundError:
        print(f"O arquivo {nome_arquivo} não foi encontrado.")
        sys.exit(1)

    return grafo


def plotar_grafo(grafo):
    # Desenhar o grafo
    plt.figure(figsize=(10, 8))
    nx.draw(
        grafo,
        with_labels=True,
        node_size=500,
        node_color="skyblue",
        font_size=10,
        font_weight="bold",
    )
    plt.title("Grafo")

    # Salvar o gráfico como um arquivo de imagem
    plt.savefig("grafo.png")  # Salva o gráfico como PNG
    print("Grafo salvo como 'grafo.png'")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python script.py <caminho_para_arquivo>")
        sys.exit(1)

    nome_arquivo = sys.argv[1]
    grafo = ler_arquivo(nome_arquivo)
    plotar_grafo(grafo)
