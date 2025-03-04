import os
import sys
import random
import numpy as np
import shutil
import networkx as nx


def listar_pastas(base_path):
    return [
        p for p in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, p))
    ]


def carregar_grafo(caminho_arquivo):
    G = nx.Graph()
    with open(caminho_arquivo, "r") as f:
        for linha in f:
            u, v = map(int, linha.strip().split())
            if u != v:  # Remover self-loops
                G.add_edge(u, v)  # NetworkX já evita múltiplas arestas
    return G


def calcular_densidade(G):
    return nx.density(G)


def selecionar_grafos(base_path, n_por_pasta_dict, output_path):
    pastas = listar_pastas(base_path)
    selecoes = {}
    os.makedirs(output_path, exist_ok=True)

    for pasta in pastas:
        caminho_pasta = os.path.join(base_path, pasta)
        arquivos = [f for f in os.listdir(caminho_pasta) if f.endswith(".txt")]
        n_por_pasta = n_por_pasta_dict.get(pasta, 0)

        if n_por_pasta == 0:
            continue

        # Calcular densidade para cada grafo
        grafos = []
        for arquivo in arquivos:
            caminho_arquivo = os.path.join(caminho_pasta, arquivo)
            G = carregar_grafo(caminho_arquivo)
            densidade = calcular_densidade(G)
            grafos.append((arquivo, densidade))

        if not grafos:
            continue

        # Ordenar pela densidade e normalizar
        grafos.sort(key=lambda x: x[1])
        densidades = np.array([d for _, d in grafos])

        # Aplicar distribuição normal apenas dentro de cada banco
        mean = np.mean(densidades)
        std = np.std(densidades)

        if std == 0:
            selecionados = random.sample(grafos, min(n_por_pasta, len(grafos)))
        else:
            probabilidades = np.exp(-((densidades - mean) ** 2) / (2 * std**2))
            probabilidades /= probabilidades.sum()
            selecionados = np.random.choice(
                len(grafos),
                size=min(n_por_pasta, len(grafos)),
                replace=False,
                p=probabilidades,
            )

        selecoes[pasta] = [grafos[i][0] for i in selecionados]

        # Criar pasta para cada banco dentro de 'instances'
        output_pasta = os.path.join(output_path, pasta)
        os.makedirs(output_pasta, exist_ok=True)

        # Copiar arquivos selecionados para a pasta de saída específica do banco
        for arquivo in selecoes[pasta]:
            origem = os.path.join(caminho_pasta, arquivo)
            destino = os.path.join(output_pasta, arquivo)
            shutil.copy(origem, destino)

    return selecoes


def main():
    if len(sys.argv) < 2:
        print("Uso: python selecionar_grafos.py <caminho_base>")
        sys.exit(1)

    base_path = sys.argv[1]
    output_path = os.path.join(base_path, "instances")

    n_por_pasta_dict = {}
    for pasta in listar_pastas(base_path):
        try:
            n = int(input(f"Quantos grafos deseja selecionar do banco '{pasta}'? "))
            n_por_pasta_dict[pasta] = n
        except ValueError:
            print(f"Entrada inválida para {pasta}, ignorando...")

    selecoes = selecionar_grafos(base_path, n_por_pasta_dict, output_path)

    for pasta, arquivos in selecoes.items():
        print(f"\nBanco: {pasta}")
        for arq in arquivos:
            print(f"  - {arq}")


if __name__ == "__main__":
    main()
