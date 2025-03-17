import os
import sys
import argparse
import shutil
import networkx as nx
import numpy as np
from scipy.stats import norm


def calcular_densidade(arquivo):
    """Lê um grafo a partir de um arquivo e calcula sua densidade."""
    G = nx.Graph()
    with open(arquivo, "r") as f:
        for line in f:
            u, v = map(int, line.split())
            G.add_edge(u, v)
    return nx.density(G)


def selecionar_grafos_por_densidade(arquivos, quantidade):
    """Seleciona uma quantidade de grafos seguindo uma distribuição gaussiana com base na densidade."""
    if len(arquivos) <= quantidade:
        return arquivos  # Se há menos ou igual arquivos do que queremos, retorna todos

    # Calcula a densidade de cada grafo
    densidades = [(arquivo, calcular_densidade(arquivo)) for arquivo in arquivos]

    # Ordena os grafos por densidade (do menor para o maior)
    densidades.sort(key=lambda x: x[1])

    # Extrai somente os valores das densidades para normalização
    valores_densidade = np.array([d[1] for d in densidades])

    # Normaliza as densidades para gerar uma distribuição gaussiana
    media = np.mean(valores_densidade)
    desvio_padrao = (
        np.std(valores_densidade) if np.std(valores_densidade) > 0 else 1
    )  # Evita divisão por zero
    probabilidades = norm.pdf(
        valores_densidade, media, desvio_padrao
    )  # Função gaussiana

    # Normaliza as probabilidades para que somem 1
    probabilidades /= probabilidades.sum()

    # Seleciona índices de acordo com a distribuição gaussiana
    indices_selecionados = np.random.choice(
        len(arquivos), size=quantidade, replace=False, p=probabilidades
    )

    # Retorna os arquivos selecionados
    return [densidades[i][0] for i in indices_selecionados]


def calcular_distribuicao(base_dir, total_instances=None):
    """Percorre as subpastas e calcula a porcentagem de grafos por base.
    Se --instances for passado, escolhe os grafos proporcionalmente com base na densidade e copia para 'instances'.
    """
    distribuicao = {}
    total_grafos = 0

    # Percorre todas as subpastas dentro do diretório fornecido
    for root, _, files in os.walk(base_dir):
        if root != base_dir:  # Garante que estamos nas subpastas
            nome_base = os.path.basename(root)
            arquivos_grafos = [
                os.path.join(root, f) for f in files if f.endswith(".txt")
            ]
            qtd_grafos = len(arquivos_grafos)

            if qtd_grafos > 0:
                distribuicao[nome_base] = arquivos_grafos
                total_grafos += qtd_grafos

    if total_grafos == 0:
        print("Nenhum grafo encontrado no diretório.")
        return

    # Exibe as estatísticas gerais
    print(f"📊 Total de grafos encontrados: {total_grafos}\n")
    print("📌 Distribuição por base:")

    proporcoes = {}
    for base, arquivos in distribuicao.items():
        qtd = len(arquivos)
        porcentagem = (qtd / total_grafos) * 100
        print(f"  - {base}: {qtd} grafos ({porcentagem:.2f}%)")
        proporcoes[base] = porcentagem

    # Se --instances foi passado, seleciona os grafos com base na densidade e copia para a pasta 'instances'
    if total_instances:
        print(
            f"\n🎯 Selecionando {total_instances} grafos de forma proporcional e copiando para a pasta 'instances/'"
        )

        # Criar diretório para armazenar instâncias
        instances_dir = os.path.join(base_dir, "instances")
        os.makedirs(instances_dir, exist_ok=True)

        # Distribuição proporcional das instâncias
        selecionados = {
            base: round((porcentagem / 100) * total_instances)
            for base, porcentagem in proporcoes.items()
        }

        # Ajuste final para garantir que o total seja exatamente `total_instances`
        ajuste = total_instances - sum(selecionados.values())
        bases_ordenadas = sorted(
            selecionados.items(), key=lambda x: -x[1]
        )  # Ordena para evitar arredondamento a zero

        for i in range(abs(ajuste)):
            if ajuste > 0:
                selecionados[bases_ordenadas[i % len(bases_ordenadas)][0]] += 1
            elif ajuste < 0:
                selecionados[bases_ordenadas[i % len(bases_ordenadas)][0]] -= 1

        # Seleciona e copia os grafos para a pasta 'instances'
        for base, qtd in selecionados.items():
            arquivos = distribuicao[base]
            arquivos_selecionados = selecionar_grafos_por_densidade(arquivos, qtd)

            # Criar subpasta para a base dentro de 'instances'
            base_instances_dir = os.path.join(instances_dir, base)
            os.makedirs(base_instances_dir, exist_ok=True)

            # Copiar os arquivos selecionados
            for arquivo in arquivos_selecionados:
                shutil.copy(arquivo, base_instances_dir)

            print(f"  - {base}: {qtd} grafos copiados para 'instances/{base}/'")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analisa a distribuição de grafos por base e seleciona uma quantidade proporcional baseada na densidade."
    )
    parser.add_argument("caminho", help="Caminho da pasta contendo os bancos de grafos")
    parser.add_argument(
        "--instances",
        type=int,
        help="Número total de grafos desejado, distribuído proporcionalmente e escolhido com base na densidade",
    )

    args = parser.parse_args()

    # Verifica se o caminho fornecido existe
    if not os.path.isdir(args.caminho):
        print(f"❌ Erro: O diretório '{args.caminho}' não existe.")
        sys.exit(1)

    # Processa e exibe a distribuição dos grafos
    calcular_distribuicao(args.caminho, args.instances)
