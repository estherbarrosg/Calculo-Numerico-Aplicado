"""
main.py

Atividade prática 01 - Métodos da Bisseção e da Falsa Posição aplicados ao
escoamento em canal aberto.

O que este script faz:
    1. Resolve f(H) = Q_Manning(H) - Q = 0 pela bisseção e pela falsa posição,
       com o mesmo intervalo inicial e o mesmo critério de parada;
    2. Calcula a velocidade média U = Q / (B*H) para cada método;
    3. Compara os dois métodos;
    4. Imprime tudo no console e salva o mesmo relatório em
       outputs/resultados.txt.

O problema físico está em problema_canal.py e os métodos em
metodos_intervalares.py.

Como executar (dentro da pasta 01_Bissecao_Falsa_Posicao):
    python main.py
"""

import os  # biblioteca padrão do Python, usada só para montar o caminho da pasta outputs

from problema_canal import funcao_objetivo, vazao_manning, vazao_conhecida, velocidade_media
from metodos_intervalares import bissecao, falsa_posicao

# ---------------------------------------------------------------------------
# Configuração dos métodos (dados de entrada da atividade)
# ---------------------------------------------------------------------------
limite_inferior_inicial = 0.5   # Hl: limite inferior inicial do intervalo [m]
limite_superior_inicial = 1.0   # Hu: limite superior inicial do intervalo [m]
tolerancia = 1e-6               # critério de parada: |f(Hr)| < tolerancia [m^3/s]

# Número máximo de iterações, para evitar laço infinito em caso de erro.
max_iteracoes = 50


def formatar_tabela_iteracoes(historico):
    """
    Converte o histórico de um método em linhas de texto alinhadas em colunas.

    Entrada:
        historico : lista de tuplas devolvida por bissecao() ou falsa_posicao()

    Saída:
        lista de strings (cabeçalho, separador e uma linha por iteração)
    """
    cabecalho = (f"{'k':>3} | {'Hl':>9} | {'Hu':>9} | {'Hr':>9} | "
                 f"{'f(Hr)':>12} | {'Erro (%)':>10}")
    linhas = [cabecalho, "-" * len(cabecalho)]

    for k, h_inferior, h_superior, h_r, f_h_r, erro in historico:
        # Na primeira iteração não existe erro, então escrevemos "---"
        if erro is None:
            texto_erro = f"{'---':>10}"
        else:
            texto_erro = f"{erro:10.4f}"

        # f(Hr) em notação científica para mostrar valores muito pequenos
        linhas.append(f"{k:3d} | {h_inferior:9.6f} | {h_superior:9.6f} | "
                      f"{h_r:9.6f} | {f_h_r:12.4e} | {texto_erro}")

    return linhas


def montar_secao_metodo(titulo, profundidade, velocidade, historico, convergiu):
    """
    Monta as linhas do relatório de UM método: tabela de iterações e resultados.

    Entradas:
        titulo      : nome do método, para o título da seção
        profundidade: raiz H encontrada (None se o intervalo inicial é inválido)
        velocidade  : velocidade média U [m/s] (None se não houve convergência)
        historico   : lista de tuplas devolvida pelo método
        convergiu   : True se o critério de parada foi atendido

    Saída:
        lista de strings
    """
    linhas = [titulo.upper(), ""]

    if profundidade is None:
        # Intervalo sem mudança de sinal: o método não pode ser aplicado
        linhas.append("ERRO: f(Hl)*f(Hu) >= 0. Não há garantia de raiz neste "
                      "intervalo; escolha outro intervalo inicial.")
        return linhas

    linhas.extend(formatar_tabela_iteracoes(historico))
    linhas.append("")

    if convergiu:
        linhas.append(f"Raiz encontrada: profundidade H = {profundidade:.6f} m")
        linhas.append(f"Total de iterações: {len(historico)}")
        linhas.append(f"Verificação: Q(H) por Manning = "
                      f"{vazao_manning(profundidade):.6f} m^3/s "
                      f"(vazão prescrita = {vazao_conhecida:.6f} m^3/s)")
        linhas.append(f"Velocidade média: U = Q/(B*H) = {velocidade:.4f} m/s")
    else:
        linhas.append(f"ATENÇÃO: o método não convergiu em {max_iteracoes} "
                      f"iterações. Última aproximação: H = {profundidade:.6f} m")

    return linhas


def main():
    """Resolve o problema do canal pelos dois métodos e gera o relatório."""
    # Valores de f nos extremos, apenas para exibir a verificação inicial
    f_inferior_inicial = funcao_objetivo(limite_inferior_inicial)
    f_superior_inicial = funcao_objetivo(limite_superior_inicial)

    # Monta o relatório como lista de linhas (assim o mesmo texto é impresso
    # no console e salvo em arquivo, sem repetir código)
    relatorio = [
        "ATIVIDADE 01 - BISSEÇÃO E FALSA POSIÇÃO",
        "Profundidade da lâmina d'água em canal aberto retangular (Manning)",
        "",
        f"Intervalo inicial: Hl = {limite_inferior_inicial:.6f} m | "
        f"Hu = {limite_superior_inicial:.6f} m",
        f"f(Hl) = {f_inferior_inicial:.6f} | f(Hu) = {f_superior_inicial:.6f} | "
        f"f(Hl)*f(Hu) = {f_inferior_inicial * f_superior_inicial:.6f}",
        f"Critério de parada: |f(Hr)| < {tolerancia:.0e}",
        "",
    ]

    # Cada método é uma tupla (nome, função do método)
    metodos = [("Método da Bisseção", bissecao),
               ("Método da Falsa Posição", falsa_posicao)]

    # Resultados dos métodos que convergiram, para a comparação final
    resultados = []

    for nome, metodo in metodos:
        profundidade, historico, convergiu = metodo(
            funcao_objetivo,
            limite_inferior_inicial,
            limite_superior_inicial,
            tolerancia,
            max_iteracoes,
        )

        velocidade = velocidade_media(profundidade) if convergiu else None
        relatorio.extend(montar_secao_metodo(nome, profundidade, velocidade,
                                             historico, convergiu))
        relatorio.append("")

        if convergiu:
            resultados.append((nome, profundidade, len(historico), velocidade))

    # Comparação final (só faz sentido se os dois métodos convergiram)
    if len(resultados) == len(metodos):
        relatorio.append("COMPARAÇÃO DOS MÉTODOS")
        relatorio.append("")
        relatorio.append(f"{'Método':<25} | {'Iterações':>9} | {'H (m)':>10} | {'U (m/s)':>9}")
        relatorio.append("-" * 61)
        for nome, profundidade, n_iteracoes, velocidade in resultados:
            relatorio.append(f"{nome:<25} | {n_iteracoes:9d} | "
                             f"{profundidade:10.6f} | {velocidade:9.6f}")

        # Diferença entre as soluções dos dois métodos
        diferenca_profundidade = abs(resultados[0][1] - resultados[1][1])
        diferenca_velocidade = abs(resultados[0][3] - resultados[1][3])
        relatorio.append("")
        relatorio.append(f"Diferença entre as profundidades: {diferenca_profundidade:.2e} m")
        relatorio.append(f"Diferença entre as velocidades:   {diferenca_velocidade:.2e} m/s")

    texto_relatorio = "\n".join(relatorio)
    print(texto_relatorio)

    # Salva o relatório em outputs/, ao lado deste script, não importa de
    # qual pasta o programa foi executado
    pasta_do_script = os.path.dirname(os.path.abspath(__file__))
    pasta_de_saida = os.path.join(pasta_do_script, "outputs")
    os.makedirs(pasta_de_saida, exist_ok=True)  # cria a pasta se ela não existir
    caminho_do_arquivo = os.path.join(pasta_de_saida, "resultados.txt")
    with open(caminho_do_arquivo, "w", encoding="utf-8") as arquivo:
        arquivo.write(texto_relatorio + "\n")


# Este bloco só executa quando o arquivo é rodado diretamente
# (python main.py). Se outro arquivo importar este módulo, por exemplo o
# graficos.py para reaproveitar a configuração, o main() não roda sozinho.
if __name__ == "__main__":
    main()
