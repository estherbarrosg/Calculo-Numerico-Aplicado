"""
main.py

Atividade prática 02 - Métodos de Newton-Raphson e da Secante aplicados ao
balanço térmico de uma superfície aquecida.

O que este script faz:
    1. Valida a função f(T) comparando-a com a Tabela 1 do enunciado;
    2. Resolve f(T) = 0 pelo Newton-Raphson (com a derivada analítica) e pela
       secante (tratando f como caixa-preta), com o mesmo critério de parada;
    3. Confere o balanço de energia na temperatura encontrada;
    4. Compara os dois métodos;
    5. Imprime tudo no console e salva o mesmo relatório em
       outputs/resultados.txt.

O problema físico está em problema_termico.py e os métodos em metodos_abertos.py.

Como executar (dentro da pasta 02_Newton_Secante):
    python main.py
"""

import os  # biblioteca padrão do Python, usada só para montar o caminho da pasta outputs

from problema_termico import (
    funcao_balanco,
    derivada_balanco,
    perda_convectiva,
    perda_radiativa,
    potencia_termica,
)
from metodos_abertos import newton_raphson, secante

# ---------------------------------------------------------------------------
# Configuração dos métodos (dados de entrada da atividade)
# ---------------------------------------------------------------------------
chute_inicial_newton = 350.0    # T0 do Newton-Raphson [K]
chute_anterior_secante = 400.0  # T0 da secante [K]
chute_atual_secante = 425.0     # T1 da secante [K]
tolerancia = 1e-6               # critério de parada: |T_novo - T_atual| < tolerancia [K]

# Número máximo de iterações, para evitar laço infinito caso o método não convirja.
max_iteracoes = 50

# ---------------------------------------------------------------------------
# Dados de referência: Tabela 1 do enunciado (resíduo do balanço térmico)
# ---------------------------------------------------------------------------
temperaturas_tabela = [350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0]  # T [K]
residuos_tabela = [-11.867325, -7.204043, -2.062000, 3.624717,
                   9.926275, 16.917092, 24.675840]                       # f(T) [W]


def validar_com_tabela():
    """
    Compara o f(T) calculado pelo programa com os valores da Tabela 1 do enunciado.

    Saída:
        lista de strings (tabela de comparação e a maior diferença encontrada)
    """
    linhas = [f"{'T (K)':>7} | {'f(T) enunciado':>15} | {'f(T) calculado':>15} | {'Diferença':>10}",
              "-" * 59]
    maior_diferenca = 0.0

    for temperatura, residuo_enunciado in zip(temperaturas_tabela, residuos_tabela):
        residuo_calculado = funcao_balanco(temperatura)
        diferenca = abs(residuo_calculado - residuo_enunciado)
        maior_diferenca = max(maior_diferenca, diferenca)
        linhas.append(f"{temperatura:7.1f} | {residuo_enunciado:15.6f} | "
                      f"{residuo_calculado:15.6f} | {diferenca:10.2e}")

    linhas.append("")
    linhas.append(f"Maior diferença: {maior_diferenca:.2e} W "
                  f"(a tabela do enunciado tem 6 casas decimais)")
    return linhas


def formatar_tabela_iteracoes(historico):
    """
    Converte o histórico de um método em linhas de texto alinhadas em colunas.

    Entrada:
        historico : lista de tuplas devolvida por newton_raphson() ou secante()

    Saída:
        lista de strings (cabeçalho, separador e uma linha por iteração)
    """
    cabecalho = (f"{'k':>3} | {'T_atual (K)':>13} | {'f(T_atual) (W)':>15} | "
                 f"{'Erro (K)':>12}")
    linhas = [cabecalho, "-" * len(cabecalho)]

    for k, temperatura_atual, f_atual, erro in historico:
        linhas.append(f"{k:3d} | {temperatura_atual:13.6f} | {f_atual:15.6e} | {erro:12.6e}")

    return linhas


def montar_secao_metodo(titulo, temperatura, historico, convergiu):
    """
    Monta as linhas do relatório de UM método: tabela de iterações e resultados.

    Entradas:
        titulo      : nome do método, para o título da seção
        temperatura : temperatura de equilíbrio encontrada [K]
        historico   : lista de tuplas devolvida pelo método
        convergiu   : True se o critério de parada foi atendido

    Saída:
        lista de strings
    """
    linhas = [titulo.upper(), ""]
    linhas.extend(formatar_tabela_iteracoes(historico))
    linhas.append("")

    if convergiu:
        calor_convectivo = perda_convectiva(temperatura)
        calor_radiativo = perda_radiativa(temperatura)
        linhas.append(f"Temperatura de equilíbrio encontrada: {temperatura:.6f} K")
        linhas.append(f"Total de iterações: {len(historico)}")
        linhas.append(f"Resíduo: f(T) = {funcao_balanco(temperatura):.2e} W")
        linhas.append(f"Balanço de energia: convecção = {calor_convectivo:.6f} W | "
                      f"radiação = {calor_radiativo:.6f} W | "
                      f"soma = {calor_convectivo + calor_radiativo:.6f} W "
                      f"(potência fornecida P = {potencia_termica:.6f} W)")
    else:
        linhas.append(f"ATENÇÃO: o método não convergiu (máximo de {max_iteracoes} "
                      f"iterações ou divisão por zero). Última aproximação: "
                      f"T = {temperatura:.6f} K")

    return linhas


def main():
    """Resolve o problema térmico pelos dois métodos e gera o relatório."""
    temperatura_newton, historico_newton, convergiu_newton = newton_raphson(
        funcao_balanco,
        derivada_balanco,
        chute_inicial_newton,
        tolerancia,
        max_iteracoes,
    )
    temperatura_secante, historico_secante, convergiu_secante = secante(
        funcao_balanco,
        chute_anterior_secante,
        chute_atual_secante,
        tolerancia,
        max_iteracoes,
    )

    # Monta o relatório como lista de linhas (assim o mesmo texto é impresso
    # no console e salvo em arquivo, sem repetir código)
    relatorio = [
        "ATIVIDADE 02 - NEWTON-RAPHSON E SECANTE",
        "Temperatura de equilíbrio de uma superfície aquecida (convecção e radiação)",
        "",
        "VALIDAÇÃO DA FUNÇÃO f(T) COM A TABELA 1 DO ENUNCIADO",
        "",
    ]
    relatorio.extend(validar_com_tabela())
    relatorio.extend([
        "",
        f"Newton-Raphson: T0 = {chute_inicial_newton:.1f} K",
        f"Secante:        T0 = {chute_anterior_secante:.1f} K e T1 = {chute_atual_secante:.1f} K",
        f"Critério de parada: |T_novo - T_atual| < {tolerancia:.0e} K",
        "",
    ])
    relatorio.extend(montar_secao_metodo("Método de Newton-Raphson", temperatura_newton,
                                         historico_newton, convergiu_newton))
    relatorio.append("")
    relatorio.extend(montar_secao_metodo("Método da Secante", temperatura_secante,
                                         historico_secante, convergiu_secante))
    relatorio.append("")

    # Comparação final (só faz sentido se os dois métodos convergiram)
    if convergiu_newton and convergiu_secante:
        relatorio.append("COMPARAÇÃO DOS MÉTODOS")
        relatorio.append("")
        relatorio.append(f"{'Método':<25} | {'Iterações':>9} | {'T (K)':>12}")
        relatorio.append("-" * 52)
        relatorio.append(f"{'Newton-Raphson':<25} | {len(historico_newton):9d} | "
                         f"{temperatura_newton:12.6f}")
        relatorio.append(f"{'Secante':<25} | {len(historico_secante):9d} | "
                         f"{temperatura_secante:12.6f}")
        relatorio.append("")
        relatorio.append(f"Diferença entre as temperaturas: "
                         f"{abs(temperatura_newton - temperatura_secante):.2e} K")

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
