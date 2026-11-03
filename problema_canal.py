"""
problema_canal.py

Definição do problema físico da atividade prática 01 (Bisseção e Falsa Posição):
escoamento uniforme de água em um canal aberto retangular, descrito pela
equação de Manning.

Este arquivo contém SOMENTE o problema (dados e equações), sem nenhum método
numérico. Assim, o mesmo problema pode ser resolvido por métodos diferentes
(bisseção, falsa posição) importando as funções daqui, sem repetir código
(princípio DRY do Manual de Boas Práticas).
"""

import math  # usado apenas para a raiz quadrada (math.sqrt)

# ---------------------------------------------------------------------------
# Parâmetros do problema (dados fornecidos no enunciado)
# ---------------------------------------------------------------------------
declividade = 0.0002        # S: inclinação do canal [adimensional, m/m]
largura_canal = 20.0        # B: largura do canal retangular [m]
rugosidade_manning = 0.03   # n: coeficiente de rugosidade de Manning [s/m^(1/3)]
vazao_conhecida = 5.0       # Q: vazão volumétrica prescrita [m^3/s]


def vazao_manning(profundidade):
    """
    Calcula a vazão Q prevista pela equação de Manning para uma profundidade H.

        Q(H) = (sqrt(S) / n) * (B*H)^(5/3) / (B + 2*H)^(2/3)

    Entrada:
        profundidade : profundidade H da lâmina d'água [m]

    Saída:
        vazão volumétrica calculada [m^3/s]
    """
    area_molhada = largura_canal * profundidade             # A_c = B*H [m^2]
    perimetro_molhado = largura_canal + 2.0 * profundidade  # P = B + 2H [m]
    coeficiente_manning = math.sqrt(declividade) / rugosidade_manning

    return (coeficiente_manning
            * area_molhada ** (5.0 / 3.0)
            / perimetro_molhado ** (2.0 / 3.0))


def funcao_objetivo(profundidade):
    """
    Função f(H) cuja raiz é a profundidade procurada.

        f(H) = Q_Manning(H) - Q

    Quando f(H) = 0, a vazão calculada pela equação de Manning é igual à
    vazão prescrita, ou seja, H é a profundidade real do escoamento.

    Entrada:
        profundidade : profundidade H da lâmina d'água [m]

    Saída:
        f(H) [m^3/s]
    """
    return vazao_manning(profundidade) - vazao_conhecida


def velocidade_media(profundidade):
    """
    Calcula a velocidade média do escoamento a partir da vazão e da área.

        U = Q / (B*H)

    Entrada:
        profundidade : profundidade H da lâmina d'água [m]

    Saída:
        velocidade média U [m/s]
    """
    return vazao_conhecida / (largura_canal * profundidade)
