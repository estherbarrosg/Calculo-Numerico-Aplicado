# -*- coding: utf-8 -*-
# ==============================================================================
# SIMULADO 3 - 1ª PROVA DE CÁLCULO NUMÉRICO APLICADO
# Aleta de um dissipador de componente eletrônico
#
# COMO LER ESTE ARQUIVO: tudo que foi corrigido em relação ao seu código está
# marcado com um comentário que começa com  ">>> CORREÇÃO n".  Procure por esse
# texto (Ctrl+F) para ver cada mudança.
#
# Índice das correções:
#   1  import sys (faltava; sem ele o aviso de pivô nulo daria erro)
#   2  valores de w e rho (estavam trocados, os erros se cancelavam)
#   3  funções do modelo: uma massa só (com L/1000), benefício e custo com c_m
#   4  rótulos da Tabela 1 (L em mm) e da impressão do sistema (b é q em W)
#   5  Q2: nome do arquivo, cabeçalho com tab, faixa 20..80 e só 2 colunas
#   6  Q2: desvio máximo comparado com o q MEDIDO (em W) e só nos 4 ensaios
#   7  Q2: erro absoluto da validação com mais casas decimais
#   8  Q3: critério de parada testado em toda iteração + aviso de máx. iterações
#   9  Q4: sem linha tracejada no cabeçalho do arquivo de dados
#   10 Q5/Q6: UMA função de razão áurea, com a comparação certa para minimizar
#   11 Q5: prints exigidos (L_B*, B, q, iterações, avaliações)
#   12 Q6(a): prints exigidos e respostas calculadas pelo programa
#   13 Q6(b): tabela usando a função + verificação por varredura (extra)
#
# Unidades: L em mm | q em W | massa em kg | B em reais | c em reais por W
# ==============================================================================

import math
import sys   # >>> CORREÇÃO 1: sem este import, o sys.exit() da Questão 1 daria NameError

# ==============================================================================
# DADOS GERAIS E CONSTANTES DO PROBLEMA
# ==============================================================================
w = 0.030      # >>> CORREÇÃO 2: largura da aleta [m]. Estava 0.003 (faltava um zero)
t = 0.0008     # espessura da aleta [m]
h = 50.0       # coeficiente de convecção [W/(m2.K)]
rho = 2700.0   # >>> CORREÇÃO 2: massa específica do alumínio [kg/m3]. Estava 27000
c_q = 0.06     # valor de cada watt dissipado [reais por W]
c_fix = 0.15   # custo fixo de fabricação por aleta [reais]
q_req = 6.0    # potência mínima exigida de cada aleta [W]
c_m = 45.0     # preço do material [reais por kg]

# ==============================================================================
# FUNÇÕES DO MODELO (usadas em várias questões)
# ==============================================================================
# OBS.: Taxa_calor usa o vetor x calculado na Questão 1. Por isso ela só pode ser
# CHAMADA depois da Questão 1 (definir a função antes não tem problema).

def Taxa_calor(L):
    """q(L) em W, com L em mm (equação 1)."""
    return x[0]*L**3 + x[1]*L**2 + x[2]*L + x[3]

# >>> CORREÇÃO 3: antes havia três funções (massa, Beneficio, custo), sendo que
# massa() não convertia L para metros, Beneficio() nunca era usada e custo() usava
# a variável m_al, que não existe (daria NameError). Agora são três funções
# limpas e todas recebem o preço do material (cm) como argumento.

def massa_aleta(L):
    """Massa de uma aleta em kg (equação 2). L chega em mm e é convertido para m."""
    return rho * w * t * (L / 1000.0)

def beneficio(L, cm):
    """Benefício líquido por aleta em reais (equação 3)."""
    return c_q * Taxa_calor(L) - cm * massa_aleta(L)

def custo_por_watt(L, cm):
    """Custo por watt dissipado em reais por W (equação 4)."""
    return (c_fix + cm * massa_aleta(L)) / Taxa_calor(L)

# ==============================================================================
# Ensaios de bancada, NA ORDEM DA TABELA 1: (L em mm, q medido em W)
# >>> CORREÇÃO 4: o comentário dizia "L em m", mas L está em mm
# ==============================================================================
tabela_1 = [
    (20.0, 3.47),  # Ensaio 1 (ajuste)
    (40.0, 5.63),  # Ensaio 2 (ajuste)
    (60.0, 6.64),  # Ensaio 3 (ajuste)
    (80.0, 7.05),  # Ensaio 4 (ajuste)
    (50.0, 6.24),  # Ensaio 5 (validação: NÃO entra no ajuste)
]

# ==============================================================================
# QUESTÃO 1 - DETERMINAÇÃO DOS COEFICIENTES DO POLINÔMIO
# Eliminação gaussiana com pivoteamento parcial
# ==============================================================================
print("=" * 75)
print("QUESTÃO 1 - DETERMINAÇÃO DOS COEFICIENTES DO POLINÔMIO")
print("=" * 75)

# Cada ensaio (L, q) vira UMA equação:  a3*L^3 + a2*L^2 + a1*L + a0 = q
# Então a linha de A é [L^3, L^2, L, 1] e o b dessa linha é o q medido.
A = []
b = []
for L_med, q_med in tabela_1[:4]:     # [:4] = só os 4 ensaios de ajuste
    A.append([L_med**3, L_med**2, L_med, 1.0])
    b.append(q_med)
n = len(b)   # número de equações (= número de incógnitas = 4)

print("Sistema montado (equações na ordem da Tabela 1):")
# >>> CORREÇÃO 4: o rótulo dizia "b (eta)", mas aqui b é o q medido, em W
print("      coef. a3     coef. a2   coef. a1  coef. a0 |  b (q em W)")
for i in range(n):
    print(f"  {A[i][0]:10.1f} {A[i][1]:12.1f} {A[i][2]:10.1f} {A[i][3]:8.1f} | {b[i]:.3f}")

# Cópias para calcular o resíduo no final (a eliminação ESTRAGA A e b)
A_orig = [linha[:] for linha in A]
b_orig = b[:]

LIMITE_PIVO = 1e-5  # pivô com valor absoluto menor que isso é tratado como nulo

# --- FASE 1: Eliminação progressiva ----------------------------------------------
print("\nEliminação progressiva (pivoteamento parcial):")
for k in range(n - 1):          # k = coluna que vamos zerar abaixo da diagonal

    # Procura, da linha k para baixo, a linha com maior valor absoluto na coluna k
    max_index = k
    maior_valor = abs(A[k][k])
    for i in range(k + 1, n):
        if abs(A[i][k]) > maior_valor:
            maior_valor = abs(A[i][k])
            max_index = i

    # Se até o MAIOR valor da coluna é (quase) zero, não tem como continuar
    if maior_valor <= LIMITE_PIVO:
        sys.exit(f"ERRO: pivô nulo na etapa k = {k + 1}. Sistema singular, programa encerrado.")

    # Troca de linhas (em A e em b) para o maior valor ficar na diagonal
    if max_index != k:
        A[k], A[max_index] = A[max_index], A[k]
        b[k], b[max_index] = b[max_index], b[k]

    # Mostra o que aconteceu nesta etapa (+1 porque Python conta a partir de 0)
    print(f"  Etapa k = {k + 1}: linha pivô = {max_index + 1} | pivô = {A[k][k]:.6e}")

    # Eliminação: zera a coluna k nas linhas abaixo
    for i in range(k + 1, n):
        m = A[i][k] / A[k][k]            # multiplicador
        for j in range(k, n):
            A[i][j] -= m * A[k][j]       # linha i = linha i - m * linha k
        b[i] -= m * b[k]

# O último pivô (A[n-1][n-1]) também precisa ser diferente de zero
if abs(A[n - 1][n - 1]) <= LIMITE_PIVO:
    sys.exit("ERRO: último pivô nulo. Sistema singular, programa encerrado.")

# --- FASE 2: Substituição regressiva (de baixo para cima) -----------------------
x = [0.0] * n
x[n - 1] = b[n - 1] / A[n - 1][n - 1]
for i in range(n - 2, -1, -1):
    soma = b[i]
    for j in range(i + 1, n):
        soma -= A[i][j] * x[j]
    x[i] = soma / A[i][i]

# --- FASE 3: Resíduo  r = A_orig*x - b_orig  e impressão ------------------------
residuo = []
for i in range(n):
    soma = 0.0
    for j in range(n):
        soma += A_orig[i][j] * x[j]
    residuo.append(soma - b_orig[i])
norma_residuo = max(abs(valor) for valor in residuo)

nomes = ["a3", "a2", "a1", "a0"]
print("\nCoeficientes encontrados:")
for i in range(n):
    print(f"  {nomes[i]} = {x[i]:.6e}")
print(f"Norma infinito do resíduo: {norma_residuo:.3e}")

# COMENTÁRIO EXIGIDO (por que pivotear mesmo sem pivô nulo): se A[k][k] for muito
# pequeno em relação aos outros elementos da coluna, o multiplicador fica enorme.
# Multiplicar a linha pivô por um número grande e subtraí-la das de baixo amplifica
# os erros de arredondamento. Aqui os elementos de A vão de 1 a mais de 5x10^5, então
# escolher sempre o MAIOR pivô disponível mantém os multiplicadores pequenos.

# ==============================================================================
# QUESTÃO 2 - ARQUIVO DE RESULTADOS E COMPARAÇÃO COM OS ENSAIOS
# ==============================================================================
print("\n" + "=" * 75)
print("QUESTÃO 2 - ARQUIVO resultado.dat E DESVIOS")
print("=" * 75)

# >>> CORREÇÃO 5: o arquivo se chama resultado.dat (estava "resultados.dat"), o
# cabeçalho precisa de \t (estava "\q", que nem é um comando válido), o laço vai de
# 20 a 80 (estava range(81), que começa em 0 e gera valores negativos de q) e o
# arquivo tem só 2 colunas (havia uma 3ª coluna, q_med, que era sobra de outro laço).
with open("resultado.dat", "w") as arquivo:
    arquivo.write("L(mm)\tq(W)\n")
    for L in range(20, 81):   # L = 20, 21, ..., 80 mm
        arquivo.write(f"{L}\t{Taxa_calor(L):.4f}\n")
print("Arquivo 'resultado.dat' gerado com sucesso!")

# >>> CORREÇÃO 6: antes a conta era abs(Taxa_calor_calc - L_med)*100, ou seja,
# comparava o q do polinômio com o COMPRIMENTO L (e não com o q medido) e ainda
# multiplicava por 100. O desvio correto é |q do polinômio - q medido|, em W, e só
# nos 4 ensaios de ajuste (tabela_1[:4]); o 5º ensaio é a validação, tratada abaixo.
maior_desvio = 0.0
for L_med, q_medido in tabela_1[:4]:
    q_calc = Taxa_calor(L_med)
    desvio = abs(q_calc - q_medido)
    if desvio > maior_desvio:
        maior_desvio = desvio
    print(f"L = {L_med:4.1f} mm | q medido: {q_medido:.4f} W | "
          f"q do polinômio: {q_calc:.4f} W | desvio: {desvio:.2e} W")
print(f"Maior desvio nos ensaios de ajuste: {maior_desvio:.2e} W")

print("-" * 75)
# Ensaio de validação (5º da lista, índice 4): NÃO foi usado no ajuste
L_val, q_val = tabela_1[4]
q_calc = Taxa_calor(L_val)
ea = abs(q_calc - q_val)
er = ea / q_val
print(f"Validação em L = {L_val:.0f} mm:")
# >>> CORREÇÃO 7: com :.1f e :.2f, o erro absoluto (0,0044 W) aparecia como 0,00.
print(f"  q do polinômio = {q_calc:.4f} W | q medido = {q_val:.2f} W")
print(f"  Erro absoluto = {ea:.4f} W | Erro relativo = {er * 100:.2f} %")
print("-" * 75)

# q(0) e comentário
q_0 = Taxa_calor(0)
print(f"Para L = 0 mm o polinômio dá q = {q_0:.4f} W")
# COMENTÁRIO: o polinômio só vale em 20 <= L <= 80 mm. Ele foi ajustado por pontos
# dentro desse intervalo (interpolação) e, fora dele, q(L) pode dar valores
# negativos, como q(0) = -0,39 W, fisicamente impossível.

# ==============================================================================
# QUESTÃO 3 - COMPRIMENTO MÍNIMO DA ALETA (falsa posição)
# Raiz de f(L) = q(L) - 6,0
# ==============================================================================
print("\n" + "=" * 75)
print("QUESTÃO 3 - COMPRIMENTO MÍNIMO (Falsa Posição)")
print("=" * 75)

def f(L):
    """f(L) = q(L) - q_req. A raiz é o comprimento em que q vale 6,0 W."""
    return Taxa_calor(L) - q_req

L_I = 20.0          # extremo inferior [mm]
L_U = 80.0          # extremo superior [mm]
tolerancia = 1e-6   # |f(L_r)| < tolerancia [W]
max_iter = 200

if f(L_I) * f(L_U) > 0:
    print("Erro: não há mudança de sinal nesse intervalo")
else:
    print("-" * 65)
    # >>> CORREÇÃO 8 (rótulo): a coluna se chamava "Hr", nome sobra de outro problema
    print(f"{'i':>2} | {'L_I':>10} | {'L_U':>10} | {'L_r':>10} | {'f(L_r)':>12} | {'Erro (%)':>10}")
    print("-" * 65)

    L_antigo = L_I
    for i in range(1, max_iter + 1):
        L_r = (L_I * f(L_U) - L_U * f(L_I)) / (f(L_U) - f(L_I))
        f_r = f(L_r)

        print(f"{i:2d} | {L_I:10.6f} | {L_U:10.6f} | {L_r:10.6f} | {f_r:12.6f} ", end="")
        if i > 1:
            ea_pct = abs((L_r - L_antigo) / L_antigo) * 100
            print(f"| {ea_pct:9.6f}%")
        else:
            print("|        ---")

        # >>> CORREÇÃO 8: antes o critério de parada ficava DENTRO do "if i > 1",
        # então só era testado a partir da 2ª iteração. Agora é testado sempre.
        if abs(f_r) < tolerancia:
            print("-" * 65)
            print(f"Raiz encontrada! Comprimento mínimo da aleta: L_s = {L_r:.6f} mm")
            print(f"Total de iterações: {i}")
            break

        # Mantém o extremo que preserva a mudança de sinal
        if f(L_I) * f_r < 0:
            L_U = L_r
        else:
            L_I = L_r
        L_antigo = L_r
    else:
        # >>> CORREÇÃO 8: o "else" do for roda se o laço terminou SEM break
        print("AVISO: número máximo de iterações atingido sem convergir.")

# ==============================================================================
# QUESTÃO 4 - BENEFÍCIO E CUSTO POR WATT (arquivo com 4 colunas)
# ==============================================================================
print("\n" + "=" * 75)
print("QUESTÃO 4 - ARQUIVO resultado.dat COM q, B E c")
print("=" * 75)

with open("resultado.dat", "w") as arquivo:
    # >>> CORREÇÃO 9: removida a linha tracejada ("-----") do cabeçalho. Linhas
    # que não são dados atrapalham quem for ler o arquivo depois (gráficos etc.).
    arquivo.write("# L(mm)   q(W)      B(reais)    c(reais/W)\n")
    for L in range(20, 81):
        q = Taxa_calor(L)
        B = beneficio(L, c_m)         # usa as funções corrigidas (CORREÇÃO 3)
        c = custo_por_watt(L, c_m)
        arquivo.write(f"  {L:<8d} {q:<9.4f} {B:<11.5f} {c:<11.5f}\n")
print("Arquivo 'resultado.dat' gerado com 4 colunas (L, q, B, c).")

# ==============================================================================
# QUESTÕES 5 E 6 - OTIMIZAÇÃO (razão áurea)
# ==============================================================================

# >>> CORREÇÃO 10: ESTE ERRO DERRUBOU AS QUESTÕES 5 E 6.
# Antes, o laço da razão áurea estava copiado duas vezes dentro do for, com a
# comparação "if f1 > f2: b = e2", que é a regra para MAXIMIZAR. Mas você estava
# MINIMIZANDO (-B e c). Para minimizar a regra é a oposta:
#       se f1 < f2, o menor valor está à esquerda  -> descarta a direita (b = x2)
#       senão                                      -> descarta a esquerda (a = x1)
# Com a regra errada o programa procurava o PIOR ponto (por isso os resultados
# ficavam colados em 20 ou 80). Agora existe UMA função, usada para B e para c.
# Para maximizar, a função minimiza o negativo (sinal = -1).
def razao_aurea(funcao, a, b, maximizar, tol=1e-6, max_iter=200):
    """Razão áurea. Devolve (x_otimo, valor_em_x_otimo, iterações, avaliações).
    funcao = função de UMA variável (passada sem parênteses).
    maximizar = True procura o máximo; False procura o mínimo."""
    phi = (math.sqrt(5.0) - 1.0) / 2.0
    if maximizar:
        sinal = -1.0      # maximizar f é o mesmo que minimizar -f
    else:
        sinal = 1.0
    x1 = b - phi * (b - a)
    x2 = a + phi * (b - a)
    f1 = sinal * funcao(x1)
    f2 = sinal * funcao(x2)
    avaliacoes = 2        # contador: as duas avaliações iniciais
    iteracoes = 0
    while (b - a) > tol and iteracoes < max_iter:
        iteracoes += 1
        # A partir daqui o objetivo é SEMPRE minimizar (o sinal já cuidou do máximo)
        if f1 < f2:
            b = x2                       # o mínimo está em [a, x2]
            x2 = x1
            f2 = f1                      # reaproveita a avaliação
            x1 = b - phi * (b - a)
            f1 = sinal * funcao(x1)      # só 1 avaliação nova
        else:
            a = x1                       # o mínimo está em [x1, b]
            x1 = x2
            f1 = f2                      # reaproveita a avaliação
            x2 = a + phi * (b - a)
            f2 = sinal * funcao(x2)      # só 1 avaliação nova
        avaliacoes += 1
    x_otimo = (a + b) / 2.0
    valor = funcao(x_otimo)
    avaliacoes += 1       # convenção: conta também a avaliação final
    return x_otimo, valor, iteracoes, avaliacoes

# Verificação por varredura (EXTRA, mas muito útil na prova): testa L = 20, 20.5,
# ..., 80 e devolve o melhor L. Se a razão áurea der um L muito diferente, tem erro.
def varredura(funcao, maximizar):
    melhor_L = 20.0
    melhor_valor = funcao(20.0)
    L = 20.0
    while L <= 80.0:
        valor = funcao(L)
        if maximizar and valor > melhor_valor:
            melhor_valor = valor
            melhor_L = L
        if (not maximizar) and valor < melhor_valor:
            melhor_valor = valor
            melhor_L = L
        L += 0.5
    return melhor_L

# Funções de UMA variável (L) com o preço c_m = 45 "embutido"
def B_base(L):
    return beneficio(L, c_m)

def c_base(L):
    return custo_por_watt(L, c_m)

# ------------------------------------------------------------------------------
# QUESTÃO 5 - comprimento de máximo benefício
# ------------------------------------------------------------------------------
print("\n" + "=" * 75)
print("QUESTÃO 5 - COMPRIMENTO DE MÁXIMO BENEFÍCIO (Razão Áurea, c_m = 45)")
print("=" * 75)
# >>> CORREÇÃO 11: o enunciado pede estes prints; antes só saía a tabela final
L_B, B_max, it_B, aval_B = razao_aurea(B_base, 20.0, 80.0, True)
print(f"L_B* = {L_B:.6f} mm")
print(f"B(L_B*) = {B_max:.6f} reais")
print(f"q(L_B*) = {Taxa_calor(L_B):.4f} W")
print(f"Iterações: {it_B} | Avaliações de B: {aval_B}")
print(f"(Conferência por varredura: melhor L = {varredura(B_base, True):.1f} mm)")

# ------------------------------------------------------------------------------
# QUESTÃO 6(a) - comprimento de custo mínimo
# ------------------------------------------------------------------------------
print("\n" + "=" * 75)
print("QUESTÃO 6(a) - COMPRIMENTO DE CUSTO MÍNIMO (Razão Áurea, c_m = 45)")
print("=" * 75)
L_C, c_min, it_C, aval_C = razao_aurea(c_base, 20.0, 80.0, False)
print(f"L_c* = {L_C:.6f} mm")
print(f"c(L_c*) = {c_min:.6f} reais por W")
print(f"q(L_c*) = {Taxa_calor(L_C):.4f} W")
print(f"Iterações: {it_C} | Avaliações de c: {aval_C}")
print(f"(Conferência por varredura: melhor L = {varredura(c_base, False):.1f} mm)")

# >>> CORREÇÃO 12: antes as respostas estavam digitadas no comentário (copiadas do
# gabarito) e NÃO batiam com o que o seu programa imprimia. Agora o programa decide.
print("\nRespostas calculadas pelo programa:")
if abs(L_C - L_B) < 1e-3:
    print("  L_c* coincide com L_B*.")
else:
    print(f"  L_c* NÃO coincide com L_B* (L_c* = {L_C:.2f} mm, L_B* = {L_B:.2f} mm).")
if Taxa_calor(L_C) >= q_req:
    print(f"  L_c* ATENDE à exigência q >= {q_req} W (q = {Taxa_calor(L_C):.4f} W).")
else:
    print(f"  L_c* NÃO atende à exigência q >= {q_req} W (q = {Taxa_calor(L_C):.4f} W).")

# COMENTÁRIOS (conceito): no mínimo de c(L), o custo marginal por watt iguala o custo
# médio por watt; no máximo de B(L), iguala o valor atribuído ao watt (c_q = 0,06).
# Como o custo marginal por watt cresce com L e c_q é maior que o custo mínimo por
# watt, L_B* fica à direita de L_c*.

# ------------------------------------------------------------------------------
# QUESTÃO 6(b) - tabela variando o preço do material
# ------------------------------------------------------------------------------
print("\n" + "=" * 90)
print("QUESTÃO 6(b) - TABELA PARA DIFERENTES PREÇOS DO MATERIAL")
print("=" * 90)
print(f"{'c_m (reais/kg)':<15} | {'L_B* (mm)':<12} | {'B(L_B*) (reais)':<16} | "
      f"{'L_c* (mm)':<12} | {'c(L_c*) (reais/W)':<16}")
print("-" * 90)

# >>> CORREÇÃO 13: o laço agora só chama razao_aurea (CORREÇÃO 10). As duas funções
# internas abaixo "embutem" o preço do teste, porque razao_aurea precisa de uma
# função de UMA variável (L).
for cm_teste in [25.0, 45.0, 90.0]:

    def B_atual(L):
        return beneficio(L, cm_teste)

    def c_atual(L):
        return custo_por_watt(L, cm_teste)

    L_B_i, B_i, _, _ = razao_aurea(B_atual, 20.0, 80.0, True)
    L_C_i, c_i, _, _ = razao_aurea(c_atual, 20.0, 80.0, False)
    print(f"{cm_teste:<15.1f} | {L_B_i:<12.3f} | {B_i:<16.6f} | {L_C_i:<12.3f} | {c_i:<16.6f}")
print("-" * 90)

# COMENTÁRIO: quanto mais caro o material, menores ficam L_B* e L_c*, porque cada
# milímetro a mais custa mais enquanto o ganho em q por milímetro é o mesmo.
# Para c_m = 90 a ordem se inverte (L_c* > L_B*): o custo mínimo por watt (0,0678)
# passa a ser maior que o valor do watt c_q (0,06).