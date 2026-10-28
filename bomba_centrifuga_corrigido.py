# -*- coding: utf-8 -*-
# ==============================================================================
# SIMULADO 1 - 1ª PROVA DE CÁLCULO NUMÉRICO APLICADO
# Seleção e operação de uma bomba centrífuga
#
# O que o programa faz (na ordem das questões):
#   Q1  Acha os coeficientes da curva de rendimento resolvendo um sistema linear
#       (eliminação gaussiana com pivoteamento parcial).
#   Q2  Grava o arquivo curva.dat e compara a curva com os ensaios medidos.
#   Q3  Acha o ponto de melhor rendimento (método da razão áurea).
#   Q4  Acha a faixa de operação aceitável (falsa posição e Bairstow).
#   Q5  Acha o ponto de operação da instalação (Newton-Raphson e secante).
#   Q6  Acha a vazão de menor consumo de energia (interpolação quadrática)
#       e repete tudo para vários valores de H0.
#
# Unidades: Q em m3/h | alturas em m | rendimento em fração (0 a 1) dentro das
# contas e em % apenas na hora de imprimir | consumo e em kWh/m3.
# ==============================================================================

import math    # funções matemáticas (raiz quadrada)
import cmath   # raiz quadrada que aceita números complexos (usada no Bairstow)
import sys     # traz o sys.exit(), que encerra o programa mostrando uma mensagem

# ==============================================================================
# DADOS GERAIS E CONSTANTES DO PROBLEMA
# ==============================================================================
H0_inicial = 25.0   # desnível estático [m]
K = 0.0065          # coeficiente de perdas da tubulação [m/(m3/h)^2]
p = 1000.0          # massa específica da água [kg/m3]
g = 9.81            # gravidade [m/s2]
n_min = 0.60        # rendimento mínimo aceitável (fração)

# Ensaios do fabricante, NA ORDEM DA TABELA 1: (Q em m3/h, rendimento em fração)
tabela_1 = [
    (0.0, 0.0),     # Ensaio 1: Q = 0  -> eta = 0,0 %
    (20.0, 0.426),  # Ensaio 2: Q = 20 -> eta = 42,6 %
    (40.0, 0.767),  # Ensaio 3: Q = 40 -> eta = 76,7 %
    (60.0, 0.511),  # Ensaio 4: Q = 60 -> eta = 51,1 %
]

# ==============================================================================
# FUNÇÕES DO MODELO (usadas em várias questões)
# ==============================================================================
# OBS.: rendimento() usa o vetor x calculado na Questão 1. Por isso ela só pode
# ser CHAMADA depois da Questão 1 (definir a função antes não tem problema).

def rendimento(Q):
    """Rendimento eta(Q) em fração. x = [a3, a2, a1, a0] (ordem do sistema da Q1)."""
    return x[0]*Q**3 + x[1]*Q**2 + x[2]*Q + x[3]

def altura_b(Q):
    """Altura fornecida pela bomba [m]."""
    return 52 - 0.0040*(Q**2)

def altura_s(Q, H0_local=25.0):
    """Altura exigida pela instalação [m]. H0_local permite mudar H0 na Questão 6."""
    return H0_local + K*(Q**2)

def consumo_e(Q, H0_local=25.0):
    """Consumo específico de energia [kWh/m3]."""
    H_s = altura_s(Q, H0_local)
    n_Q = rendimento(Q)
    if n_Q <= 0:
        return float('inf')   # rendimento nulo: consumo "infinito" (evita dividir por zero)
    return (p*g*H_s)/(n_Q*3.6e6)

# ==============================================================================
# QUESTÃO 1: COEFICIENTES DA CURVA DE RENDIMENTO
# Eliminação gaussiana com pivoteamento parcial
# ==============================================================================
print("=" * 75)
print("QUESTÃO 1 - COEFICIENTES DA CURVA DE RENDIMENTO")
print("=" * 75)

# --- Montagem do sistema A*x = b a partir da Tabela 1 ---------------------------
# Cada ensaio (Q, eta) vira UMA equação:  a3*Q^3 + a2*Q^2 + a1*Q + a0 = eta
# Então a linha de A é [Q^3, Q^2, Q, 1] e o b dessa linha é o eta medido.
# Montar com um laço evita erro de digitação (como o 6400 no lugar de 64000).
A = []
b = []
for Q_med, eta_med in tabela_1:
    A.append([Q_med**3, Q_med**2, Q_med, 1.0])
    b.append(eta_med)
n = len(b)   # número de equações (= número de incógnitas = 4)

print("Sistema montado (equações na ordem da Tabela 1):")
print("      coef. a3     coef. a2   coef. a1  coef. a0 |  b (eta)")
for i in range(n):
    print(f"  {A[i][0]:10.1f} {A[i][1]:12.1f} {A[i][2]:10.1f} {A[i][3]:8.1f} | {b[i]:.3f}")

# Cópias para calcular o resíduo no final (a eliminação ESTRAGA A e b)
A_orig = [linha[:] for linha in A]
b_orig = b[:]

LIMITE_PIVO = 1e-12   # pivô com valor absoluto menor que isso é tratado como nulo

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

# COMENTÁRIO EXIGIDO: com as equações na ordem da Tabela 1, a primeira linha é a de
# Q = 0, que fica [0, 0, 0, 1]. Então A[0][0] = 0 e a eliminação ingênua tenta
# dividir por zero ao calcular o multiplicador m = A[i][0]/A[0][0] (ZeroDivisionError).
# O pivoteamento parcial resolve trocando essa linha com a de maior pivô (Q = 60).
# O coeficiente a0 sai 0 porque o ensaio Q = 0 mede eta = 0.

# ==============================================================================
# QUESTÃO 2: ARQUIVO DE RESULTADOS E COMPARAÇÃO COM OS ENSAIOS
# ==============================================================================
print("\n" + "=" * 75)
print("QUESTÃO 2 - ARQUIVO curva.dat E DESVIOS")
print("=" * 75)

with open("curva.dat", "w") as arquivo:
    # Cabeçalho: nome de cada coluna com a unidade entre parênteses
    arquivo.write("Q(m3/h)\teta(%)\tHb(m)\tHs(m)\n")
    for Q in range(66):   # Q = 0, 1, ..., 65
        # rendimento(Q) é FRAÇÃO; multiplica por 100 para a coluna ficar em %
        arquivo.write(f"{Q}\t{rendimento(Q)*100:.4f}\t{altura_b(Q):.4f}\t{altura_s(Q):.4f}\n")
print("Arquivo 'curva.dat' gerado com sucesso!")

maior_desvio = 0.0
for Q_medido, n_medido in tabela_1:
    n_calc = rendimento(Q_medido)
    desvio_pp = abs(n_calc - n_medido) * 100    # em pontos percentuais (p.p.)
    if desvio_pp > maior_desvio:
        maior_desvio = desvio_pp
    print(f"Q = {Q_medido:4.1f} m3/h | medido: {n_medido*100:6.2f} % | "
          f"calculado: {n_calc*100:6.2f} % | desvio: {desvio_pp:.2e} p.p.")
print(f"Maior desvio observado: {maior_desvio:.2e} pontos percentuais")

# ==============================================================================
# QUESTÃO 3: PONTO DE MELHOR RENDIMENTO (razão áurea)
# Queremos o MÁXIMO de eta(Q) em [10, 65]
# ==============================================================================
print("\n" + "=" * 75)
print("QUESTÃO 3 - PONTO DE MELHOR RENDIMENTO (Razão Áurea)")
print("=" * 75)

a_ra, b_ra = 10.0, 65.0                    # extremos do intervalo
phi = (math.sqrt(5.0) - 1.0) / 2.0         # razão áurea (0,618...)
Q1_ra = b_ra - phi * (b_ra - a_ra)         # ponto interno da esquerda
Q2_ra = a_ra + phi * (b_ra - a_ra)         # ponto interno da direita
f1_ra = rendimento(Q1_ra)
f2_ra = rendimento(Q2_ra)

iter_ra = 0
aval_ra = 2          # contador de avaliações: já fizemos duas acima

while (b_ra - a_ra) > 1e-6 and iter_ra < 200:
    iter_ra += 1
    if f1_ra > f2_ra:
        # O máximo está à esquerda de Q2_ra: descarta o lado direito
        b_ra = Q2_ra
        Q2_ra = Q1_ra
        f2_ra = f1_ra                      # reaproveita a avaliação (não recalcula)
        Q1_ra = b_ra - phi * (b_ra - a_ra)
        f1_ra = rendimento(Q1_ra)          # só 1 avaliação nova
    else:
        # O máximo está à direita de Q1_ra: descarta o lado esquerdo
        a_ra = Q1_ra
        Q1_ra = Q2_ra
        f1_ra = f2_ra                      # reaproveita a avaliação
        Q2_ra = a_ra + phi * (b_ra - a_ra)
        f2_ra = rendimento(Q2_ra)          # só 1 avaliação nova
    aval_ra += 1

Q_BEP = (a_ra + b_ra) / 2.0
n_BEP = rendimento(Q_BEP)
aval_ra += 1         # convenção: conta também a avaliação final em Q_BEP

print(f"Q_BEP (vazão de melhor rendimento): {Q_BEP:.6f} m3/h")
print(f"Rendimento máximo:                  {n_BEP * 100:.6f} %")
print(f"Iterações: {iter_ra} | Avaliações de eta: {aval_ra}")

# ==============================================================================
# QUESTÃO 4: FAIXA DE OPERAÇÃO ACEITÁVEL (falsa posição e Bairstow)
# Raízes de f(Q) = eta(Q) - 0,60
# ==============================================================================
print("\n" + "=" * 75)
print("QUESTÃO 4 - FAIXA DE OPERAÇÃO (Falsa Posição e Bairstow)")
print("=" * 75)

def func_faixa(Q):
    """f(Q) = eta(Q) - 0,60. A raiz é a vazão em que o rendimento vale 60 %."""
    return rendimento(Q) - n_min

def falsa_posicao(a_fp, b_fp, tol=1e-6, max_iter=50):
    """Falsa posição. Imprime o histórico e devolve (raiz, número de iterações).
    Se não houver mudança de sinal no intervalo, devolve (None, 0)."""
    # Verificação inicial: f(a) e f(b) precisam ter sinais opostos
    if func_faixa(a_fp) * func_faixa(b_fp) > 0:
        print(f"ERRO: não há mudança de sinal em [{a_fp}, {b_fp}].")
        return None, 0

    print(f"  Intervalo [{a_fp}, {b_fp}]")
    print(f"  {'k':>3} | {'Q_l':>10} | {'Q_u':>10} | {'Q_r':>12} | {'f(Q_r)':>12}")
    for k in range(1, max_iter + 1):
        fa = func_faixa(a_fp)
        fb = func_faixa(b_fp)
        Q_r = (a_fp * fb - b_fp * fa) / (fb - fa)   # fórmula da falsa posição
        f_r = func_faixa(Q_r)
        print(f"  {k:3d} | {a_fp:10.6f} | {b_fp:10.6f} | {Q_r:12.8f} | {f_r:12.4e}")
        if abs(f_r) < tol:
            return Q_r, k
        # Mantém o extremo que continua com sinal oposto ao de f(Q_r)
        if fa * f_r < 0:
            b_fp = Q_r
        else:
            a_fp = Q_r
    print("AVISO: número máximo de iterações atingido sem convergir.")
    return Q_r, max_iter

# --- (a) Falsa posição ------------------------------------------------------------
print("(a) Falsa posição - limite inferior da faixa:")
Q_inf, it_inf = falsa_posicao(20.0, 35.0)
print("(a) Falsa posição - limite superior da faixa:")
Q_sup, it_sup = falsa_posicao(45.0, 65.0)
if Q_inf is None or Q_sup is None:
    sys.exit("Não foi possível determinar a faixa de operação. Programa encerrado.")
print(f"Q_inf = {Q_inf:.4f} m3/h ({it_inf} iterações) | Q_sup = {Q_sup:.4f} m3/h ({it_sup} iterações)")

# --- (b) Bairstow --------------------------------------------------------------------
def bairstow(a, r, s, tol=1e-10, max_iter=100):
    """Todas as raízes do polinômio a[0] + a[1]*x + ... + a[n]*x^n.
    a = lista de coeficientes COMEÇANDO PELO TERMO INDEPENDENTE.
    r, s = valores iniciais do fator quadrático x^2 - r*x - s.
    Devolve (lista de raízes, total de iterações)."""
    a = a[:]                    # cópia, para não alterar a lista original
    raizes = []
    n = len(a) - 1              # grau do polinômio
    total_iter = 0

    while n >= 3:               # enquanto sobrar grau >= 3, extrai um fator quadrático
        ri = r
        si = s
        convergiu = False
        for it in range(max_iter):
            total_iter += 1
            b = [0.0] * (n + 1)
            c = [0.0] * (n + 1)
            # Divisão sintética pelo fator quadrático (recorrência dos b)
            b[n] = a[n]
            b[n - 1] = a[n - 1] + ri * b[n]
            for i in range(n - 2, -1, -1):
                b[i] = a[i] + ri * b[i + 1] + si * b[i + 2]
            # Derivadas parciais (recorrência dos c)
            c[n] = b[n]
            c[n - 1] = b[n - 1] + ri * c[n]
            for i in range(n - 2, 0, -1):
                c[i] = b[i] + ri * c[i + 1] + si * c[i + 2]
            # Sistema 2x2 resolvido por Cramer:
            #   c1*dr + c2*ds = -b0
            #   c2*dr + c3*ds = -b1
            det = c[1] * c[3] - c[2] * c[2]
            if det == 0.0:
                sys.exit("Bairstow: determinante nulo. Tente outros r e s iniciais.")
            dr = (-b[0] * c[3] + b[1] * c[2]) / det
            ds = (-b[1] * c[1] + b[0] * c[2]) / det
            ri = ri + dr
            si = si + ds
            if abs(dr) < tol and abs(ds) < tol:
                convergiu = True
                break
        if not convergiu:
            sys.exit("Bairstow não convergiu. Tente outros r e s iniciais.")
        # Raízes do fator quadrático x^2 - ri*x - si = 0 (cmath aceita raiz de negativo)
        d = cmath.sqrt(ri * ri + 4.0 * si)
        raizes.append((ri + d) / 2.0)
        raizes.append((ri - d) / 2.0)
        # O que sobrou da divisão (b2, b3, ..., bn) é o novo polinômio, de grau n-2
        a = b[2:]
        n = n - 2

    # Sobrou grau 2 ou grau 1: resolve direto
    if n == 2:
        d = cmath.sqrt(a[1] ** 2 - 4.0 * a[2] * a[0])
        raizes.append((-a[1] + d) / (2.0 * a[2]))
        raizes.append((-a[1] - d) / (2.0 * a[2]))
    elif n == 1:
        raizes.append(complex(-a[0] / a[1], 0.0))
    return raizes, total_iter

# Polinômio P(Q) = eta(Q) - 0,60 = a3*Q^3 + a2*Q^2 + a1*Q + (a0 - 0,60).
# ATENÇÃO: x está em ordem [a3, a2, a1, a0], mas o Bairstow quer o termo
# independente PRIMEIRO. Por isso a lista abaixo inverte a ordem de x.
coef_poly = [x[3] - n_min, x[2], x[1], x[0]]

raizes_b, it_b = bairstow(coef_poly, 0.1, 0.1)
raizes_b.sort(key=lambda z: z.real)    # organiza da menor para a maior raiz

print(f"\n(b) Bairstow (r0 = s0 = 0,1) - {it_b} iterações. Raízes de P(Q) = eta(Q) - 0,60:")
for z in raizes_b:
    if abs(z.imag) < 1e-9:
        if z.real > 0:
            situacao = "real e positiva -> fisicamente admissível"
        else:
            situacao = "real e NEGATIVA -> sem sentido físico (vazão não pode ser negativa)"
        print(f"  Q = {z.real:10.4f} m3/h  ({situacao})")
    else:
        print(f"  Q = {z.real:.4f} {z.imag:+.4f}i  (complexa -> sem sentido físico)")
print(f"A faixa aceitável é, portanto, [{Q_inf:.2f}; {Q_sup:.2f}] m3/h.")

# ==============================================================================
# QUESTÃO 5: PONTO DE OPERAÇÃO DA INSTALAÇÃO (Newton-Raphson e secante)
# Raiz de f(Q) = Hb(Q) - Hs(Q) com H0 = 25 m
# ==============================================================================
print("\n" + "=" * 75)
print("QUESTÃO 5 - PONTO DE OPERAÇÃO (Newton-Raphson e Secante)")
print("=" * 75)

def func_operacao(Q, H0_local=25.0):
    """f(Q) = Hb(Q) - Hs(Q). No ponto de operação, f = 0."""
    return altura_b(Q) - altura_s(Q, H0_local)

def derivada_operacao(Q):
    """f'(Q) deduzida à mão: Hb' = -2*0,0040*Q  e  Hs' = 2*K*Q."""
    return -2 * 0.0040 * Q - 2 * K * Q

def newton_raphson(f, df, Q0, tol=1e-8, max_iter=100, mostrar=True):
    """Newton-Raphson. f e df são FUNÇÕES passadas como argumento.
    Devolve (raiz, iterações)."""
    Q = Q0
    if mostrar:
        print(f"  {'k':>3} | {'Q_k':>14} | {'f(Q_k)':>12} | {'erro':>12}")
    for k in range(1, max_iter + 1):
        Q_novo = Q - f(Q) / df(Q)
        erro = abs(Q_novo - Q)
        if mostrar:
            print(f"  {k:3d} | {Q_novo:14.8f} | {f(Q_novo):12.4e} | {erro:12.4e}")
        Q = Q_novo
        if erro < tol:
            return Q, k
    print("AVISO: Newton não convergiu no máximo de iterações.")
    return Q, max_iter

def secante(f, Q0, Q1, tol=1e-8, max_iter=100, mostrar=True):
    """Método da secante (usa só f, sem derivada). Devolve (raiz, iterações)."""
    if mostrar:
        print(f"  {'k':>3} | {'Q_k':>14} | {'f(Q_k)':>12} | {'erro':>12}")
    for k in range(1, max_iter + 1):
        f0 = f(Q0)
        f1 = f(Q1)
        Q2 = Q1 - f1 * (Q1 - Q0) / (f1 - f0)
        erro = abs(Q2 - Q1)
        if mostrar:
            print(f"  {k:3d} | {Q2:14.8f} | {f(Q2):12.4e} | {erro:12.4e}")
        Q0 = Q1
        Q1 = Q2
        if erro < tol:
            return Q2, k
    print("AVISO: secante não convergiu no máximo de iterações.")
    return Q1, max_iter

print("(a) Newton-Raphson (Q0 = 30):")
Q_nr, it_nr = newton_raphson(func_operacao, derivada_operacao, 30.0)
print("(b) Secante (Q0 = 30, Q1 = 40):")
Q_sec, it_sec = secante(func_operacao, 30.0, 40.0)

print(f"\nNewton-Raphson: Q_op = {Q_nr:.8f} m3/h em {it_nr} iterações")
print(f"Secante:        Q_op = {Q_sec:.8f} m3/h em {it_sec} iterações")
print(f"Diferença entre os dois métodos: {abs(Q_nr - Q_sec):.2e} m3/h")

Q_op = Q_nr   # usamos o resultado do Newton daqui para frente
print(f"\nVazão de operação (Q_op): {Q_op:.4f} m3/h")
print(f"Altura de operação (H):   {altura_b(Q_op):.4f} m")
print(f"Rendimento eta(Q_op):     {rendimento(Q_op) * 100:.2f} %")
print(f"Consumo e(Q_op):          {consumo_e(Q_op):.6f} kWh/m3")

if Q_inf <= Q_op <= Q_sup:
    print("-> Q_op está DENTRO da faixa aceitável.")
else:
    print("-> Q_op está FORA da faixa aceitável.")

# Posição de Q_op em relação ao Q_BEP (o programa decide, não escrevemos "direita" fixo)
if Q_op > Q_BEP:
    lado = "à direita"
else:
    lado = "à esquerda"
print(f"-> Q_op está {abs(Q_op - Q_BEP):.2f} m3/h {lado} de Q_BEP.")

# ==============================================================================
# QUESTÃO 6: CONSUMO ESPECÍFICO MÍNIMO E ESTUDO PARAMÉTRICO
# ==============================================================================
print("\n" + "=" * 75)
print("QUESTÃO 6 - CONSUMO MÍNIMO (Interpolação Quadrática) E ESTUDO EM H0")
print("=" * 75)

def interpolacao_quadratica(f, x1, x2, x3, tol=1e-6, max_iter=200):
    """Minimiza f pela interpolação quadrática, com a tripla (x1 < x2 < x3).
    Devolve (x_min, f_min, iterações, avaliações, proteções) ou None se a tripla
    inicial não for válida."""
    f1 = f(x1)
    f2 = f(x2)
    f3 = f(x3)
    avaliacoes = 3

    # Verificação da tripla: o ponto do meio precisa ter o MENOR valor de f
    if not (f2 < f1 and f2 < f3):
        print(f"ERRO: tripla inicial inválida: f({x1}) = {f1:.6f}, f({x2}) = {f2:.6f}, f({x3}) = {f3:.6f}")
        return None

    tau = (math.sqrt(5.0) - 1.0) / 2.0
    iteracoes = 0
    protecoes = 0              # quantas vezes foi preciso usar o passo de razão áurea
    x_novo_anterior = None     # guarda o vértice da iteração anterior

    # Para quando o intervalo ficar menor que tol (critério principal)
    while (x3 - x1) > tol and iteracoes < max_iter:

        # Parábola por diferenças divididas: p(x) = A*x^2 + B*x + D
        d1 = (f2 - f1) / (x2 - x1)
        A_par = ((f3 - f2) / (x3 - x2) - d1) / (x3 - x1)
        B_par = d1 - A_par * (x1 + x2)

        # Decide se dá para confiar na parábola
        usar_parabola = True
        if A_par <= 1e-14:                   # A <= 0 ou quase zero: parábola não tem mínimo
            usar_parabola = False
        else:
            x_novo = -B_par / (2.0 * A_par)  # vértice da parábola
            perto_de_algum = min(abs(x_novo - x1), abs(x_novo - x2), abs(x_novo - x3)) < 1e-12
            if x_novo <= x1 or x_novo >= x3 or perto_de_algum:
                usar_parabola = False        # vértice fora do intervalo ou colado num ponto

        if not usar_parabola:
            # Proteção: um passo de razão áurea dentro do maior sub-intervalo
            protecoes += 1
            if (x2 - x1) > (x3 - x2):
                x_novo = x2 - (1.0 - tau) * (x2 - x1)
            else:
                x_novo = x2 + (1.0 - tau) * (x3 - x2)

        f_novo = f(x_novo)       # UMA única avaliação nova por iteração
        avaliacoes += 1
        iteracoes += 1

        # Atualiza a tripla, sempre mantendo o menor valor no meio
        if x_novo < x2:
            if f_novo < f2:
                x3, f3 = x2, f2             # o antigo meio vira o extremo direito
                x2, f2 = x_novo, f_novo     # o novo ponto vira o meio
            else:
                x1, f1 = x_novo, f_novo     # o novo ponto vira o extremo esquerdo
        else:
            if f_novo < f2:
                x1, f1 = x2, f2             # o antigo meio vira o extremo esquerdo
                x2, f2 = x_novo, f_novo     # o novo ponto vira o meio
            else:
                x3, f3 = x_novo, f_novo     # o novo ponto vira o extremo direito

        # Critério secundário: o vértice quase não mudou entre duas iterações
        if x_novo_anterior is not None and abs(x_novo - x_novo_anterior) < tol:
            break
        x_novo_anterior = x_novo

    return x2, f2, iteracoes, avaliacoes, protecoes

# --- (a) H0 = 25 m: consumo_e já usa H0 = 25 por padrão ----------------------------
print("(a) H0 = 25 m, tripla inicial (20, 40, 60):")
resultado = interpolacao_quadratica(consumo_e, 20.0, 40.0, 60.0)
if resultado is None:
    sys.exit("Escolha outra tripla inicial. Programa encerrado.")
Q_min25, e_min25, it_q, aval_q, prot_q = resultado
print(f"  Q* = {Q_min25:.6f} m3/h | e(Q*) = {e_min25:.6f} kWh/m3")
print(f"  Iterações: {it_q} | Avaliações de e: {aval_q} | Proteções acionadas: {prot_q}")

# --- (b) Tabela paramétrica ---------------------------------------------------------
resultados_tabela = []
for H0_atual in [15.0, 25.0, 35.0, 45.0]:

    # Funções de UMA variável (Q) com o H0 deste passo "embutido"
    def f_op(Q):
        return func_operacao(Q, H0_atual)

    def f_consumo(Q):
        return consumo_e(Q, H0_atual)

    Q_op_i, _ = newton_raphson(f_op, derivada_operacao, 30.0, mostrar=False)
    resultado = interpolacao_quadratica(f_consumo, 20.0, 40.0, 60.0)
    if resultado is None:
        sys.exit(f"Tripla inválida para H0 = {H0_atual}. Programa encerrado.")
    Q_est, e_est, _, _, _ = resultado
    resultados_tabela.append((H0_atual, Q_op_i, rendimento(Q_op_i) * 100, Q_est, e_est))

print("\n(b) Tabela paramétrica:")
print(f"| {'H0 (m)':^6} | {'Q_op (m3/h)':^12} | {'eta(Q_op) (%)':^13} | {'Q* (m3/h)':^10} | {'e(Q*) (kWh/m3)':^15} |")
print("-" * 71)
for H0_i, Qop_i, eta_i, Qe_i, ee_i in resultados_tabela:
    print(f"| {H0_i:^6.0f} | {Qop_i:^12.4f} | {eta_i:^13.2f} | {Qe_i:^10.4f} | {ee_i:^15.6f} |")
print("-" * 71)

# --- Comentários pedidos: as respostas saem dos CÁLCULOS, não de texto fixo ----------
print("\nComentário 1 - Q* é maior ou menor que Q_BEP?")
for H0_i, Qop_i, eta_i, Qe_i, ee_i in resultados_tabela:
    if Qe_i < Q_BEP:
        relacao = "MENOR"
    else:
        relacao = "MAIOR"
    print(f"  H0 = {H0_i:.0f} m: Q* = {Qe_i:.2f} é {relacao} que Q_BEP = {Q_BEP:.2f}")
# Explicação: Hs cresce com Q. Minimizar Hs/eta pede um ponto onde eta ainda está
# subindo (derivada positiva), isto é, à ESQUERDA do máximo de eta (Q_BEP).

print("\nComentário 2 - Para quais H0 o ponto de operação cai fora da faixa aceitável?")
for H0_i, Qop_i, eta_i, Qe_i, ee_i in resultados_tabela:
    if Qop_i < Q_inf or Qop_i > Q_sup:
        print(f"  H0 = {H0_i:.0f} m: Q_op = {Qop_i:.2f} m3/h está FORA de [{Q_inf:.2f}; {Q_sup:.2f}] "
              f"(eta = {eta_i:.2f} %)")
