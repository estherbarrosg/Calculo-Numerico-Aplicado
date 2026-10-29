# -*- coding: utf-8 -*-
# ==============================================================================
# SIMULADO 2 - 1ª PROVA DE CÁLCULO NUMÉRICO APLICADO
# Aleta de resfriamento de um componente eletrônico
#
#   Q1  Perfil de temperatura: sistema tridiagonal (Thomas e Gauss-Seidel)
#   Q2  Arquivo perfil.dat e validação com a solução analítica
#   Q3  Comprimento mínimo para a ponta fria (falsa posição e Newton-Raphson)
#   Q4  Espessura ótima (razão áurea)
#   Q5  Comprimento economicamente ótimo (interpolação quadrática x razão áurea)
#   Q6  Estudo paramétrico em h
#
# Unidades: comprimentos em m (só imprimimos em mm quando o enunciado pede)
#           temperaturas em °C | theta = T - T_inf em K | q em W | custo em reais
# ==============================================================================

import math
import sys

# ==============================================================================
# DADOS DO PROBLEMA (Tabela 1 do enunciado)
# ==============================================================================
k_al = 200.0      # condutividade térmica do alumínio [W/(m.K)]
rho = 2700.0      # massa específica do alumínio [kg/m3]
w = 0.030         # largura da aleta [m]
t0 = 0.0008       # espessura da aleta do caso base [m]
L0 = 0.040        # comprimento da aleta do caso base [m]
h0 = 50.0         # coeficiente de convecção do caso base [W/(m2.K)]
T_inf = 25.0      # temperatura do ar [°C]
T_b = 85.0        # temperatura da base [°C] (prescrita nas Questões 1 a 4)
theta_b = T_b - T_inf   # excesso de temperatura na base [K]
N = 20            # número de intervalos da malha
T_seg = 50.0      # temperatura máxima admissível na ponta [°C]
A_p = t0 * L0     # área de perfil FIXA da Questão 4 [m2] (= 3,2e-5)

# Dados da Questão 5 (dissipador completo)
N_a = 12          # número de aletas
P_comp = 40.0     # potência dissipada pelo componente [W]
A_exp = 6.0e-4    # área de base exposta [m2]
c_m = 40.0        # preço do material [reais por kg]
c_T = 0.15        # custo térmico [reais por K]

# ==============================================================================
# FUNÇÕES DO MODELO ANALÍTICO (equações 1 a 4 do enunciado)
# ==============================================================================
def parametros_aleta(h, t):
    """Devolve (P, Ac, m): perímetro [m], área da seção [m2] e m [1/m]."""
    P = 2.0 * (w + t)
    Ac = w * t
    m = math.sqrt(h * P / (k_al * Ac))
    return P, Ac, m

def theta_exata(x, L, h, t):
    """Excesso de temperatura theta(x) [K], equação (3). x e L em m."""
    P, Ac, m = parametros_aleta(h, t)
    r = h / (m * k_al)                       # razão h/(m.k), sem unidade
    numerador = math.cosh(m * (L - x)) + r * math.sinh(m * (L - x))
    denominador = math.cosh(m * L) + r * math.sinh(m * L)
    return theta_b * numerador / denominador

def q_base_exata(L, h, t):
    """Taxa de calor que entra pela base da aleta [W], equação (4)."""
    P, Ac, m = parametros_aleta(h, t)
    r = h / (m * k_al)
    numerador = math.sinh(m * L) + r * math.cosh(m * L)
    denominador = math.cosh(m * L) + r * math.sinh(m * L)
    return math.sqrt(h * P * k_al * Ac) * theta_b * numerador / denominador

def T_ponta_exata(L, h, t):
    """Temperatura na ponta (x = L) [°C]: equação (3) com x = L."""
    P, Ac, m = parametros_aleta(h, t)
    r = h / (m * k_al)
    return T_inf + theta_b / (math.cosh(m * L) + r * math.sinh(m * L))

# ==============================================================================
# SISTEMA TRIDIAGONAL: MONTAGEM, THOMAS E GAUSS-SEIDEL
# ==============================================================================
def monta_sistema(N, L, h, t):
    """Monta os vetores e (subdiagonal), f (diagonal), g (superdiagonal) e r (lado
    direito) das equações (5) e (6). Incógnitas: theta_1 ... theta_N.
    ATENÇÃO AOS ÍNDICES: o nó i fica na posição p = i - 1 (nó 1 -> posição 0)."""
    P, Ac, m = parametros_aleta(h, t)
    dx = L / N
    beta = (m * dx) ** 2
    e = [0.0] * N
    f = [0.0] * N
    g = [0.0] * N
    r = [0.0] * N
    for p in range(N - 1):          # nós internos: -theta(i-1) + (2+beta)*theta(i) - theta(i+1) = 0
        f[p] = 2.0 + beta
        g[p] = -1.0
    for p in range(1, N):
        e[p] = -1.0
    f[N - 1] = 1.0 + h * dx / k_al + beta / 2.0     # nó da ponta (equação 6)
    r[0] = theta_b                                  # theta_0 conhecido vai para o lado direito
    return e, f, g, r, dx

def thomas(e, f, g, r):
    """Algoritmo de Thomas. Trabalha em cópias e devolve a solução."""
    n = len(f)
    e = e[:]
    f = f[:]
    g = g[:]
    r = r[:]
    for p in range(1, n):                       # decomposição LU
        e[p] = e[p] / f[p - 1]
        f[p] = f[p] - e[p] * g[p - 1]
    for p in range(1, n):                       # substituição progressiva
        r[p] = r[p] - e[p] * r[p - 1]
    solucao = [0.0] * n
    solucao[n - 1] = r[n - 1] / f[n - 1]        # substituição regressiva
    for p in range(n - 2, -1, -1):
        solucao[p] = (r[p] - g[p] * solucao[p + 1]) / f[p]
    return solucao

def gauss_seidel(e, f, g, r, tol=1e-8, max_iter=100000):
    """Gauss-Seidel para o sistema tridiagonal, com theta inicial nulo.
    Devolve (solução, iterações). Critério: maior variação < tol."""
    n = len(f)
    theta = [0.0] * n
    for it in range(1, max_iter + 1):
        maior_variacao = 0.0
        for p in range(n):
            soma = r[p]
            if p > 0:
                soma -= e[p] * theta[p - 1]     # usa o valor JÁ atualizado
            if p < n - 1:
                soma -= g[p] * theta[p + 1]     # usa o valor da iteração anterior
            novo = soma / f[p]
            variacao = abs(novo - theta[p])
            if variacao > maior_variacao:
                maior_variacao = variacao
            theta[p] = novo
        if maior_variacao < tol:
            return theta, it
    print("AVISO: Gauss-Seidel não convergiu no máximo de iterações.")
    return theta, max_iter

# ==============================================================================
# MÉTODOS NUMÉRICOS GENÉRICOS (recebem a função como argumento)
# ==============================================================================
def falsa_posicao(func, a, b, tol=1e-6, max_iter=200, mostrar=True):
    """Falsa posição. Devolve (raiz, iterações) ou (None, 0) se não houver mudança de sinal."""
    if func(a) * func(b) > 0:
        print(f"ERRO: f não muda de sinal em [{a}, {b}].")
        return None, 0
    if mostrar:
        print(f"  {'k':>3} | {'L_l':>10} | {'L_u':>10} | {'L_r':>12} | {'f(L_r)':>12}")
    for k in range(1, max_iter + 1):
        fa = func(a)
        fb = func(b)
        L_r = b - fb * (b - a) / (fb - fa)
        f_r = func(L_r)
        if mostrar:
            print(f"  {k:3d} | {a:10.6f} | {b:10.6f} | {L_r:12.8f} | {f_r:12.4e}")
        if abs(f_r) < tol:
            return L_r, k
        if fa * f_r < 0:
            b = L_r
        else:
            a = L_r
    print("AVISO: falsa posição atingiu o máximo de iterações.")
    return L_r, max_iter

def newton_raphson(func, dfunc, x0, tol=1e-10, max_iter=100, mostrar=True):
    """Newton-Raphson. Devolve (raiz, iterações)."""
    x = x0
    if mostrar:
        print(f"  {'k':>3} | {'L_k':>12} | {'f(L_k)':>12} | {'dL':>12}")
    for k in range(1, max_iter + 1):
        x_novo = x - func(x) / dfunc(x)
        erro = abs(x_novo - x)
        if mostrar:
            print(f"  {k:3d} | {x_novo:12.8f} | {func(x_novo):12.4e} | {erro:12.4e}")
        x = x_novo
        if erro < tol:
            return x, k
    print("AVISO: Newton atingiu o máximo de iterações.")
    return x, max_iter

def razao_aurea(func, a, b, maximizar, tol=1e-6, max_iter=200):
    """Razão áurea. Devolve (x_otimo, iterações, avaliações).
    Convenção: contamos 2 avaliações iniciais + 1 por iteração."""
    phi = (math.sqrt(5.0) - 1.0) / 2.0
    if maximizar:
        sinal = -1.0         # maximizar f = minimizar -f
    else:
        sinal = 1.0
    x1 = b - phi * (b - a)
    x2 = a + phi * (b - a)
    f1 = sinal * func(x1)
    f2 = sinal * func(x2)
    avaliacoes = 2
    iteracoes = 0
    while (b - a) > tol and iteracoes < max_iter:
        iteracoes += 1
        if f1 < f2:          # (objetivo é sempre MINIMIZAR; o sinal já tratou o máximo)
            b = x2
            x2 = x1
            f2 = f1
            x1 = b - phi * (b - a)
            f1 = sinal * func(x1)
        else:
            a = x1
            x1 = x2
            f1 = f2
            x2 = a + phi * (b - a)
            f2 = sinal * func(x2)
        avaliacoes += 1
    return (a + b) / 2.0, iteracoes, avaliacoes

def tripla_valida(func, x1, x2, x3):
    """O ponto do meio precisa ter o MENOR valor de func."""
    return func(x2) < func(x1) and func(x2) < func(x3)

def interpolacao_quadratica(func, x1, x2, x3, tol=1e-6, max_iter=200):
    """Minimiza func por interpolação quadrática. Devolve
    (x_min, f_min, iterações, avaliações, proteções) ou None se a tripla for inválida."""
    f1 = func(x1)
    f2 = func(x2)
    f3 = func(x3)
    avaliacoes = 3
    if not (f2 < f1 and f2 < f3):
        print("ERRO: tripla inicial inválida.")
        return None
    tau = (math.sqrt(5.0) - 1.0) / 2.0
    iteracoes = 0
    protecoes = 0
    x_novo_anterior = None
    while (x3 - x1) > tol and iteracoes < max_iter:
        # parábola por diferenças divididas
        d1 = (f2 - f1) / (x2 - x1)
        A_par = ((f3 - f2) / (x3 - x2) - d1) / (x3 - x1)
        B_par = d1 - A_par * (x1 + x2)
        usar_parabola = True
        if A_par <= 1e-14:
            usar_parabola = False
        else:
            x_novo = -B_par / (2.0 * A_par)
            perto = min(abs(x_novo - x1), abs(x_novo - x2), abs(x_novo - x3)) < 1e-12
            if x_novo <= x1 or x_novo >= x3 or perto:
                usar_parabola = False
        if not usar_parabola:                 # proteção: passo de razão áurea
            protecoes += 1
            if (x2 - x1) > (x3 - x2):
                x_novo = x2 - (1.0 - tau) * (x2 - x1)
            else:
                x_novo = x2 + (1.0 - tau) * (x3 - x2)
        f_novo = func(x_novo)
        avaliacoes += 1
        iteracoes += 1
        if x_novo < x2:                       # atualiza a tripla
            if f_novo < f2:
                x3, f3 = x2, f2
                x2, f2 = x_novo, f_novo
            else:
                x1, f1 = x_novo, f_novo
        else:
            if f_novo < f2:
                x1, f1 = x2, f2
                x2, f2 = x_novo, f_novo
            else:
                x3, f3 = x_novo, f_novo
        if x_novo_anterior is not None and abs(x_novo - x_novo_anterior) < tol:
            break                              # critério secundário: vértice parou de mudar
        x_novo_anterior = x_novo
    return x2, f2, iteracoes, avaliacoes, protecoes

# ==============================================================================
# QUESTÃO 1 - PERFIL DE TEMPERATURA: SISTEMA LINEAR
# ==============================================================================
print("=" * 75)
print("QUESTÃO 1 - PERFIL DE TEMPERATURA (Thomas e Gauss-Seidel)")
print("=" * 75)

e, f, g, r, dx = monta_sistema(N, L0, h0, t0)

# (a) Thomas
theta_thomas = thomas(e, f, g, r)

# Resíduo ||r||_inf = max|A*theta - b| usando os vetores ORIGINAIS e, f, g, r
maior_residuo = 0.0
for p in range(N):
    soma = f[p] * theta_thomas[p]
    if p > 0:
        soma += e[p] * theta_thomas[p - 1]
    if p < N - 1:
        soma += g[p] * theta_thomas[p + 1]
    residuo_p = abs(soma - r[p])
    if residuo_p > maior_residuo:
        maior_residuo = residuo_p

# Vetor de temperaturas dos nós 0..N em °C (nó 0 = base)
T_num = [T_b]
for p in range(N):
    T_num.append(T_inf + theta_thomas[p])

print("(a) Thomas:")
print(f"  T(L/2) = {T_num[N // 2]:.4f} °C   (nó i = N/2 = {N // 2})")
print(f"  T(L)   = {T_num[N]:.4f} °C")
print(f"  Resíduo ||r||_inf = {maior_residuo:.3e}")
# COMENTÁRIO: o Thomas faz só 3 laços simples (custo proporcional a N), enquanto a
# eliminação gaussiana geral faz laços aninhados (custo proporcional a N^3) e ainda
# ignora que quase todos os elementos da matriz são zero.

# (b) Gauss-Seidel
theta_gs, it_gs = gauss_seidel(e, f, g, r)
maior_dif = 0.0
for p in range(N):
    dif = abs(theta_gs[p] - theta_thomas[p])
    if dif > maior_dif:
        maior_dif = dif
print("(b) Gauss-Seidel (theta inicial = 0, tol = 1e-8 K):")
print(f"  Iterações: {it_gs} | Maior diferença em relação a Thomas: {maior_dif:.2e} K")
# COMENTÁRIO: converge porque a matriz tem dominância diagonal estrita: |2+beta| > 2
# nos nós internos e |1 + h*dx/k + beta/2| > 1 na ponta. Thomas resolve de uma vez
# (direto); Gauss-Seidel precisa de ~2000 varreduras sobre os N nós.

# ==============================================================================
# QUESTÃO 2 - ARQUIVO perfil.dat E VALIDAÇÃO
# ==============================================================================
print("\n" + "=" * 75)
print("QUESTÃO 2 - ARQUIVO perfil.dat E VALIDAÇÃO")
print("=" * 75)

maior_erro = 0.0
no_maior_erro = 0
with open("perfil.dat", "w") as arquivo:
    arquivo.write("x(mm)\tT_num(C)\tT_exata(C)\terro(K)\n")
    for i in range(N + 1):
        x_i = i * dx                                           # posição do nó i [m]
        T_ex = T_inf + theta_exata(x_i, L0, h0, t0)            # equação (3)
        erro = abs(T_num[i] - T_ex)
        arquivo.write(f"{x_i * 1000:.4f}\t{T_num[i]:.6f}\t{T_ex:.6f}\t{erro:.3e}\n")
        if erro > maior_erro:
            maior_erro = erro
            no_maior_erro = i
print("Arquivo 'perfil.dat' gerado.")
print(f"Maior erro: {maior_erro:.3e} K no nó {no_maior_erro} (x = {no_maior_erro * dx * 1000:.1f} mm)")
print(f"T(L): numérica = {T_num[N]:.4f} °C | exata = {T_ponta_exata(L0, h0, t0):.4f} °C")

P0, Ac0, m0 = parametros_aleta(h0, t0)
q_exata = q_base_exata(L0, h0, t0)
q_1a = k_al * Ac0 * (theta_b - theta_thomas[0]) / dx                                   # 1ª ordem
q_2a = -k_al * Ac0 * (-3.0 * theta_b + 4.0 * theta_thomas[0] - theta_thomas[1]) / (2.0 * dx)  # 2ª ordem
print(f"q_b exata = {q_exata:.4f} W")
print(f"q_b numérica (1ª ordem) = {q_1a:.4f} W | erro relativo = {(q_1a - q_exata) / q_exata * 100:.2f} %")
print(f"q_b numérica (2ª ordem) = {q_2a:.4f} W | erro relativo = {(q_2a - q_exata) / q_exata * 100:.2f} %")
# COMENTÁRIO: a derivada de 1ª ordem na base tem erro da ordem de dx, e a temperatura
# tem erro da ordem de dx^2. Por isso o erro de q_b é bem maior que o da temperatura.

# ==============================================================================
# QUESTÃO 3 - COMPRIMENTO MÍNIMO PARA A PONTA FRIA
# f(L) = T_ponta(L) - 50
# ==============================================================================
print("\n" + "=" * 75)
print("QUESTÃO 3 - COMPRIMENTO MÍNIMO (T_ponta = 50 °C)")
print("=" * 75)

def f_ponta(L):
    return T_ponta_exata(L, h0, t0) - T_seg

def df_ponta(L):
    """Derivada de f. Com D(L) = cosh(mL) + (h/mk)*sinh(mL), T_ponta = T_inf + theta_b/D, então
    f'(L) = -theta_b * D'(L) / D(L)^2, com D'(L) = m*sinh(mL) + (h/k)*cosh(mL).
    O sinal é NEGATIVO: quanto maior L, mais fria a ponta."""
    P, Ac, m = parametros_aleta(h0, t0)
    r = h0 / (m * k_al)
    D = math.cosh(m * L) + r * math.sinh(m * L)
    dD = m * math.sinh(m * L) + (h0 / k_al) * math.cosh(m * L)
    return -theta_b * dD / D ** 2

print("(a) Falsa posição em [0,02; 0,10] m:")
L_s, it_fp = falsa_posicao(f_ponta, 0.02, 0.10, tol=1e-6)
if L_s is None:
    sys.exit("Sem mudança de sinal. Programa encerrado.")
print(f"  L_s = {L_s * 1000:.4f} mm em {it_fp} iterações")

print("(b) Newton-Raphson com L0 = 0,04 m:")
L_s_nr, it_nr = newton_raphson(f_ponta, df_ponta, 0.04, tol=1e-10)
print(f"  L_s = {L_s_nr * 1000:.4f} mm em {it_nr} iterações")
# COMENTÁRIO: Newton converge quadraticamente (poucas iterações) porque usa a derivada;
# a falsa posição é linear e mantém sempre o mesmo extremo preso, por isso precisa de mais.

# ==============================================================================
# QUESTÃO 4 - ESPESSURA ÓTIMA (máximo de q_b com área de perfil fixa)
# ==============================================================================
print("\n" + "=" * 75)
print("QUESTÃO 4 - ESPESSURA ÓTIMA (Razão Áurea)")
print("=" * 75)

def q_da_espessura(t):
    """q_b(t) com h = 50. O comprimento muda com t: L(t) = A_p / t."""
    return q_base_exata(A_p / t, h0, t)

t_otimo, it_t, aval_t = razao_aurea(q_da_espessura, 0.0002, 0.0030, True, tol=1e-7)
L_otimo = A_p / t_otimo
q_otimo = q_da_espessura(t_otimo)
P_o, Ac_o, m_o = parametros_aleta(h0, t_otimo)
q_base_caso = q_base_exata(L0, h0, t0)
print(f"t* = {t_otimo * 1000:.4f} mm | L* = {L_otimo * 1000:.2f} mm")
print(f"q_b(t*) = {q_otimo:.4f} W | mL* = {m_o * L_otimo:.3f}")
print(f"Iterações: {it_t} | Avaliações: {aval_t}")
print(f"Ganho sobre a aleta base ({q_base_caso:.4f} W): {(q_otimo - q_base_caso) / q_base_caso * 100:.2f} %")
# COMENTÁRIO: aleta fina e longa tem pouca seção para conduzir calor até a ponta;
# aleta grossa e curta tem pouca superfície para trocar calor. O máximo está no meio.

# ==============================================================================
# QUESTÃO 5 - COMPRIMENTO ECONOMICAMENTE ÓTIMO
# ==============================================================================
print("\n" + "=" * 75)
print("QUESTÃO 5 - COMPRIMENTO ECONOMICAMENTE ÓTIMO")
print("=" * 75)

def theta_base_balanco(L, h):
    """theta_b(L) [K] pelo balanço P_comp = N_a*q_b + h*A_exp*theta_b."""
    P, Ac, m = parametros_aleta(h, t0)
    r = h / (m * k_al)
    G = (math.sinh(m * L) + r * math.cosh(m * L)) / (math.cosh(m * L) + r * math.sinh(m * L))
    S = math.sqrt(h * P * k_al * Ac)
    return P_comp / (N_a * S * G + h * A_exp)

def custo_total(L, h):
    """C(L) em reais: custo do material + custo térmico."""
    material = c_m * rho * N_a * w * t0 * L
    termico = c_T * theta_base_balanco(L, h)
    return material + termico

def custo_h0(L):
    return custo_total(L, h0)

print("(a) Verificação da tripla (0,03; 0,05; 0,09) m:")
valores = [custo_h0(0.03), custo_h0(0.05), custo_h0(0.09)]
print(f"  C = {valores[0]:.4f} | {valores[1]:.4f} | {valores[2]:.4f}  -> válida: {tripla_valida(custo_h0, 0.03, 0.05, 0.09)}")

print("(b) Interpolação quadrática:")
res = interpolacao_quadratica(custo_h0, 0.03, 0.05, 0.09, tol=1e-6)
if res is None:
    sys.exit("Tripla inválida. Programa encerrado.")
L_c, C_min, it_q, aval_q, prot_q = res
print(f"  L_c* = {L_c * 1000:.4f} mm | C(L_c*) = {C_min:.4f} reais | T_b(L_c*) = {T_inf + theta_base_balanco(L_c, h0):.2f} °C")
print(f"  Iterações: {it_q} | Avaliações: {aval_q} | Proteções: {prot_q}")

print("(c) Razão áurea (minimização):")
L_c_ra, it_ra, aval_ra = razao_aurea(custo_h0, 0.01, 0.12, False, tol=1e-6)
print(f"  L_c* = {L_c_ra * 1000:.4f} mm | Iterações: {it_ra} | Avaliações: {aval_ra}")
# COMENTÁRIO: a interpolação quadrática usou bem menos avaliações que a razão áurea, pois
# a função é suave e a parábola aproxima bem o mínimo.

# ==============================================================================
# QUESTÃO 6 - ESTUDO PARAMÉTRICO EM h
# ==============================================================================
print("\n" + "=" * 75)
print("QUESTÃO 6 - ESTUDO PARAMÉTRICO EM h")
print("=" * 75)

def tripla_automatica(func):
    """Varre L = 0,01; 0,02; ...; 0,12 m e devolve o menor valor com seus vizinhos."""
    pontos = []
    for i in range(1, 13):
        pontos.append(i * 0.01)
    melhor = 0
    for i in range(len(pontos)):
        if func(pontos[i]) < func(pontos[melhor]):
            melhor = i
    if melhor == 0:
        melhor = 1
    if melhor == len(pontos) - 1:
        melhor = len(pontos) - 2
    return pontos[melhor - 1], pontos[melhor], pontos[melhor + 1]

resultados = []
for h_atual in [25.0, 50.0, 100.0, 200.0]:

    def f_ponta_h(L):
        return T_ponta_exata(L, h_atual, t0) - T_seg

    def q_da_espessura_h(t):
        return q_base_exata(A_p / t, h_atual, t)

    def custo_h(L):
        return custo_total(L, h_atual)

    L_s_h, _ = falsa_posicao(f_ponta_h, 0.005, 0.20, tol=1e-6, mostrar=False)
    t_h, _, _ = razao_aurea(q_da_espessura_h, 0.0002, 0.0030, True, tol=1e-7)
    q_h = q_da_espessura_h(t_h)
    P_h, Ac_h, m_h = parametros_aleta(h_atual, t_h)
    mL_h = m_h * A_p / t_h
    tripla_base_ok = tripla_valida(custo_h, 0.03, 0.05, 0.09)
    tripla = tripla_automatica(custo_h)
    resultado_c = interpolacao_quadratica(custo_h, tripla[0], tripla[1], tripla[2], tol=1e-6)
    L_c_h = resultado_c[0]
    resultados.append((h_atual, L_s_h, t_h, q_h, L_c_h, mL_h, tripla_base_ok, tripla))

print(f"| {'h (W/m2K)':^10} | {'L_s (mm)':^9} | {'t* (mm)':^8} | {'q_b(t*) (W)':^11} | {'L_c* (mm)':^9} |")
print("-" * 62)
for h_i, Ls_i, t_i, q_i, Lc_i, mL_i, ok_i, tr_i in resultados:
    print(f"| {h_i:^10.0f} | {Ls_i * 1000:^9.2f} | {t_i * 1000:^8.4f} | {q_i:^11.3f} | {Lc_i * 1000:^9.2f} |")
print("-" * 62)

print("\nComentários (calculados pelo programa):")
print("(i) A tripla (0,03; 0,05; 0,09) é válida para h =", end="")
for h_i, Ls_i, t_i, q_i, Lc_i, mL_i, ok_i, tr_i in resultados:
    if ok_i:
        print(f" {h_i:.0f}", end="")
print(" e INVÁLIDA para os demais.")
print("    Triplas automáticas usadas:", [tuple(round(v, 2) for v in r_[7]) for r_ in resultados])
print("(iii) mL* em cada caso:", [round(r_[5], 3) for r_ in resultados])
# COMENTÁRIO (ii): quanto maior h, maior m (m^2 = hP/(kAc)) e mais rápido theta cai ao longo
# da aleta, então L_s e L_c* diminuem. O q_b ótimo aumenta, e a espessura ótima aumenta
# porque o ótimo mantém mL perto de 1,4 e m cresce com raiz de h/t.