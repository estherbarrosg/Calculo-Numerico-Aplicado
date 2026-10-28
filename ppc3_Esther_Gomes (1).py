"""
PPC3 - Cálculo Numérico Aplicado (Prof. Rafael Gabler Gontijo)

Condução de calor transiente 1D em uma placa de combustível nuclear (UO2)
com geração interna de energia.

Equação governante (meia placa, 0 <= x <= L):
    dT/dt = alpha * d2T/dx2 + q / (rho * cp)

Condições de contorno:
    x = 0 (centro) : simetria adiabática (dT/dx = 0)
    x = L (face)   : convecção com o fluido (h, T_inf)

Método numérico:
    - diferenças finitas na malha nodal (nó 0 = centro, nó N-1 = superfície)
    - esquema IMPLÍCITO no tempo (estável para qualquer passo de tempo)
    - sistema tridiagonal resolvido com o algoritmo de Thomas (decomposição L.U)

Validações:
    1) q = 0  -> comparação com a solução exata por séries
    2) q > 0  -> comparação com a solução analítica do regime estacionário
    3) estudo de refinamento da malha e do passo de tempo

Organização do arquivo:
    1. Parâmetros do problema
    2. Algoritmo de Thomas
    3. Sistema linear e evolução no tempo
    4. Soluções analíticas (usadas só para validar)
    5. Validação
    6. Gráficos
    7. Programa principal

Execução:  python ppc3.py
Bibliotecas: math (funções matemáticas) e matplotlib/numpy (somente gráficos).
O método numérico em si é todo implementado aqui, com listas, laços e funções.
"""

import math
import os

import numpy as np                # usado apenas para montar os gráficos
import matplotlib.pyplot as plt


# ==============================================================================
# 1. PARÂMETROS DO PROBLEMA
# ==============================================================================
# Estas variáveis ficam fora das funções (são "globais"). As funções podem LER
# esses valores diretamente. O que muda de um teste para outro (q, temperatura
# inicial, malha, passo de tempo) é passado como argumento da função.

# --- Material: combustível cerâmico UO2 (valores típicos de reator a água) ---
k_cond = 4.0                 # condutividade térmica [W/(m.K)]
rho = 10500.0                # massa específica [kg/m3]
cp = 300.0                   # calor específico [J/(kg.K)]
q_ponto = 2.0e8              # taxa de geração interna [W/m3] (200 MW/m3)
alpha = k_cond / (rho * cp)  # difusividade térmica [m2/s]

# --- Geometria e contorno ---
L = 0.005                    # meia espessura da placa [m] (5 mm)
h = 25000.0                  # coeficiente de convecção da água [W/(m2.K)]
T_inf = 300.0                # temperatura do fluido [°C]
T_inicial = 300.0            # temperatura inicial uniforme [°C]
Bi_L = h * L / k_cond        # Biot da PLACA (usado na solução exata)

# --- Discretização (caso principal) ---
N = 21                       # número de nós da malha
dt = 0.1                     # passo de tempo [s]
tempo_total = 60.0           # tempo físico simulado [s]

# --- Saída e solução exata ---
N_TERMOS = 100               # termos da série da solução exata
MOSTRAR_FIGURAS = True       # False: só salva os PNGs, sem abrir janelas
PASTA_RESULTADOS = "resultados"


def verificar_entradas(n_nos, passo_dt, t_final, q):
    """Interrompe o programa (raise) se algum dado for inválido."""
    for nome, valor in (("k_cond", k_cond), ("rho", rho), ("cp", cp),
                        ("L", L), ("h", h)):
        if valor <= 0:
            raise ValueError(f"O parâmetro '{nome}' deve ser positivo.")
    if n_nos < 3 or int(n_nos) != n_nos:
        raise ValueError("O número de nós deve ser um inteiro >= 3.")
    if passo_dt <= 0:
        raise ValueError("O passo de tempo dt deve ser positivo.")
    if t_final < passo_dt:
        raise ValueError("O tempo total não pode ser menor que dt.")
    if q < 0:
        raise ValueError("A geração interna q não pode ser negativa.")


def numeros_da_malha(n_nos, passo_dt):
    """Calcula dx, Fourier da malha e Biot da malha.

    Atenção: o Biot que aparece na matriz [A] é o da MALHA (h*dx/k);
    o Biot da solução exata é o da PLACA (Bi_L = h*L/k). Não confundir.
    """
    dx = L / (n_nos - 1)
    Fo = alpha * passo_dt / dx ** 2     # Fo = alpha*dt/dx^2
    Bi_malha = h * dx / k_cond
    return dx, Fo, Bi_malha


# ==============================================================================
# 2. ALGORITMO DE THOMAS (decomposição L.U para matrizes tridiagonais)
# ==============================================================================
# Notação do enunciado: e = diagonal inferior, f = diagonal principal,
# g = diagonal superior, todos com n posições (e[0] e g[n-1] não são usados).
# A matriz [A] deste problema não muda no tempo. Por isso a decomposição L.U é
# feita uma única vez (thomas_fatorar) e, a cada passo de tempo, só repetimos as
# substituições progressiva e regressiva (thomas_resolver).

def thomas_fatorar(e, f, g):
    """Decomposição L.U

    Devolve dois vetores: e_lu (fatores de [L]) e f_lu (diagonal de [U]).
    Os vetores originais não são alterados.
    """
    n = len(f)
    if len(e) != n or len(g) != n:
        raise ValueError("e, f e g devem ter o mesmo tamanho.")

    e = e[:]      # cópias, para não alterar os dados de fora da função
    f = f[:]

    for k in range(1, n):
        if f[k - 1] == 0.0:
            raise ZeroDivisionError(f"Pivô nulo na linha {k - 1}.")
        e[k] = e[k] / f[k - 1]
        f[k] = f[k] - e[k] * g[k - 1]
    if f[n - 1] == 0.0:
        raise ZeroDivisionError("Pivô nulo na última linha: matriz singular.")
    return e, f


def thomas_resolver(e_lu, f_lu, g, r):
    """Resolve [L][U]{x} = {r} (2º e 3º laços do pseudocódigo).

    e_lu, f_lu: resultado de thomas_fatorar; g: diagonal superior original;
    r: termos independentes (o vetor original não é alterado).
    """
    n = len(f_lu)
    if len(r) != n:
        raise ValueError("O vetor r deve ter o mesmo tamanho da matriz.")
    r = r[:]
    x = [0.0] * n

    # Substituição progressiva: [L]{d} = {r}   (d fica guardado em r)
    for k in range(1, n):
        r[k] = r[k] - e_lu[k] * r[k - 1]

    # Substituição regressiva: [U]{x} = {d}
    x[n - 1] = r[n - 1] / f_lu[n - 1]
    for k in range(n - 2, -1, -1):
        x[k] = (r[k] - g[k] * x[k + 1]) / f_lu[k]
    return x


# ==============================================================================
# 3. SISTEMA LINEAR E EVOLUÇÃO NO TEMPO (esquema implícito)
# ==============================================================================

def montar_matriz(n_nos, Fo, Bi):
    """Diagonais e, f, g da matriz [A] (equações 69, 70 e 71 do enunciado)."""
    e = [-Fo] * n_nos
    f = [1 + 2 * Fo] * n_nos
    g = [-Fo] * n_nos

    # Nó 0 (centro, adiabático): o vizinho "espelhado" dobra o termo da direita
    e[0] = 0.0
    g[0] = -2 * Fo
    # Nó N-1 (superfície, convecção): balanço de energia no meio volume
    e[n_nos - 1] = -2 * Fo
    f[n_nos - 1] = 1 + 2 * Fo + 2 * Fo * Bi
    g[n_nos - 1] = 0.0
    return e, f, g


def montar_b(T, Fo, Bi, q, passo_dt):
    """Vetor de termos independentes {b}, calculado com as temperaturas T atuais."""
    n_nos = len(T)
    # Termo de geração por passo. No PDF aparece A = q/(rho*cp), mas, pela
    # eq. (66) multiplicada por dt, o correto é A = q*dt/(rho*cp) (unidade: °C).
    A = q * passo_dt / (rho * cp)

    b = [0.0] * n_nos
    for i in range(n_nos):
        b[i] = T[i] + A
    # Nó da superfície: T + 2*Fo*Bi*T_inf + Fo*q*dx^2/k, e Fo*q*dx^2/k = A
    b[n_nos - 1] = b[n_nos - 1] + 2 * Fo * Bi * T_inf
    return b


def simular(q, T_ini, n_nos, passo_dt, t_final):
    """Evolui a temperatura de t = 0 até t_final.

    Devolve três coisas (x, tempos, historico):
        x[j]            posição do nó j [m]
        tempos[i]       instante do passo i [s]
        historico[i][j] temperatura no instante i, nó j [°C]
    """
    verificar_entradas(n_nos, passo_dt, t_final, q)
    dx, Fo, Bi = numeros_da_malha(n_nos, passo_dt)

    e, f, g = montar_matriz(n_nos, Fo, Bi)
    e_lu, f_lu = thomas_fatorar(e, f, g)       # fatora uma única vez

    # round() evita o erro de arredondamento que int(t_final/passo_dt) poderia dar
    n_passos = max(1, int(round(t_final / passo_dt)))

    T = [T_ini] * n_nos                         # condição inicial
    historico = [T[:]]                          # guarda t = 0

    for passo in range(n_passos):
        b = montar_b(T, Fo, Bi, q, passo_dt)
        T = thomas_resolver(e_lu, f_lu, g, b)
        for valor in T:                         # detecta divergência numérica
            if math.isnan(valor) or math.isinf(valor):
                raise FloatingPointError(f"Solução divergiu no passo {passo + 1}.")
        historico.append(T[:])

    x = [j * dx for j in range(n_nos)]
    tempos = [i * passo_dt for i in range(n_passos + 1)]
    return x, tempos, historico


# ==============================================================================
# 4. SOLUÇÕES ANALÍTICAS (usadas só para validar o código numérico)
# ==============================================================================

def raizes_autovalores(Bi, n_termos):
    """Primeiras raízes de lambda*tan(lambda) = Bi  (equivale a lambda = Bi*cot(lambda)).

    A n-ésima raiz fica no intervalo (n*pi, n*pi + pi/2), n = 0, 1, 2, ...
    Usamos F(lambda) = lambda*sin(lambda) - Bi*cos(lambda), que não tem
    singularidade, e o método da bisseção.
    """
    def F(lam):
        return lam * math.sin(lam) - Bi * math.cos(lam)

    raizes = []
    for n in range(n_termos):
        a = n * math.pi + 1e-12
        b = n * math.pi + math.pi / 2 - 1e-12
        Fa = F(a)
        for _ in range(200):
            meio = 0.5 * (a + b)
            Fm = F(meio)
            if Fa * Fm <= 0.0:
                b = meio
            else:
                a = meio
                Fa = Fm
            if b - a < 1e-14:
                break
        raizes.append(0.5 * (a + b))
    return raizes


def temperatura_exata(x, t, T_ini, lam):
    """Solução exata por séries (eq. 59) para q = 0, no instante t > 0.

    theta = soma_n [4 sin(l_n) / (2 l_n + sin(2 l_n))] cos(l_n x*) exp(-l_n^2 Fo)
    com theta = (T - T_inf)/(T_ini - T_inf), x* = x/L, Fo = alpha*t/L^2.
    Exige T_ini diferente de T_inf. Não vale em t = 0.
    """
    Fo_t = alpha * t / L ** 2
    coef = [4 * math.sin(l) / (2 * l + math.sin(2 * l)) for l in lam]

    T = [0.0] * len(x)
    for j in range(len(x)):
        soma = 0.0
        for n in range(len(lam)):
            soma = soma + coef[n] * math.cos(lam[n] * x[j] / L) * math.exp(-lam[n] ** 2 * Fo_t)
        T[j] = T_inf + soma * (T_ini - T_inf)
    return T


def temperatura_estacionaria(x, q):
    """Regime permanente com geração: T(x) = T_inf + qL/h + q(L^2 - x^2)/(2k)."""
    T = [0.0] * len(x)
    for j in range(len(x)):
        T[j] = T_inf + q * L / h + q * (L ** 2 - x[j] ** 2) / (2 * k_cond)
    return T


# ==============================================================================
# 5. VALIDAÇÃO
# ==============================================================================

def validar_sem_geracao():
    """q = 0: compara o código completo (com o termo de geração) à solução exata."""
    T_ini_val = T_inf + 300.0      # com q = 0 precisamos de T_ini != T_inf
    x, tempos, hist = simular(0.0, T_ini_val, N, dt, tempo_total)
    lam = raizes_autovalores(Bi_L, N_TERMOS)

    erro_max = 0.0
    t_pior = 0.0
    erro_max_5s = 0.0
    for i in range(1, len(tempos)):               # começa em 1: exata não vale em t = 0
        exata = temperatura_exata(x, tempos[i], T_ini_val, lam)
        for j in range(len(x)):
            erro = abs(hist[i][j] - exata[j])
            if erro > erro_max:
                erro_max = erro
                t_pior = tempos[i]
            if tempos[i] >= 5.0 and erro > erro_max_5s:
                erro_max_5s = erro

    _, Fo, _ = numeros_da_malha(N, dt)
    print("\n[Validação 1] q = 0: numérico x solução exata (séries)")
    print(f"  Bi_L = {Bi_L:.2f} | Fo = {Fo:.3f} | N = {N} | dt = {dt} s")
    print(f"  Erro máximo (t > 0): {erro_max:.3f} °C "
          f"({100 * erro_max / (T_ini_val - T_inf):.2f} % de Ti - T_inf), em t = {t_pior:.2f} s")
    print(f"  Erro máximo para t >= 5 s: {erro_max_5s:.4f} °C")
    return x, tempos, hist, T_ini_val


def validar_estacionario():
    """q > 0: em t grande, a solução numérica deve coincidir com a analítica."""
    t_longo = max(tempo_total, 20 * L ** 2 / alpha)   # ~20 constantes de tempo
    x, tempos, hist = simular(q_ponto, T_inicial, N, dt, t_longo)
    analitica = temperatura_estacionaria(x, q_ponto)

    erro_max = 0.0
    for j in range(len(x)):
        erro = abs(hist[-1][j] - analitica[j])
        if erro > erro_max:
            erro_max = erro

    print("\n[Validação 2] q > 0: regime estacionário numérico x analítico")
    print(f"  t = {tempos[-1]:.0f} s | T_centro numérico = {hist[-1][0]:.4f} °C | "
          f"analítico = {analitica[0]:.4f} °C")
    print(f"  Erro máximo: {erro_max:.2e} °C")


def estudo_refinamento():
    """Erro (q = 0, em t = 10 s) ao refinar a malha e o passo de tempo."""
    T_ini_val = T_inf + 300.0
    lam = raizes_autovalores(Bi_L, N_TERMOS)
    print("\n[Validação 3] Refinamento (q = 0, erro máximo em t = 10 s)")
    print("    N     dt [s]    erro máx [°C]")
    for fator in (1, 2, 4):
        n_nos = (N - 1) * fator + 1        # fator 2: dx cai à metade
        passo = dt / fator ** 2            # fator 2: dt cai a um quarto
        x, tempos, hist = simular(0.0, T_ini_val, n_nos, passo, 10.0)
        exata = temperatura_exata(x, tempos[-1], T_ini_val, lam)
        erro_max = 0.0
        for j in range(len(x)):
            erro_max = max(erro_max, abs(hist[-1][j] - exata[j]))
        print(f"  {n_nos:4d}   {passo:8.5f}   {erro_max:10.5f}")


# ==============================================================================
# 6. GRÁFICOS
# ==============================================================================

def finalizar_figura(nome_arquivo):
    """Salva a figura em PNG (e mostra na tela, se MOSTRAR_FIGURAS = True)."""
    plt.tight_layout()
    os.makedirs(PASTA_RESULTADOS, exist_ok=True)
    plt.savefig(os.path.join(PASTA_RESULTADOS, nome_arquivo), dpi=150)
    if MOSTRAR_FIGURAS:
        plt.show()
    plt.close()


def indice_do_tempo(t_alvo, passo_dt, n_passos):
    """Índice do passo mais próximo de t_alvo; None se passar do tempo simulado."""
    i = int(round(t_alvo / passo_dt))
    if i > n_passos:
        return None
    return i


def em_mm(x):
    """Converte uma lista de posições de metros para milímetros."""
    return [xi * 1e3 for xi in x]


def grafico_validacao(x, tempos, hist, T_ini_val):
    """Perfis numéricos (pontos) sobre os exatos (linhas), caso q = 0."""
    lam = raizes_autovalores(Bi_L, N_TERMOS)
    plt.figure(figsize=(9, 5.5))
    for t_alvo in (0.5, 2, 5, 10, 20):
        i = indice_do_tempo(t_alvo, dt, len(tempos) - 1)
        if i is None or i == 0:
            continue
        exata = temperatura_exata(x, tempos[i], T_ini_val, lam)
        linha, = plt.plot(em_mm(x), exata, "-", linewidth=1.5)
        plt.plot(em_mm(x), hist[i], "o", markersize=4, color=linha.get_color(),
                 label=f"t = {tempos[i]:.1f} s")
    plt.title("Validação (q = 0): numérico (pontos) x exato (linhas)")
    plt.xlabel("Posição a partir do centro (mm)")
    plt.ylabel("Temperatura (°C)")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.6)
    finalizar_figura("validacao_sem_geracao.png")


def grafico_perfis(x, tempos, hist):
    """Perfis em vários instantes e evolução da temperatura no centro/superfície."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    estacionaria = temperatura_estacionaria(x, q_ponto)
    n_passos = len(tempos) - 1

    for t_alvo in (0, 1, 3, 5, 10, 30, tempo_total):
        i = indice_do_tempo(t_alvo, dt, n_passos)
        if i is not None:
            ax1.plot(em_mm(x), hist[i], marker=".", label=f"t = {tempos[i]:.1f} s")
    ax1.plot(em_mm(x), estacionaria, "k--", linewidth=1.2, label="estacionário (analítico)")
    ax1.set_title("Perfis de temperatura no combustível (UO$_2$)")
    ax1.set_xlabel("Posição a partir do centro (mm)")
    ax1.set_ylabel("Temperatura (°C)")
    ax1.legend()
    ax1.grid(True, linestyle="--", alpha=0.6)

    T_centro = [linha[0] for linha in hist]
    T_superficie = [linha[-1] for linha in hist]
    ax2.plot(tempos, T_centro, label="centro (x = 0)")
    ax2.plot(tempos, T_superficie, label="superfície (x = L)")
    ax2.axhline(estacionaria[0], color="k", linestyle="--", linewidth=1,
                label="estacionário (centro)")
    ax2.set_title("Evolução temporal")
    ax2.set_xlabel("Tempo (s)")
    ax2.set_ylabel("Temperatura (°C)")
    ax2.legend()
    ax2.grid(True, linestyle="--", alpha=0.6)
    finalizar_figura("perfis_com_geracao.png")


def grafico_campo(x, tempos, hist):
    """Mapa de cores: posição x tempo."""
    plt.figure(figsize=(9, 5.5))
    X, Tempo = np.meshgrid(em_mm(x), tempos)     # grades com as mesmas dimensões de hist
    mapa = plt.contourf(X, Tempo, np.array(hist), levels=50, cmap="inferno")
    plt.colorbar(mapa).set_label("Temperatura (°C)")
    plt.title("Campo térmico: posição x tempo")
    plt.xlabel("Posição a partir do centro (mm)")
    plt.ylabel("Tempo (s)")
    finalizar_figura("campo_termico.png")


# ==============================================================================
# 7. PROGRAMA PRINCIPAL
# ==============================================================================

def main():
    dx, Fo, Bi_malha = numeros_da_malha(N, dt)
    print("=== PPC3: condução de calor transiente 1D em combustível UO2 ===")
    print(f"alpha = {alpha:.3e} m2/s | dx = {dx * 1e3:.3f} mm | dt = {dt} s | "
          f"Fo = {Fo:.3f} | Bi_malha = {Bi_malha:.3f}")
    print("(Esquema implícito: estável para qualquer Fo; o explícito exigiria Fo <= 0,5.)")

    # --- Validações ---
    x_v, t_v, hist_v, T_ini_val = validar_sem_geracao()
    validar_estacionario()
    estudo_refinamento()

    # --- Caso com geração interna ---
    x, tempos, hist = simular(q_ponto, T_inicial, N, dt, tempo_total)
    print("\n[Resultado] Com geração interna")
    print(f"  t = {tempos[-1]:.0f} s | T_centro = {hist[-1][0]:.2f} °C | "
          f"T_superfície = {hist[-1][-1]:.2f} °C")

    # --- Visualização ---
    grafico_validacao(x_v, t_v, hist_v, T_ini_val)
    grafico_perfis(x, tempos, hist)
    grafico_campo(x, tempos, hist)
    print(f"\nFiguras salvas na pasta '{PASTA_RESULTADOS}/'.")


if __name__ == "__main__":
    main()