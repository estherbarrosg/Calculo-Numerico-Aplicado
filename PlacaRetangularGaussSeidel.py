import numpy as np


# ---------------------------------------------------------------------------
# 1. Malha e condições de contorno
# ---------------------------------------------------------------------------
def criar_placa(h, Lx=2.0, Ly=1.0, T_superior=50.0, T_inferior=0.0, T_esquerda=100.0, T_direita=0.0, T_inicial=0.0):
    nx = round(Lx / h) -1 #nós internos em x
    ny = round(Ly / h) -1 #nós internos em y
    
    
    T = np.full((ny+2, nx+2), T_inicial)
    T[0, 1:nx + 1]       = T_inferior   # linha de baixo
    T[ny + 1, 1:nx + 1]  = T_superior   # linha de cima
    T[1:ny + 1, 0]       = T_esquerda   # coluna da esquerda
    T[1:ny + 1, nx + 1]  = T_direita    # coluna da direita
    return T
    


# ---------------------------------------------------------------------------
# 2. Gauss-Seidel sem montar matriz (equação de Laplace)
#    T[i, j] = ( T[i-1, j] + T[i+1, j] + T[i, j-1] + T[i, j+1] ) / 4
# ---------------------------------------------------------------------------
def gauss_seidel_malha_2d(T, tol=1e-6, max_iter=100000):
    """
    Relaxa os nós INTERNOS de T (já com as bordas preenchidas).
    T é modificada no lugar.

    Retorna: iteracoes, erro_final, convergiu
    """
    N_lin, N_col = T.shape

    iteracoes = 0
    erro = np.inf

    while erro > tol and iteracoes < max_iter:
        erro = 0.0

        for j in range(1, N_col - 1):          # só nós internos
            for i in range(1, N_lin - 1):
                antigo = T[i, j]

                T[i, j] = (T[i - 1, j] + T[i + 1, j]
                           + T[i, j - 1] + T[i, j + 1]) / 4.0

                erro = max(erro, abs(T[i, j] - antigo))

        iteracoes += 1

    convergiu = erro <= tol
    return iteracoes, erro, convergiu
 
 
# ---------------------------------------------------------------------------
# 3. Conferência adicional: max |r_ij|, r_ij = 4*T_ij - (soma dos 4 vizinhos)
# ---------------------------------------------------------------------------
def residuo_maximo_2d(T):
    N_lin, N_col = T.shape
    max_res = 0.0

    for j in range(1, N_col - 1):
        for i in range(1, N_lin - 1):
            res = abs(4 * T[i, j] - (T[i - 1, j] + T[i + 1, j]
                                     + T[i, j - 1] + T[i, j + 1]))
            if res > max_res:
                max_res = res

    return max_res



# ---------------------------------------------------------------------------
# 4. Apresentar os resultados
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    valores_h = [0.25, 0.125, 0.0625]
    resultados = {}

    for h in valores_h:
        T = criar_placa(h)                                        # monta a placa
        iters, erro_final, convergiu = gauss_seidel_malha_2d(T)   # resolve
        max_res = residuo_maximo_2d(T)                            # confere
        resultados[h] = T

        print(f"\n--- Resultados para h = {h} ---")
        print(f"Convergiu: {convergiu}")
        print(f"Iterações necessárias: {iters}")
        print(f"Erro final: {erro_final:.2e}")
        print(f"Máximo resíduo |r_ij|: {max_res:.2e}")

        if h == 0.25:
            # T(4,4) do PDF (base 1) = T[3, 3] em Python (base 0)
            i_p = round(0.5 / h)   # linha  -> y = 0,5 m
            j_p = round(1.0 / h)   # coluna -> x = 1,0 m
            print(f"Temperatura em (x=1,0 m; y=0,5 m): {T[i_p, j_p]:.4f} °C")
        # item (a): salvar a matriz (incluindo as bordas) em arquivo de dados
        np.savetxt(f"temperaturas_h{h}.txt", T, fmt="%.4f")
