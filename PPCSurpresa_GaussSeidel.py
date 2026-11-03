import numpy as np
import matplotlib.pyplot as plt

def resolver_placa_gauss_seidel(n, tol=1e-6, max_iter=100000):
    """
    Resolve a distribuição de temperatura na placa 2D usando Gauss-Seidel.
    """
    # 1. Definir a malha e as condições de contorno
    # Criar matriz com (n+2) x (n+2) inicializada a 20°C
    T = np.full((n+2, n+2), 20.0)
    
    # Condição de contorno: borda inferior a 100°C.
    # Em Python (0-indexado), se j é a linha (eixo y) e i é a coluna (eixo x),
    # definimos a linha 0 (base) como 100°C, entre i=1 e i=n.
    T[0, 1:n+1] = 100.0
    
    iteracoes = 0
    erro = float('inf')
    
    # 2 e 3. Implementar a atualização iterativa e o critério de paragem
    while erro > tol and iteracoes < max_iter:
        erro = 0.0
        # Percorrer os nós internos
        for j in range(1, n+1):
            for i in range(1, n+1):
                antigo = T[j, i]
                # Atualização imediata (Gauss-Seidel)
                T[j, i] = (T[j, i-1] + T[j, i+1] + T[j-1, i] + T[j+1, i]) / 4.0
                
                # Calcular a diferença máxima
                diferenca = abs(T[j, i] - antigo)
                if diferenca > erro:
                    erro = diferenca
                    
        iteracoes += 1
        
    # Conferência adicional: calcular max |r_i,j|
    max_res = 0.0
    for j in range(1, n+1):
        for i in range(1, n+1):
            res = abs(4 * T[j, i] - (T[j, i-1] + T[j, i+1] + T[j-1, i] + T[j+1, i]))
            if res > max_res:
                max_res = res
                
    return T, iteracoes, erro, max_res

# ---------------------------------------------------------
# 4. Apresentar os resultados
# ---------------------------------------------------------
valores_n = [5, 10, 20]
resultados = {}

for n in valores_n:
    T_final, iters, erro_final, max_res = resolver_placa_gauss_seidel(n)
    resultados[n] = T_final
    
    print(f"\n--- Resultados para n = {n} ---")
    print(f"Iterações necessárias: {iters}")
    print(f"Erro final: {erro_final:.2e}")
    print(f"Máximo resíduo |r_ij|: {max_res:.2e}")
    
    # Nó central para n=5
    # Num array de 7x7 (índices 0 a 6), o nó interno central está no índice 3,3.
    # Isto corresponde ao T(4,4) da notação 1-indexada da folha de exercícios.
    if n == 5:
        print(f"Temperatura do nó central T(3,3) em Python [equivale a T(4,4)]: {T_final[3, 3]:.4f} °C")

# Representação gráfica (Gráficos de calor)
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

for idx, n in enumerate(valores_n):
    ax = axes[idx]
    T = resultados[n]
    
    # origin='lower' garante que a borda de 100ºC (j=0) fique na base do gráfico
    im = ax.imshow(T, origin='lower', cmap='jet', extent=[0, 1, 0, 1], vmin=20, vmax=100)
    ax.set_title(f'Distribuição de Temperatura (n = {n})')
    ax.set_xlabel('x')
    if idx == 0:
        ax.set_ylabel('y')

# Barra de cor global
cbar = fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.02, pad=0.04)
cbar.set_label('Temperatura (°C)')

plt.suptitle('Método de Gauss-Seidel: Condução de Calor numa Placa 2D', fontsize=16)
plt.show()


# ==============================================================================
# ALGORITMO GENÉRICO DO MÉTODO DE GAUSS-SEIDEL (COM E SEM RELAXAMENTO)
# ==============================================================================

def gauss_seidel(A, b, tol=1e-5, max_iter=100, lambda_val=1.0):
    """
    Resolve o sistema linear Ax = b utilizando o método iterativo de Gauss-Seidel.
    
    Parâmetros:
    A          : Matriz dos coeficientes (lista de listas).
    b          : Vetor de termos independentes (lista).
    tol        : Tolerância para o critério de parada (erro relativo).
    max_iter   : Número máximo de iterações permitidas.
    lambda_val : Fator de relaxamento. 
                 Se = 1.0, é o Gauss-Seidel padrão[cite: 18].
                 Se > 1.0, é sobre-relaxamento (SOR)[cite: 18].
                 Se < 1.0, é sub-relaxamento[cite: 18].
    """
    n = len(b)
    
    # 1. Chute inicial: Inicializamos o vetor solução (x) com zeros.
    x = [0.0] * n 
    
    print("-" * 60)
    print(f"Iniciando Gauss-Seidel (lambda = {lambda_val})")
    print("-" * 60)
    
    # 2. Laço principal de iterações
    for k in range(max_iter):
        
        # Guardamos uma cópia exata do vetor 'x' da iteração anterior.
        # Isso é fundamental para calcularmos o erro e aplicarmos a relaxação[cite: 18].
        x_ant = x[:]
        erro_max = 0.0
        
        # 3. Laço para calcular cada incógnita x_i da iteração atual
        for i in range(n):
            
            # Prevenção: O método falha se o elemento da diagonal principal for zero
            if A[i][i] == 0.0:
                print(f"Erro: O elemento da diagonal principal A[{i}][{i}] é nulo.")
                return None
            
            soma = b[i]
            
            # 4. Somatório dos termos (A_ij * x_j) que passam para o lado direito da equação[cite: 18]
            for j in range(n):
                if i != j:
                    # A grande vantagem do Gauss-Seidel: o vetor 'x' está a ser atualizado 
                    # em tempo real. Logo, x[j] utiliza SEMPRE a estimativa mais recente disponível[cite: 18].
                    soma -= A[i][j] * x[j]
            
            # 5. Cálculo do novo valor da variável isolada (antes da relaxação)
            x_calc = soma / A[i][i]
            
            # 6. Aplicação do parâmetro de relaxação[cite: 18]
            # Mistura o valor recém-calculado com o valor da iteração anterior.
            x[i] = lambda_val * x_calc + (1.0 - lambda_val) * x_ant[i]
            
            # 7. Cálculo do erro relativo para a variável atual
            # Evitamos divisão por zero caso a estimativa atual seja exatemente nula
            if x[i] != 0.0:
                erro_rel = abs((x[i] - x_ant[i]) / x[i])
                
                # O erro da iteração é o MAIOR erro relativo entre todas as variáveis
                if erro_rel > erro_max:
                    erro_max = erro_rel
                    
        # Exibição do progresso no terminal
        vetor_format = [f"{val:8.4f}" for val in x]
        print(f"Iteração {k+1:2d} | x = {vetor_format} | Erro Relativo Máx = {erro_max*100:5.2f}%")
        
        # 8. Critério de Parada
        # Se o maior erro relativo for menor que a tolerância exigida, o sistema convergiu!
        if erro_max < tol:
            print("-" * 60)
            print(f">>> Convergência atingida com sucesso em {k+1} iterações!")
            return x
            
    # Se o loop terminar sem atingir o 'return' da convergência:
    print("-" * 60)
    print("Aviso: Número máximo de iterações atingido sem convergência.")
    print("Verifique se a matriz possui dominância diagonal.")
    return x


# ==============================================================================
# BLOCO PRINCIPAL (EXEMPLO DE COMO USAR A FUNÇÃO)
# ==============================================================================
if __name__ == "__main__":
    
    # Exemplo genérico de um sistema 3x3 com dominância diagonal
    A_exemplo = [
        [ 8.0, -1.0,  2.0],
        [-2.0, 10.0, -3.0],
        [ 1.0,  2.0, 12.0]
    ]
    b_exemplo = [20.0, 33.0, 40.0]
    
    # Chama a função estipulando uma tolerância de 5% (0.05) e relaxamento padrão (1.0)[cite: 19]
    solucao = gauss_seidel(A_exemplo, b_exemplo, tol=0.05, lambda_val=1.0)
    
    if solucao:
        print("\nVetor Solução Final:")
        for i in range(len(solucao)):
            print(f"x_{i+1} = {solucao[i]:.6f}")










# ---------------------------------------------------------
# 5. Questões rápidas (Respostas em forma de comentários)
# ---------------------------------------------------------
"""
a) Por que as temperaturas das bordas não são atualizadas?
R: As temperaturas das bordas representam as condições de contorno de Dirichlet 
   do problema. São valores fixos conhecidos desde o início (neste caso, a fonte 
   térmica de 100ºC e as margens a 20ºC). O objetivo é descobrir a distribuição 
   apenas nos nós internos.

b) O que diferencia a atualização de Gauss-Seidel da de Jacobi?
R: No método de Jacobi, a atualização de todos os nós baseia-se estritamente 
   nos valores da iteração anterior. No método de Gauss-Seidel (utilizado acima), 
   os valores que acabaram de ser atualizados num passo iterativo são utilizados 
   imediatamente no cálculo dos nós subsequentes da mesma iteração, o que 
   habitualmente acelera a convergência do sistema.

c) Como o número de iterações muda quando a malha é refinada?
R: O número de iterações aumenta drasticamente à medida que a malha é refinada 
   (n maior). Como a distância entre os nós (delta x e delta y) diminui, são 
   necessários muito mais passos iterativos para que a energia e a informação 
   das condições de contorno se propaguem por todo o domínio até o sistema 
   atingir o estado estacionário dentro da tolerância exigida.
"""




