import math

# Matriz A (3x3)
A = [
    [2000.0, -800.0,  0.0, 0.0, 0.0],  # Linha 0
    [-800.0, 2300.0, -1500.0, 0.0, 0.0],  # Linha 1
    [0.0, -1500,  2500.0, -1000.0, 0.0],   # Linha 2
    [0.0, 0.0, -1000.0 , 2800.0, -1800.0], # Linha 3
    [0.0, 0.0, 0.0, -1800.0, 2700.0], # Linha 4
]
   
n = 5

# Passo 1: Inicializar L como Identidade e U como cópia de A
L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
U = [linha[:] for linha in A]

# Passo 2: Laços de eliminação (Decomposição LU)
for k in range(n - 1):          # k é a coluna do pivô atual (vai de 0 até n-2)
    for i in range(k + 1, n):   # i percorre as linhas ABAIXO do pivô (de k+1 até n-1)
        
        # Calcula o multiplicador que vai zerar U[i][k]
        fator = U[i][k] / U[k][k]
        
        # Guarda esse multiplicador na mesma posição da matriz L
        L[i][k] = fator
        
        # Atualiza o restante da linha i na matriz U
        for j in range(k + 1, n):
            U[i][j] = U[i][j] - fator * U[k][j]
            
        # Zera explicitamente o elemento abaixo do pivô em U
        U[i][k] = 0.0

# Passo 3: Imprimir os resultados
print("Matriz L:")
for linha in L:
    print([f"{val:8.4f}" for val in linha])

print("\nMatriz U:")
for linha in U:
    print([f"{val:8.4f}" for val in linha])