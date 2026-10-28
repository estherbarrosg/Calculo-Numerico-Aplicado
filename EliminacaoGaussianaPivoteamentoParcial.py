import math

#1. Definicao do Sistema Linear (Exemplo: Ax = b)
A = [
    [ 0.0,  2.0,  1.0],  # O pivo inicial nulo faria o metodo ingenuo falhar imediatamente
    [ 1.0, -2.0,  3.0],
    [-1.0,  1.0,  2.0]
]

b = [3.0, 2.0, 2.0]

n = len(b)  # Dimensao do sistema (n = 3)

# =================================================================
# FASE 1: ELIMINACAO PROGRESSIVA COM PIVOTEAMENTO PARCIAL
# =================================================================

for k in range(n - 1):
    
    # --- ETAPA DE PIVOTEAMENTO ---
    # Encontrar o indice da linha com o maior valor absoluto na coluna k
    max_index = k
    maior_valor = abs(A[k][k])
    
    for i in range(k + 1, n):
        if abs(A[i][k]) > maior_valor:
            maior_valor = abs(A[i][k])
            max_index = i
            
    # Se o maior valor absoluto for zero, a matriz e singular (nao ha solucao unica)
    if maior_valor == 0.0:
        print(f"Erro: O sistema e singular (coluna {k} nula da diagonal para baixo).")
        break
        
    # Trocar as linhas na matriz A e no vetor b, se o maior elemento nao estiver na diagonal
    if max_index != k:
        A[k], A[max_index] = A[max_index], A[k]
        b[k], b[max_index] = b[max_index], b[k]
    # -----------------------------
    
    # --- ETAPA DE ELIMINACAO ---
    # Com a garantia de que A[k][k] e o melhor pivo possivel, fazemos o escalonamento
    for i in range(k + 1, n):
        m = A[i][k] / A[k][k]
        
        # Atualiza a linha i da matriz A
        for j in range(k, n):
            A[i][j] = A[i][j] - m * A[k][j]
            
        # Atualiza o vetor independente b
        b[i] = b[i] - m * b[k]

print("Matriz A escalonada apos o Pivoteamento:")
for linha in A:
    print([f"{val:8.4f}" for val in linha])

# =================================================================
# FASE 2: SUBSTITUICAO REGRESSIVA
# =================================================================

# Cria o vetor solucao x cheio de zeros
x = [0.0 for _ in range(n)]

# Verifica se a ultima equacao e viavel
if A[n - 1][n - 1] == 0.0:
    print("\nErro: O sistema e singular (ultimo pivo nulo). Nao e possivel retroceder.")
else:
    # A ultima equacao so tem uma incognita, isolamos diretamente
    x[n - 1] = b[n - 1] / A[n - 1][n - 1]

    # Percorre as linhas de baixo para cima (da penultima ate a primeira)
    for i in range(n - 2, -1, -1):
        soma = b[i]
        
        # Subtrai os termos cujas incognitas ja foram descobertas
        for j in range(i + 1, n):
            soma = soma - A[i][j] * x[j]
            
        # Divide pelo pivo da linha atual
        x[i] = soma / A[i][i]

    # Imprime o vetor solucao final
    print("\nVetor Solucao (x) com Pivoteamento Parcial:")
    for i in range(n):
        print(f"x[{i}] = {x[i]:.6f}")
        
        
        
# 1. Parametros geometricos (Alterar theta_1=0.01 e theta_2=0.02 para o Cenario C)
theta_1 = 15.0
theta_2 = 65.0
phi = 45.0

sin_t1 = math.sin(math.radians(theta_1))
sin_t2 = math.sin(math.radians(theta_2))
cos_t1 = math.cos(math.radians(theta_1))
cos_t2 = math.cos(math.radians(theta_2))
sin_p = math.sin(math.radians(phi))
cos_p = math.cos(math.radians(phi))

# 2. Montagem da Matriz A (5x5) e do vetor de referencia F_ref
A = [
    [sin_t1, sin_t2,  1.0,    0.0,  0.0],
    [cos_t1, cos_t2,  0.0,    0.0,  0.0],
    [   0.0,    0.0, -1.0,  cos_p,  0.0],
    [   0.0,    0.0,  0.0, -sin_p, -1.0],
    [   0.0,    0.0,  0.0,    0.0,  1.0]
]

F_ref = [100.0, -80.0, 20.0, 35.0, 15.0]
n = len(F_ref)

# Calculo do vetor de cargas consistente: b = A * F_ref
b = [0.0 for _ in range(n)]
for i in range(n):
    soma = 0.0
    for j in range(n):
        soma = soma + A[i][j] * F_ref[j]
    b[i] = soma

# Salvaguarda das copias originais intactas para calculo de residuos e perturbacao
A_orig = [linha[:] for linha in A]
b_orig = b[:]

print("Vetor de cargas b gerado:")
print([f"{valor:.6f}" for valor in b])

# ==============================================================================
# FASE 1: ELIMINACAO PROGRESSIVA COM PIVOTEAMENTO PARCIAL
# ==============================================================================
limite_pivo = 1e-12

for k in range(n - 1):
    # Etapa de Pivoteamento: busca pelo maior valor absoluto na coluna k
    max_index = k
    maior_valor = abs(A[k][k])
    
    for i in range(k + 1, n):
        if abs(A[i][k]) > maior_valor:
            maior_valor = abs(A[i][k])
            max_index = i
            
    # Verificacao de singularidade numerica
    if maior_valor <= limite_pivo:
        print(f"Erro: Sistema singular ou quase singular na coluna {k}.")
        break
        
    # Troca de linhas em A e em b se necessario
    if max_index != k:
        A[k], A[max_index] = A[max_index], A[k]
        b[k], b[max_index] = b[max_index], b[k]
        
    # Etapa de Eliminacao abaixo do pivo
    for i in range(k + 1, n):
        m = A[i][k] / A[k][k]
        for j in range(k, n):
            A[i][j] = A[i][j] - m * A[k][j]
        b[i] = b[i] - m * b[k]

# ==============================================================================
# FASE 2: SUBSTITUICAO REGRESSIVA
# ==============================================================================
x = [0.0 for _ in range(n)]

if abs(A[n - 1][n - 1]) <= limite_pivo:
    print("Erro: Ultimo pivo nulo. Nao e possivel realizar a substituicao.")
else:
    x[n - 1] = b[n - 1] / A[n - 1][n - 1]
    for i in range(n - 2, -1, -1):
        soma = b[i]
        for j in range(i + 1, n):
            soma = soma - A[i][j] * x[j]
        x[i] = soma / A[i][i]

F_calc = x
print("\nVetor de Forcas Calculado (F_calc):")
for i in range(n):
    print(f"F[{i+1}] = {F_calc[i]:.6f} N")

# ==============================================================================
# FASE 3: VERIFICACAO DE ERRO E RESIDUO (r = A_orig * F_calc - b_orig)
# ==============================================================================
erro = [abs(F_calc[i] - F_ref[i]) for i in range(n)]
residuo = [0.0 for _ in range(n)]

for i in range(n):
    soma = 0.0
    for j in range(n):
        soma = soma + A_orig[i][j] * F_calc[j]
    residuo[i] = soma - b_orig[i]

norma_erro = max(erro)
norma_residuo = max(abs(val) for val in residuo)

print(f"\nErro Maximo face a F_ref: {norma_erro:.6e} N")
print(f"Norma Infinito do Residuo: {norma_residuo:.6e}")

# ==============================================================================
# FASE 4: ANALISE DE SENSIBILIDADE (PERTURBACAO EM b)
# ==============================================================================
copia_b = b_orig[:]
copia_b[1] = copia_b[1] + 0.01  # Perturbacao de +0.01 N na segunda componente

# Copias de trabalho para nao destruir A_orig nem copia_b
A_pert = [linha[:] for linha in A_orig]
b_para_resolver = copia_b[:]

for k in range(n - 1):
    max_index = k
    maior_valor = abs(A_pert[k][k])
    for i in range(k + 1, n):
        if abs(A_pert[i][k]) > maior_valor:
            maior_valor = abs(A_pert[i][k])
            max_index = i
    if max_index != k:
        A_pert[k], A_pert[max_index] = A_pert[max_index], A_pert[k]
        b_para_resolver[k], b_para_resolver[max_index] = b_para_resolver[max_index], b_para_resolver[k]
    for i in range(k + 1, n):
        m = A_pert[i][k] / A_pert[k][k]
        for j in range(k, n):
            A_pert[i][j] = A_pert[i][j] - m * A_pert[k][j]
        b_para_resolver[i] = b_para_resolver[i] - m * b_para_resolver[k]

F_pert = [0.0 for _ in range(n)]
F_pert[n - 1] = b_para_resolver[n - 1] / A_pert[n - 1][n - 1]
for i in range(n - 2, -1, -1):
    soma = b_para_resolver[i]
    for j in range(i + 1, n):
        soma = soma - A_pert[i][j] * F_pert[j]
    F_pert[i] = soma / A_pert[i][i]

delta_F = [abs(F_pert[i] - F_ref[i]) for i in range(n)]
residuo_pert = [0.0 for _ in range(n)]
for i in range(n):
    soma = 0.0
    for j in range(n):
        soma = soma + A_orig[i][j] * F_pert[j]
    residuo_pert[i] = soma - copia_b[i]

print(f"\n--- Resultado apos Perturbacao de +0.01 N ---")
print(f"Norma Infinito do Residuo 2 (perturbado): {max(abs(v) for v in residuo_pert):.6e}")
print(f"Norma Infinito de Delta F: {max(delta_F):.6f} N")



# ==============================================================================
# 10. Questões para discussão (Respostas)
# ==============================================================================

# a) Qual a relação entre a diferença angular das barras e o determinante do sistema 2x2?
# O módulo do determinante é dado por |sin(θ1 - θ2)|[cite: 40]. Quando as barras estão 
# com ângulos muito próximos, a diferença angular (θ1 - θ2) aproxima-se de zero, o que 
# faz com que o seno dessa diferença, e consequentemente o determinante da matriz, também 
# se aproxime de zero[cite: 40]. Um determinante muito próximo de zero indica que a matriz 
# está mal condicionada, tornando o sistema extremamente sensível a pequenas perturbações 
# nas cargas[cite: 39, 40]. Se as barras forem exatamente paralelas (θ1 = θ2), o 
# determinante será nulo, e o sistema não terá solução única[cite: 40].

# b) No sistema 5x5, quais linhas foram escolhidas como pivô em cada cenário? Por quê?
# A escolha da linha pivô no método com pivoteamento parcial é feita selecionando, na 
# coluna atual, o elemento de maior valor absoluto[cite: 42]. No cenário A (barras 
# afastadas com θ1=15°), o primeiro elemento A[0][0] = sin(15°) já possui um valor 
# razoável, embora o A[1][0] = cos(15°) seja maior. Logo, a troca de linhas ocorrerá 
# para maximizar esse valor e evitar divisão por números pequenos. No cenário C (barras 
# quase paralelas, θ1=0.01°), o primeiro elemento A[0][0] = sin(0.01°) é extremamente 
# próximo de zero[cite: 40, 41]. Sem pivoteamento, o algoritmo usaria este valor ínfimo 
# como divisor, gerando multiplicadores imensos e amplificando erros de arredondamento 
#[cite: 39, 42]. Com o pivoteamento parcial, a linha 1 (índice 0) será trocada pela 
# linha 2 (índice 1), pois cos(0.01°) é quase 1, selecionando o maior pivô possível 
# na coluna k=1 e garantindo a estabilidade numérica da eliminação[cite: 41, 42].

# c) O pivoteamento parcial alterou o sistema físico ou apenas a ordem de processamento das equações?
# O pivoteamento parcial altera única e exclusivamente a ordem de processamento matemático 
# das equações (troca de linhas da matriz A simultaneamente com a troca no vetor b)[cite: 42]. 
# Como as equações representam balanços de forças estáticas que devem ser satisfeitos 
# simultaneamente[cite: 39, 41], a ordem em que essas equações são escritas não altera 
# em absolutamente nada o sistema físico, a geometria da estrutura ou a resposta real 
# dos esforços nas barras[cite: 39].

# d) Um resíduo pequeno garante forças calculadas próximas das forças exatas em um problema mal condicionado? Justifique.
# Não. Em sistemas mal condicionados (como o cenário C com barras quase paralelas), a 
# matriz é extremamente sensível[cite: 39]. Isso significa que mesmo gerando um vetor 
# residual numérico muito pequeno (r = A*F_calc - b ≈ 0)[cite: 43], as forças calculadas 
# (F_calc) podem desviar drasticamente das forças de referência exatas (F_ref) se houver o 
# mínimo erro de arredondamento[cite: 39]. O pequeno resíduo apenas indica que a 
# matemática da máquina fechou a conta do sistema perturbado, mas não garante a precisão 
# física da solução[cite: 44].

# e) Qual é a diferença entre estabilidade do algoritmo e condicionamento do problema?
# O condicionamento é uma propriedade intrínseca do problema físico e matemático 
# (depende exclusivamente da geometria das barras na matriz A) e dita como pequenas 
# variações nas cargas afetam os esforços axiais[cite: 39, 40]. Um problema mal 
# condicionado (barras paralelas) continuará mal condicionado independentemente do 
# método usado para o resolver[cite: 39, 43]. Já a estabilidade do algoritmo 
# refere-se à capacidade do método numérico (como a eliminação gaussiana) não 
# amplificar erros durante os cálculos[cite: 39, 42]. Ao usar pivoteamento parcial, 
# conferimos estabilidade ao algoritmo, mas não resolvemos a sensibilidade física 
# de um problema mal condicionado[cite: 39, 43].

# f) O que ocorreria se θ1 = θ2 exatamente?
# Se os ângulos forem idênticos (θ1 = θ2), as duas barras estariam perfeitamente sobrepostas, 
# tornando-se geometricamente colineares[cite: 39, 40, 44]. Matematicamente, a primeira e a 
# segunda coluna da matriz A passariam a ser idênticas[cite: 41], tornando a matriz 
# singular (determinante exatamente zero). O sistema seria matematicamente insolúvel (sem 
# solução única), e o código computacional falharia ao encontrar um pivô máximo igual a zero 
# em uma das colunas, forçando a interrupção do programa[cite: 42].