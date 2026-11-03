import math

# Definindo a equação do problema
def f(H):
    numerador = 0.471405 * (20.0 * H) ** (5 / 3)
    denominador = (20.0 + 2*H) ** (2 / 3)
    return (numerador / denominador) - 5

# Configuração do método
Hl = 0.5
Hu = 1.0
tolerancia = 1e-6
max_iter = 50

# Verificando os limites iniciais (função com f minúsculo e sem break)
if f(Hl) * f(Hu) > 0:
    print("Erro: não há mudança de sinal nesse intervalo")
else:
    print("Iniciando o Método da Falsa Posição")
    print("-" * 65)
    # Cabeçalho para manter as colunas alinhadas
    print(f"{'i':>2} | {'Hl':>8} | {'Hu':>8} | {'Hr':>8} | {'f(Hr)':>12} | {'Erro (%)':>10}")
    print("-" * 65)
    
    h_antigo = Hl
    
    # O laço de repetição (agora perfeitamente alinhado)
    for i in range(1, max_iter + 1):
        Hr = (Hl * f(Hu) - Hu * f(Hl)) / (f(Hu) - f(Hl))
        
        # Imprime os dados (sem quebrar a linha)
        print(f"{i:2d} | {Hl:.6f} | {Hu:.6f} | {Hr:.6f} | {f(Hr):12.6f} ", end="")
   
        # Calculando o erro apenas a partir da rodada 2
        if i > 1:
            ea = abs((Hr - h_antigo) / h_antigo) * 100
            print(f"| {ea:10.6f}%")
            
            # Critério de parada
            if abs(f(Hr)) < tolerancia:
                print("-" * 65)
                print(f"Raiz Encontrada! Profundidade H = {Hr:.6f} metros.")
                print(f"Total de iterações: {i}")
                break
        else:
            # Esse else pertence ao if i > 1
            print("|        ---")
            
        # Teorema de Bolzano (agora fora do if i > 1, mas dentro do for)
        if f(Hl) * f(Hr) < 0:
            Hu = Hr
        else:
            Hl = Hr
            
        h_antigo = Hr
            
# Cálculo final da velocidade
# Ele fica colado na margem esquerda para rodar apenas depois que o laço terminar
Q = 5.0
B = 20.0
U = Q / (B * Hr)
print(f"\nVelocidade média U = {U:.4f} m/s")


import math

def falsa_posicao(f, a, b, tol=1e-6, max_iter=100):
    """
    Encontra a raiz de f(x) no intervalo [a, b] usando o Método da Falsa Posição.
    """
    # 1. Teorema de Bolzano: Verificar se existe uma raiz no intervalo
    if f(a) * f(b) >= 0:
        print("Erro: O intervalo [a, b] é inválido. f(a) e f(b) devem ter sinais opostos.")
        return None
    
    # Cabeçalho da tabela de iterações
    print(f"{'Iter':>4} | {'a':>10} | {'b':>10} | {'x':>10} | {'f(x)':>14} | {'Erro':>10}")
    print("-" * 75)
    
    x_anterior = a # Apenas para inicializar o cálculo do erro
    
    for i in range(1, max_iter + 1):
        fa = f(a)
        fb = f(b)
        
        # 2. Fórmula da Falsa Posição (cruzamento da reta secante com o eixo x)
        x = (a * fb - b * fa) / (fb - fa)
        fx = f(x)
        
        # Cálculo do erro (distância entre a nova estimativa e a anterior)
        erro = abs(x - x_anterior)
        
        # Imprime a linha da iteração atual
        print(f"{i:4d} | {a:10.6f} | {b:10.6f} | {x:10.6f} | {fx:14.6e} | {erro:10.6e}")
        
        # 3. Critérios de Parada
        # Para a função se |f(x)| for praticamente zero OU se a estimativa x parar de mudar
        if abs(fx) < tol or erro < tol:
            print("-" * 75)
            print(f"Raiz encontrada com sucesso na iteração {i}!")
            return x
            
        # 4. Atualização do intervalo [a, b]
        if fa * fx < 0:
            b = x  # A raiz está no subintervalo esquerdo [a, x]
        else:
            a = x  # A raiz está no subintervalo direito [x, b]
            
        x_anterior = x
        
    print("\nAviso: Número máximo de iterações atingido sem convergir totalmente.")
    return x


# ==============================================================================
# EXEMPLO DE USO GENÉRICO
# ==============================================================================

# Defina aqui a equação da qual quer encontrar a raiz
def equacao_exemplo(x):
    # Exemplo: Encontrar a raiz de f(x) = x^3 - x - 2
    return x**3 - x - 2

# Limites do intervalo inicial (escolhidos de forma a que f(a) e f(b) tenham sinais opostos)
limite_a = 1.0
limite_b = 2.0

print("INICIANDO O MÉTODO DA FALSA POSIÇÃO...")
raiz_estimada = falsa_posicao(equacao_exemplo, limite_a, limite_b, tol=1e-6)

if raiz_estimada is not None:
    print(f"A raiz da função é aproximadamente: x = {raiz_estimada:.8f}")

if raiz_estimada is not None:
    print(f"A raiz da função é aproximadamente: x = {raiz_estimada:.8f}")
