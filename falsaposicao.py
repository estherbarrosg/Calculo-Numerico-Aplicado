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
    print("Iniciando o Método da Bissecção")
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