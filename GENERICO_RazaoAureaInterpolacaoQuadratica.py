import math

def razao_aurea(f, a, b, tol=1e-6, max_iter=100, minimizacao=True):
    """
    Encontra o ponto ótimo (mínimo ou máximo) de uma função unimodal f(x)
    no intervalo [a, b] usando o Método da Razão Áurea.
    
    f           : Função objetivo
    a, b        : Limites do intervalo inicial
    tol         : Tolerância para o critério de paragem (largura do intervalo)
    minimizacao : Booleano (True para achar o Mínimo, False para o Máximo)
    """
    # 1. Definição da constante da Razão Áurea (phi)
    # O valor (sqrt(5) - 1) / 2 equivale a aproximadamente 0.618.
    # O corte nesta proporção exata é o que garante que possamos reaproveitar 
    # pontos de avaliações passadas e poupar esforço computacional a cada ciclo.
    phi = (math.sqrt(5.0) - 1.0) / 2.0
    
    # 2. Cálculo dos dois pontos interiores iniciais (x1 e x2)
    # x1 fica do lado esquerdo do intervalo e x2 do lado direito
    x1 = b - phi * (b - a)
    x2 = a + phi * (b - a)
    
    # 3. Primeira avaliação "pesada" da função
    # Só calculamos a função f() do zero para pontos novos
    f1 = f(x1)
    f2 = f(x2)
    
    print(f"{'Iter':>4} | {'a':>10} | {'b':>10} | {'x_otimo':>10} | {'f(x_otimo)':>14} | {'Erro (b-a)':>14}")
    print("-" * 75)
    
    # 4. Ciclo Iterativo
    # O ciclo repete-se enquanto o "funil" do intervalo for maior que a tolerância
    for i in range(1, max_iter + 1):
        
        # A estimativa atual é sempre declarada como o meio do intervalo restante
        x_atual = (a + b) / 2.0
        f_atual = f(x_atual)
        largura = b - a
        
        print(f"{i:4d} | {a:10.6f} | {b:10.6f} | {x_atual:10.6f} | {f_atual:14.6e} | {largura:14.6e}")
        
        # Critério de Paragem
        if largura < tol:
            print("-" * 75)
            print(f"Convergência atingida com sucesso na iteração {i}!")
            return x_atual
        
        # 5. Lógica de Atualização do Intervalo [a, b]
        # Aqui decidimos qual lado do intervalo será descartado.
        if minimizacao:
            # Se for um problema de Custo, Energia, Perdas, etc.
            condicao = f1 < f2 
        else:
            # Se for um problema de Lucro, Rendimento, Resistência, etc.
            condicao = f1 > f2 
            
        if condicao:
            # O ponto ótimo está "preso" no subintervalo esquerdo: [a, x2]
            b = x2                 # O limite direito recua até x2
            
            # MAGIA DA RAZÃO ÁUREA: O antigo x1 passa a ser o novo x2
            # Copiamos o valor de f1 para f2 sem mandar o computador recalcular nada!
            x2 = x1
            f2 = f1
            
            # Só precisamos de calcular 1 ponto novo: o novo x1
            x1 = b - phi * (b - a)
            f1 = f(x1)
            
        else:
            # O ponto ótimo está "preso" no subintervalo direito: [x1, b]
            a = x1                 # O limite esquerdo avança até x1
            
            # MAGIA DA RAZÃO ÁUREA: O antigo x2 passa a ser o novo x1
            # Copiamos o valor sem recalcular
            x1 = x2
            f1 = f2
            
            # Só precisamos de calcular 1 ponto novo: o novo x2
            x2 = a + phi * (b - a)
            f2 = f(x2)
            
    print("\nAviso: Número máximo de iterações atingido sem convergir totalmente.")
    return (a + b) / 2.0


# ==============================================================================
# EXEMPLO DE USO GENÉRICO
# ==============================================================================

# Defina a equação física/matemática que pretende otimizar
def equacao_exemplo(x):
    # Exemplo: Encontrar o mínimo da função f(x) = (x - 2)^2 + 1
    # Sabemos de antemão que o mínimo físico desta parábola ocorre em x = 2
    return (x - 2)**2 + 1

# Limites do intervalo inicial [a, b]
limite_a = 0.0
limite_b = 5.0

print("INICIANDO O MÉTODO DA RAZÃO ÁUREA...\n")
# Se pretender maximizar a função em vez de minimizar, basta alterar minimizacao=False
ponto_otimo_estimado = razao_aurea(equacao_exemplo, limite_a, limite_b, tol=1e-5, minimizacao=True)

if ponto_otimo_estimado is not None:
    print(f"\nO ponto ótimo da função ocorre em: x = {ponto_otimo_estimado:.8f}")
    
    
    
    
    
import math

def interpolacao_quadratica(f, x0, x1, x2, tol=1e-6, max_iter=100):
    """
    Encontra o ponto de mínimo de uma função unimodal f(x) 
    usando o Método da Interpolação Quadrática.
    
    f           : Função objetivo a ser minimizada
    x0, x1, x2  : Três pontos iniciais (chutes). Devem englobar o mínimo.
    """
    # 1. Ordenação inicial dos pontos
    # Para a lógica de atualização funcionar, precisamos garantir que 
    # os pontos estão rigorosamente em ordem crescente no eixo X (x0 < x1 < x2).
    pontos = sorted([x0, x1, x2])
    x0, x1, x2 = pontos[0], pontos[1], pontos[2]
    
    # Variável de memória para o critério de paragem
    x3_anterior = None
    
    print(f"{'Iter':>4} | {'x0':>8} | {'x1':>8} | {'x2':>8} | {'x3_novo':>10} | {'f(x3)':>12} | {'Erro':>10}")
    print("-" * 75)
    
    # 2. Ciclo Iterativo
    for i in range(1, max_iter + 1):
        
        # Avaliar a função nos três pontos atuais
        f0 = f(x0)
        f1 = f(x1)
        f2 = f(x2)
        
        # 3. Fórmula do Vértice da Parábola (O Coração do Método)
        # O numerador e denominador calculam o ponto exato (x3) onde a 
        # derivada da parábola traçada sobre os 3 pontos é igual a zero.
        num = f0*(x1**2 - x2**2) + f1*(x2**2 - x0**2) + f2*(x0**2 - x1**2)
        den = 2.0 * (f0*(x1 - x2) + f1*(x2 - x0) + f2*(x0 - x1))
        
        # Proteção matemática: Se o denominador for zero, os 3 pontos formam 
        # uma reta perfeitamente horizontal (colineares) e não uma parábola.
        if den == 0:
            print(f"Erro na iteração {i}: Pontos colineares (denominador nulo).")
            return None
            
        # Calcula o novo ponto estimado (o vértice)
        x3 = num / den
        f3 = f(x3)
        
        # 4. Cálculo do Erro e Critério de Paragem
        # O erro é a distância entre o vértice que acabámos de calcular (x3)
        # e o vértice que calculámos na rodada passada (x3_anterior).
        if x3_anterior is not None:
            erro = abs(x3 - x3_anterior)
        else:
            erro = float('inf') # Na primeira iteração o erro é infinito
            
        print(f"{i:4d} | {x0:8.4f} | {x1:8.4f} | {x2:8.4f} | {x3:10.6f} | {f3:12.4e} | {erro:10.4e}")
        
        # Se o vértice parou de se mover, atingimos o fundo da curva real!
        if erro < tol:
            print("-" * 75)
            print(f"Convergência atingida com sucesso na iteração {i}!")
            return x3
            
        x3_anterior = x3
        
        # 5. Lógica de Atualização dos Pontos
        # O objetivo aqui é descartar o pior ponto e manter os 3 pontos 
        # mais próximos do mínimo para desenhar uma nova parábola mais precisa.
        if x3 > x1:
            # O novo ponto caiu do lado DIREITO do ponto central
            if f3 < f1:
                # E o novo ponto é mais baixo (melhor).
                # Descartamos o limite esquerdo x0. O novo trio é [x1, x3, x2]
                x0, x1, x2 = x1, x3, x2
            else:
                # O novo ponto é mais alto (pior). 
                # Descartamos o limite direito x2. O novo trio é [x0, x1, x3]
                x0, x1, x2 = x0, x1, x3
        else: 
            # O novo ponto caiu do lado ESQUERDO do ponto central (x3 < x1)
            if f3 < f1:
                # E o novo ponto é mais baixo (melhor).
                # Descartamos o limite direito x2. O novo trio é [x0, x3, x1]
                x0, x1, x2 = x0, x3, x1
            else:
                # O novo ponto é mais alto (pior).
                # Descartamos o limite esquerdo x0. O novo trio é [x3, x1, x2]
                x0, x1, x2 = x3, x1, x2

    print("\nAviso: Número máximo de iterações atingido sem convergir.")
    return x3


# ==============================================================================
# EXEMPLO DE USO GENÉRICO
# ==============================================================================

# Defina a equação física/matemática que pretende otimizar
def equacao_exemplo(x):
    # Exemplo: Encontrar o mínimo da função f(x) = x^2 - 4x + 4 (cujo mínimo é em x=2)
    return x**2 - 4*x + 4

# Três pontos iniciais quaisquer perto da região de interesse
chute_x0 = 0.0
chute_x1 = 1.0
chute_x2 = 5.0

print("INICIANDO O MÉTODO DA INTERPOLAÇÃO QUADRÁTICA...\n")
ponto_otimo_estimado = interpolacao_quadratica(equacao_exemplo, chute_x0, chute_x1, chute_x2, tol=1e-5)

if ponto_otimo_estimado is not None:
    print(f"\nO ponto ótimo da função ocorre em: x = {ponto_otimo_estimado:.8f}")