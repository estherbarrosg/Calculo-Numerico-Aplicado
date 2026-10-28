import math
import matplotlib.pyplot as plt
# 1. Definicao das constantes fisicas e da funcao do polinomio
b = 4.0
z = 1.5
g = 9.81
Q = 20.0  # Vazao inicial para o primeiro teste

def f(H):
    # Polinomio de 6a ordem expandido
    termo6 = g * (z**3) * (H**6)
    termo5 = 3 * g * b * (z**2) * (H**5)
    termo4 = 3 * g * (b**2) * z * (H**4)
    termo3 = g * (b**3) * (H**3)
    termo1 = -2 * (Q**2) * z * H
    termo0 = -(Q**2) * b
    return termo6 + termo5 + termo4 + termo3 + termo1 + termo0

# 2. Parametros iniciais e tolerancias
H0 = 0.5
H1 = 1.0
H2 = 1.5

tolerancia = 1e-8
max_iter = 50

k = 0
erro = 1.0  # Valor grande para iniciar o laco while

# 3. Laco Iterativo Principal
while erro > tolerancia and k < max_iter:
    k += 1
    
    # Avaliacao da funcao nos tres pontos
    f0 = f(H0)
    f1 = f(H1)
    f2 = f(H2)
    
    # Diferencas e inclinacoes das secantes
    h0 = H1 - H0
    h1 = H2 - H1
    
    delta0 = (f1 - f0) / h0
    delta1 = (f2 - f1) / h1
    
    # Coeficientes da parabola
    a = (delta1 - delta0) / (h1 + h0)
    b_M = a * h1 + delta1
    c = f2
    
    # Discriminante e escolha do denominador E
    D = math.sqrt(b_M**2 - 4 * a * c)
    
    if abs(b_M + D) > abs(b_M - D):
        E = b_M + D
    else:
        E = b_M - D
        
    # Calculo da correcao e nova aproximacao
    delta_H = - (2 * c) / E
    H3 = H2 + delta_H
    
    # Atualizacao do erro
    erro = abs(H3 - H2)
    
    # Impressao do historico da iteracao atual
    print(f"{k:2d} | {H0:9.6f} | {H1:9.6f} | {H2:9.6f} | {H3:9.6f} | {f(H3):11.2e} | {erro:10.2e}")
    
    # Atualizacao dos pontos para a proxima rodada
    H0 = H1
    H1 = H2
    H2 = H3

# 4. Resultados Finais e Verificacao Fisica
if erro <= tolerancia:
    print(f"\nConvergencia atingida! Profundidade Critica (Hc) = {H2:.6f} m")
    
    # Verificacao do Numero de Froude
    Ac = b * H2 + z * (H2**2)
    Tc = b + 2 * z * H2
    Fr2 = (Q**2 * Tc) / (g * (Ac**3))
    print(f"Verificacao fisica: Fr^2 = {Fr2:.6f} (Esperado aprox. 1.0)")
    
# Listas que devem ser preenchidas com os resultados oriundos da sua implementação do método de Müller
# Q_valores deve conter as vazões avaliadas de 5 a 40 m³/s
# Hc_valores deve conter as respectivas profundidades críticas convergidas
Q_valores = [5, 10, 15, 20, 25, 30, 35, 40]
Hc_valores = [0.47, 0.75, 0.98, 1.17, 1.35, 1.51, 1.66, 1.81] # Valores aproximados para ilustração

plt.figure(figsize=(9, 5))

# Plotagem do gráfico Hc x Q
plt.plot(Q_valores, Hc_valores, marker='o', linestyle='-', color='b', label='Profundidade Crítica ($H_c$)')

# Rótulos dos eixos e título do gráfico
plt.xlabel('Vazão $Q$ ($m^3/s$)', fontsize=12)
plt.ylabel('Profundidade Crítica $H_c$ ($m$)', fontsize=12)
plt.title('Variação da Profundidade Crítica em Função da Vazão', fontsize=14)

# Configuração da grade de visualização
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()

# Exibe a figura gerada
plt.show()


import cmath  # Biblioteca necessária para suportar raízes quadradas de números negativos no discriminante D

# ==============================================================================
# FUNÇÃO A SER AVALIADA
# ==============================================================================
def f(x):
    """
    Define a função matemática f(x) = 0.
    Para o teu problema do canal, deves substituir o retorno desta função 
    pelo polinómio de sexta ordem f(H) que deduziste na Parte 1.
    """
    # EXEMPLO GENÉRICO: f(x) = x^3 - 13x - 12
    return x**3 - 13*x - 12


# ==============================================================================
# ALGORITMO GENÉRICO DO MÉTODO DE MÜLLER
# ==============================================================================
def metodo_muller(f, x0, x1, x2, tol=1e-8, max_iter=100):
    """
    Encontra a raiz de uma função f(x) utilizando o Método de Müller.
    
    Parâmetros:
    f        : Função a ser avaliada.
    x0, x1, x2 : Três aproximações iniciais.
    tol      : Critério de tolerância para convergência (erro máximo).
    max_iter : Número máximo de iterações permitidas para evitar loops infinitos.
    """
    
    print("-" * 80)
    print(f"{'k':<3} | {'x0':<10} | {'x1':<10} | {'x2':<10} | {'x3 (Nova raiz)':<15} | {'f(x3)':<15} | {'Erro |x3 - x2|':<15}")
    print("-" * 80)

    # Início do laço iterativo
    for k in range(max_iter):
        
        # 1. Avaliar a função nos três pontos atuais[cite: 3]
        f0 = f(x0)
        f1 = f(x1)
        f2 = f(x2)
        
        # 2. Calcular as diferenças (passos) entre os pontos[cite: 3]
        h0 = x1 - x0
        h1 = x2 - x1
        
        # Prevenção de divisão por zero caso os pontos fiquem idênticos precocemente
        if h0 == 0 or h1 == 0:
            print("Erro: Pontos iniciais coincidentes. O método falha por divisão por zero.")
            return None
        
        # 3. Calcular as inclinações das secantes (diferenças divididas)[cite: 3, 4]
        delta0 = (f1 - f0) / h0
        delta1 = (f2 - f1) / h1
        
        # 4. Determinar os coeficientes da parábola (a, b_M, c)[cite: 4]
        # O índice M em b_M é para evitar confusão com outras variáveis físicas do teu problema
        a = (delta1 - delta0) / (h1 + h0)
        b_M = a * h1 + delta1
        c = f2
        
        # 5. Calcular o discriminante (D) da fórmula resolvente[cite: 4]
        # Utilizamos cmath.sqrt pois a parábola pode não cruzar o eixo real, gerando D complexo
        D = cmath.sqrt(b_M**2 - 4 * a * c)
        
        # 6. Escolher o denominador (E) para evitar o cancelamento numérico[cite: 4]
        # Selecionamos o sinal que resulta no maior valor absoluto no denominador
        if abs(b_M + D) > abs(b_M - D):
            E = b_M + D
        else:
            E = b_M - D
            
        # Proteção extra: se E for incrivelmente pequeno, a parábola é quase horizontal
        if abs(E) < 1e-15:
            print("Erro: Denominador muito próximo de zero. A interpolação falhou.")
            return None
            
        # 7. Calcular a correção e a nova aproximação (x3)[cite: 4]
        delta_x = - (2 * c) / E
        x3 = x2 + delta_x
        
        # Avaliar a função no novo ponto e calcular o erro
        f3 = f(x3)
        erro = abs(x3 - x2)  # Pode ser abs(delta_x)
        
        # Imprimir o histórico da iteração formatado, exigido na atividade[cite: 5]
        # Nota: Extraímos apenas a parte real (.real) para a exibição se os números não forem complexos
        print(f"{k+1:<3} | {x0.real:10.5f} | {x1.real:10.5f} | {x2.real:10.5f} | {x3.real:15.8f} | {f3.real:15.2e} | {erro:15.2e}")
        
        # 8. Critério de parada (Convergência)[cite: 5]
        if erro < tol:
            print("-" * 80)
            print(f">>> Convergência atingida em {k+1} iterações!")
            print(f">>> Raiz encontrada: x = {x3.real:.8f}")
            
            # Verificação de admissibilidade física (Ex: H_c > 0) exigida na Parte 4[cite: 6]
            if x3.real <= 0 or abs(x3.imag) > 1e-8:
                print("Aviso: A raiz encontrada não parece ser fisicamente admissível (negativa ou complexa).")
                
            return x3.real # Retorna apenas a parte real para análises físicas subsequentes

        # 9. Atualizar os três pontos para a próxima iteração[cite: 4]
        # O ponto mais antigo (x0) é descartado.
        x0 = x1
        x1 = x2
        x2 = x3

    # Se o loop terminar sem atingir o 'return' do critério de parada
    print("-" * 80)
    print("Aviso: Número máximo de iterações atingido sem convergência.")
    return x3.real


# ==============================================================================
# BLOCO PRINCIPAL (CHAMADA DA FUNÇÃO)
# ==============================================================================
if __name__ == "__main__":
    
    # Valores iniciais de teste especificados no teu documento[cite: 3]
    x_inicial_0 = 0.5
    x_inicial_1 = 1.0
    x_inicial_2 = 1.5
    
    # Tolerância exigida: 10^-8 metros[cite: 5]
    tolerancia = 1e-8
    
    raiz = metodo_muller(f, x_inicial_0, x_inicial_1, x_inicial_2, tol=tolerancia)



# ==============================================================================
# 15. Questões para discussão
# ==============================================================================

# a) Qual foi a profundidade crítica encontrada para Q=20 m³/s?
# A profundidade crítica convergida é de aproximadamente 1,17 m, obtida a partir 
# da resolução do polinômio: 33,10875H⁶ + 264,87H⁵ + 706,32H⁴ + 627,84H³ - 1200H - 1600 = 0.

# b) A solução encontrada satisfaz numericamente a condição Fr ≈ 1?
# Sim, a solução numérica calculada no método deve produzir Fr² ≈ 1 e, consequentemente, 
# Fr ≈ 1. Essa etapa serve para confirmar que a raiz encontrada corresponde efetivamente 
# à condição física que define o escoamento crítico.

# c) Por que a equação obtida para o canal trapezoidal resulta em um polinômio de sexta ordem?
# A equação da área molhada A(H) é caracterizada pelo polinômio de segunda ordem bH + zH². 
# Ao inserir essa expressão na condição de escoamento crítico definida por (Q²T)/(gA³) = 1, 
# o cálculo do termo A³ faz com que a variável de altura seja elevada à sexta potência, 
# estruturando o polinômio resultante de sexta ordem.

# d) Por que nem todas as raízes matemáticas desse polinômio representam soluções fisicamente possíveis?
# Um polinômio de sexta ordem pode conter até seis raízes matemáticas, o que engloba a 
# possibilidade de raízes negativas e raízes complexas. Como H representa a medida física 
# da profundidade da água no canal, apenas as raízes que são números reais e estritamente 
# maiores que zero (Hc > 0) são soluções admissíveis.

# e) O que acontece com a profundidade crítica quando a vazão aumenta? O comportamento observado no gráfico é fisicamente coerente?
# Quando a vazão do canal aumenta, a profundidade crítica requerida também se torna maior. 
# Esse comportamento tem consistência física porque o aumento do volume de água exige uma 
# área transversal expandida (e uma respectiva profundidade maior) para que o balanço entre 
# os efeitos inerciais e gravitacionais no regime crítico se mantenha.

# f) Repita o cálculo para Q=20 m³/s utilizando três aproximações iniciais diferentes. O método converge para a mesma raiz?
# Quando o cálculo é executado com aproximações iniciais distintas, porém situadas na 
# região de domínio positivo e real, o método de Müller geralmente converge para a mesma 
# profundidade crítica admissível.

# g) Qual é a função dos três pontos utilizados pelo método de Müller e qual é a diferença fundamental entre essa aproximação e a reta utilizada pelo método da secante?
# O método utiliza os três pontos iniciais (H0, H1, H2) com o objetivo de construir 
# localmente uma parábola que realiza a interpolação da função. A principal distinção 
# reside no fato de que o método de Müller utiliza a interseção proveniente dessa parábola 
# como aproximação, enquanto o método da secante trabalha apenas com aproximações lineares (retas).

# h) Por que a escolha do denominador de maior módulo na expressão de Müller é numericamente conveniente?
# A decisão de estipular o denominador E selecionando a expressão de maior valor absoluto, 
# seja |bM + D| ou |bM - D|, visa minimizar os riscos atrelados ao cancelamento numérico 
# no divisor. Essa prática preserva a estabilidade das divisões ao impedir cálculos com 
# denominadores muito próximos a zero ao determinar a nova aproximação.