# =============================================================================
# Universidade de Brasília (UnB) - Engenharia Mecânica
# Disciplina: Cálculo Numérico Aplicado (ENM0227)
# Professor: Rafael Gabler Gontijo
# Aluna: Esther Barros Gomes
# 
# Programa para Casa #1 (PPC1)
# Método de Runge-Kutta de 4ª Ordem para Sedimentação de Partícula
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# 1. FUNÇÕES DO PROBLEMA E SOLUÇÕES ANALÍTICAS
# =============================================================================

def f_derivada(t, v, St, Re_s):
    """
    Equação diferencial adimensional do movimento da partícula.
    dv*/dt* = (1/St) * (1 - v* - (3/8)*Re_s*(v*)^2)
    """
    return (1.0 / St) * (1.0 - v - (3.0 / 8.0) * Re_s * (v**2))

def analitica_stokes(t, St):
    """
    Solução analítica para Re -> 0 (Fatores Integrantes).
    """
    return 1.0 - np.exp(-t / St)

def analitica_oseen(t, St, Re_s):
    """
    Solução analítica para Re =/= 0 (Equação de Riccati).
    """
    if Re_s == 0:
        return analitica_stokes(t, St)
        
    v_star = (-1.0 + np.sqrt(1.0 + 1.5 * Re_s)) / (0.75 * Re_s)
    Q = (3.0 / (8.0 * St)) * Re_s
    P = -(3.0 / (4.0 * St)) * Re_s * v_star - (1.0 / St)
    
    # Prevenção de divisão por zero ou overflow
    termo_exponencial = np.exp(-P * t)
    v_z = v_star + 1.0 / ((Q / P) + ((-1.0 / v_star) - (Q / P)) * termo_exponencial)
    return v_z

# =============================================================================
# 2. IMPLEMENTAÇÃO DO MÉTODO DE RUNGE-KUTTA (RK4)
# =============================================================================

def integrador_rk4(St, Re_s, h, t_final):
    """
    Implementação do método RK4 conforme o pseudocódigo planejado.
    Retorna as listas de tempo e velocidade.
    """
    t_atual = 0.0
    v_atual = 0.0
    
    lista_t = [t_atual]
    lista_v = [v_atual]
    
    while t_atual < t_final:
        # Etapas do RK4
        k1 = f_derivada(t_atual, v_atual, St, Re_s)
        k2 = f_derivada(t_atual + 0.5 * h, v_atual + 0.5 * k1 * h, St, Re_s)
        k3 = f_derivada(t_atual + 0.5 * h, v_atual + 0.5 * k2 * h, St, Re_s)
        k4 = f_derivada(t_atual + h, v_atual + k3 * h, St, Re_s)
        
        # Atualização da velocidade e tempo
        v_atual = v_atual + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        t_atual = t_atual + h
        
        # Armazenamento (Append)
        lista_t.append(t_atual)
        lista_v.append(v_atual)
        
    return np.array(lista_t), np.array(lista_v)

# =============================================================================
# 3. ANÁLISES REQUISITADAS
# =============================================================================

def analise_1_stokes_varios_St():
    """Tarefa 1: Comparação Numérico x Analítico para Re -> 0 com diferentes St"""
    print("Gerando Análise 1: Re -> 0 para diferentes St...")
    Re_s = 0.0
    h = 0.1
    t_final = 20.0
    valores_St = [0.5, 1.0, 2.0]
    
    plt.figure(figsize=(10, 6))
    for St in valores_St:
        t_num, v_num = integrador_rk4(St, Re_s, h, t_final)
        v_ana = analitica_stokes(t_num, St)
        
        plt.plot(t_num, v_num, 'o', markevery=10, label=f'RK4 (St={St})')
        plt.plot(t_num, v_ana, '-', label=f'Analítica (St={St})')
        
    plt.title('Sedimentação de Partícula: Numérico x Analítico (Re -> 0)')
    plt.xlabel('Tempo Adimensional ($t^*$)')
    plt.ylabel('Velocidade Adimensional ($v_z^*$)')
    plt.legend()
    plt.grid(True)
    plt.savefig('grafico_1_stokes_St.png')
    plt.show()

def analise_2_refinamento_temporal():
    """Tarefa 2: Variação do passo de tempo (h)"""
    print("Gerando Análise 2: Refinamento Temporal...")
    St = 1.0
    Re_s = 0.0
    t_final = 10.0
    passos_h = [1.0, 0.5, 0.1]
    
    plt.figure(figsize=(10, 6))
    t_ref = np.linspace(0, t_final, 500)
    plt.plot(t_ref, analitica_stokes(t_ref, St), 'k-', label='Solução Exata')
    
    for h in passos_h:
        t_num, v_num = integrador_rk4(St, Re_s, h, t_final)
        plt.plot(t_num, v_num, '--', marker='x', markevery=int(1/h), label=f'RK4 (h={h})')
        
    plt.title('Estudo do Refinamento Temporal (h)')
    plt.xlabel('Tempo Adimensional ($t^*$)')
    plt.ylabel('Velocidade Adimensional ($v_z^*$)')
    plt.legend()
    plt.grid(True)
    plt.savefig('grafico_2_refinamento_h.png')
    plt.show()

def analise_3_4_oseen_validacao():
    """Tarefas 3 e 4: Validação com efeito inercial (Re =/= 0)"""
    print("Gerando Análises 3 e 4: Validação Oseen (Re =/= 0)...")
    St = 1.0
    Re_s = 1.0
    h = 0.05
    t_final = 10.0
    
    t_num, v_num = integrador_rk4(St, Re_s, h, t_final)
    v_ana = analitica_oseen(t_num, St, Re_s)
    
    plt.figure(figsize=(10, 6))
    plt.plot(t_num, v_num, 'ro', markevery=5, label=f'Numérico RK4 (Re_s={Re_s})')
    plt.plot(t_num, v_ana, 'b-', label=f'Analítica Riccati (Re_s={Re_s})')
    
    plt.title('Validação do Modelo com Inércia ($Re_s \neq 0$)')
    plt.xlabel('Tempo Adimensional ($t^*$)')
    plt.ylabel('Velocidade Adimensional ($v_z^*$)')
    plt.legend()
    plt.grid(True)
    plt.savefig('grafico_3_4_validacao_oseen.png')
    plt.show()

def analise_5_diferentes_reynolds():
    """Tarefa 5: Comportamento para diferentes valores de Re_s"""
    print("Gerando Análise 5: Comparação de diferentes números de Reynolds...")
    St = 1.0
    h = 0.05
    t_final = 15.0
    valores_Re = [0.0, 0.5, 1.0, 5.0]
    
    plt.figure(figsize=(10, 6))
    for Re_s in valores_Re:
        t_num, v_num = integrador_rk4(St, Re_s, h, t_final)
        plt.plot(t_num, v_num, '-', linewidth=2, label=f'Re_s = {Re_s}')
        
    plt.title('Desvio do Regime Assintótico para diferentes $Re_s$')
    plt.xlabel('Tempo Adimensional ($t^*$)')
    plt.ylabel('Velocidade Adimensional ($v_z^*$)')
    plt.legend()
    plt.grid(True)
    plt.savefig('grafico_5_varios_Re.png')
    plt.show()

# =============================================================================
# BLOCO PRINCIPAL
# =============================================================================
if __name__ == "__main__":
    print("Iniciando bateria de análises - PPC1...")
    analise_1_stokes_varios_St()
    analise_2_refinamento_temporal()
    analise_3_4_oseen_validacao()
    analise_5_diferentes_reynolds()
    print("Análises concluídas! Gráficos salvos no diretório atual.")
    