# Cálculo Numérico Aplicado

Repositório de códigos da disciplina, com as atividades práticas e os Programas para Casa (PPC) desenvolvidos ao longo do semestre.

| | |
|---|---|
| **Instituição** | Universidade de Brasília (UnB) |
| **Faculdade** | Faculdade de Tecnologia |
| **Departamento** | Departamento de Engenharia Mecânica |
| **Disciplina** | ENM0227 - Cálculo Numérico Aplicado |
| **Semestre** | 2026/2 |
| **Docente** | Prof. Rafael Gabler Gontijo |
| **Discente** | Esther Barros Gomes |

## Por que estudar métodos numéricos?

Poucos problemas reais têm solução analítica. Equações como a de Manning para escoamento em canais ou o balanço térmico com radiação (que envolve $T^4$) não podem ser resolvidas isolando a incógnita com álgebra. Os métodos numéricos substituem a solução exata por uma sequência de aproximações que se aproximam dela, com erro controlado. Estudar esses métodos "por dentro", escrevendo cada algoritmo, ensina não só a obter uma resposta, mas a julgar se ela é confiável: quão rápido o método converge, quanto custa cada iteração e em que condições ele falha.

Este repositório reúne essas implementações, cada uma acompanhada de documentação com o problema, a fundamentação matemática, o algoritmo e a validação dos resultados.

## Topologia do repositório

```
Calculo-Numerico-Aplicado/
├── README.md                        Este arquivo (README externo)
├── .gitignore                       Arquivos que o Git deve ignorar (cache do Python)
├── Atividades_Praticas/             Atividades de prática propostas em aula (não são PPCs)
│   ├── 01_Bissecao_Falsa_Posicao/   Raízes de funções por métodos intervalares
│   └── 02_Newton_Secante/           Raízes de funções por métodos abertos
└── PPC_XX/                          Um diretório por Programa para Casa (PPC_01, PPC_02, ...),
                                     adicionados à medida que cada PPC é entregue
```

Cada diretório de atividade ou PPC é isolado e contém o código-fonte, um `README.md` interno (problema, matemática, dicionário de variáveis, execução e validação) e a pasta `outputs/` com os resultados e figuras gerados.

### Atividades práticas

| Pasta | Problema físico | Métodos | Resultado principal |
|---|---|---|---|
| [`01_Bissecao_Falsa_Posicao`](Atividades_Praticas/01_Bissecao_Falsa_Posicao) | Profundidade da lâmina d'água em um canal aberto (equação de Manning) | Bisseção e falsa posição | $H = 0{,}702293$ m; 22 e 7 iterações |
| [`02_Newton_Secante`](Atividades_Praticas/02_Newton_Secante) | Temperatura de equilíbrio de uma superfície com convecção e radiação | Newton-Raphson e secante | $T = 409{,}361573$ K; 5 e 5 iterações |

### Guia rápido: onde está cada método

| Método | Arquivo | Função | Precisa de |
|---|---|---|---|
| Bisseção | `Atividades_Praticas/01_Bissecao_Falsa_Posicao/metodos_intervalares.py` | `bissecao(funcao, limite_inferior, limite_superior, tolerancia, max_iteracoes)` | intervalo com mudança de sinal |
| Falsa posição | `Atividades_Praticas/01_Bissecao_Falsa_Posicao/metodos_intervalares.py` | `falsa_posicao(funcao, limite_inferior, limite_superior, tolerancia, max_iteracoes)` | intervalo com mudança de sinal |
| Newton-Raphson | `Atividades_Praticas/02_Newton_Secante/metodos_abertos.py` | `newton_raphson(funcao, derivada, chute_inicial, tolerancia, max_iteracoes)` | chute inicial e derivada analítica |
| Secante | `Atividades_Praticas/02_Newton_Secante/metodos_abertos.py` | `secante(funcao, chute_anterior, chute_atual, tolerancia, max_iteracoes)` | dois chutes iniciais |

Todos os métodos são genéricos (servem para qualquer $f(x)$) e devolvem três valores: a aproximação final, o histórico de iterações e um indicador de convergência.

## Como executar

Requisitos:

- **Python 3** (testado com Python 3.12). Os métodos e as simulações usam apenas a biblioteca padrão.
- **matplotlib**, somente para os scripts `graficos.py` que geram as figuras (`pip install matplotlib`).

Cada atividade é executada de dentro da sua pasta (ou da raiz, informando o caminho):

```bash
cd Atividades_Praticas/01_Bissecao_Falsa_Posicao
python main.py        # simulação e relatório em outputs/resultados.txt
python graficos.py    # figuras em outputs/
```

As instruções detalhadas de cada programa estão no `README.md` da respectiva pasta.

## Convenções de código

Seguindo o Manual de Boas Práticas da disciplina:

- Somente funcionalidades matemáticas básicas e gráficos (`matplotlib`); nenhuma função pronta que resolva diretamente o problema numérico.
- Nomes de variáveis autoexplicativos, código comentado e sem variáveis não utilizadas.
- Sem repetição de lógica: o problema físico fica em um módulo próprio e os métodos em outro, de modo que o mesmo problema é resolvido por métodos diferentes sem duplicar código.
- Cada método é uma função genérica, separada do problema, e todos os resultados são verificados contra uma referência (solução analítica, dados do enunciado ou outro método) na seção de validação do README de cada pasta.

## Bibliografia geral

- CHAPRA, Steven C.; CANALE, Raymond P. *Métodos numéricos para engenharia*. McGraw-Hill.
- UNIVERSIDADE DE BRASÍLIA. Departamento de Engenharia Mecânica. *Plano de Ensino: ENM0227 - Cálculo Numérico Aplicado*. Semestre 2026/2.
- UNIVERSIDADE DE BRASÍLIA. Departamento de Engenharia Mecânica. *Manual de boas práticas: Cálculo Numérico Aplicado*. Brasília, 2026.
