# Universidade de Brasília (UnB)
## Faculdade de Tecnologia (FT) - Departamento de Engenharia Mecânica (ENM)
**Disciplina:** ENM0227 - Cálculo Numérico Aplicado
**Semestre:** 2026/2
**Professor:** Rafael Gabler Gontijo
**Discente:** Esther Barros Gomes

### Motivação
Este repositório armazena as implementações computacionais desenvolvidas ao longo da disciplina. O objetivo central é a transição rigorosa entre a fundamentação matemática teórica e a construção de algoritmos numéricos eficientes para a modelagem e solução de problemas práticos de engenharia. O repositório funciona como uma apostila didática digital dos métodos numéricos codificados.

### Topologia do Repositório
O repositório adota um particionamento lógico estrito para segregar os Programas para Casa (PPCs) avaliativos das implementações de treino e exercícios de fixação. Cada diretório de prática possui sua própria documentação interna (`README.md`), rotinas executáveis em Python e os respectivos dados de saída:

```text
Calculo-Numerico-Aplicado/
├── README.md (Metadados globais, motivação e topologia da disciplina)
├── Treino_Raizes_de_Funcoes/
│   ├── README.md (Dicionário de variáveis, I/O e escopo do método)
│   ├── bissecao.py (Algoritmo numérico do método da bissecção)
│   ├── problema_canal.py (Aplicação prática do método numérico)
│   └── outputs/
│       └── resultados_bissecao.txt (Resultados numéricos exportados)
└── PPC_01/
    ├── README.md (Documentação e fundamentação da atividade avaliativa)
    ├── ppc1_esther_gomes.py (Script de processamento principal)
    └── outputs/
        ├── grafico_1_stokes_St.png
        ├── grafico_2_refinamento_h.png
        ├── grafico_3_4_validacao_oseen.png
        └── grafico_5_varios_Re.png
```

### Bibliografia Geral
1. CHAPRA, S. C.; CANALE, R. P. Métodos Numéricos para Engenharia. 5. ed. McGraw Hill, 2008.
