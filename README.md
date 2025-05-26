Trabalho 2: Otimização de Manutenção para Sistemas Reparáveis

### Introdução

Este projeto aborda a análise de confiabilidade e a otimização de políticas de manutenção para sistemas reparáveis, especificamente aplicando o modelo de **Processo de Lei de Potência (PLP)** para sistemas que sofrem **Reparo Mínimo (MR)** após falhas e Manutenção Preventiva (PM) perfeita em intervalos programados. A metodologia e os dados base são inspirados no artigo "Optimal Maintenance Time for Repairable Systems" de Gilardoni & Colosimo (2007) e nas planilhas fornecidas pelo professor, que utilizam dados de falha de transformadores de potência.

O objetivo principal é determinar os parâmetros do modelo PLP a partir dos dados históricos de falha e, subsequentemente, calcular o intervalo ótimo para a realização de manutenções preventivas que minimize o custo total esperado por unidade de tempo.

### Estrutura do Trabalho e Metodologia

O trabalho foi dividido em duas partes principais, implementadas em R:

**Parte a: Estimação dos Parâmetros do Processo de Lei de Potência (PLP) e Avaliação do Modelo**

1.  **Modelo de Falha:**
    *   Assume-se que o processo de falha dos transformadores entre manutenções preventivas (que renovam o sistema) segue um Processo de Poisson Não Homogêneo (NHPP), especificamente o Processo de Lei de Potência (PLP).
    *   A função de intensidade de falha (ROCOF - Rate of Occurrence of Failures) para o PLP é dada por:
        ρ(t | β, θ) = (β/θ) \* (t/θ)^(β-1)
    *   A função de intensidade acumulada (MCF - Mean Cumulative Function) ou número esperado de falhas até o tempo `t` é:
        M(t | β, θ) = (t/θ)^β
    *   Os parâmetros a serem estimados são:
        *   **β (Beta):** Parâmetro de forma do PLP. Indica a natureza da taxa de falha (β > 1: crescente/desgaste; β = 1: constante/HPP; β < 1: decrescente).
        *   **θ (Theta):** Parâmetro de escala do PLP.

2.  **Preparação dos Dados:**
    *   Os dados de entrada consistem em uma tabela processada (derivada da Tabela 1 do artigo de Gilardoni & Colosimo), onde cada linha representa um "período de operação".
    *   Cada período é definido pelo seu tempo total de observação (Tᵢ) e um indicador se o período terminou com uma falha em Tᵢ (`Event_final = 1`) ou com uma censura/PM em Tᵢ (`Event_final = 0`).
    *   *Nota: Para esta análise, foi utilizada uma função de verossimilhança simplificada que considera apenas o status no final de cada período Tᵢ, não incluindo explicitamente falhas que possam ter ocorrido *antes* de Tᵢ dentro do mesmo período.*

3.  **Estimação por Máxima Verossimilhança (MLE):**
    *   Uma função de log-verossimilhança para o PLP (baseada na interpretação simplificada dos dados de entrada) foi implementada.
    *   A função `optim()` do R foi utilizada para encontrar as estimativas de máxima verossimilhança (β̂ e θ̂) que maximizam a função de log-verossimilhança.

4.  **Avaliação da Adequação do Modelo:**
    *   O valor do parâmetro de forma estimado (β̂) é analisado: um β̂ > 1 sugere uma taxa de falha crescente, o que é consistente com a premissa de envelhecimento do sistema e justifica a manutenção preventiva.
    *   Um gráfico comparativo é gerado, mostrando:
        *   A **MCF Paramétrica** estimada: M̂(t) = (t/θ̂)^β̂ (linha vermelha).
        *   A **MCF Não-Paramétrica** (Nelson-Aalen), calculada de forma agregada a partir dos dados dos períodos de observação (pontos/degraus azuis).
    *   A proximidade visual entre a MCF paramétrica e a não-paramétrica ajuda a avaliar a adequação do modelo PLP ajustado.
    *   O Critério de Informação de Akaike (AIC) também é calculado para o modelo ajustado.

**Parte b: Determinação da Política Ótima de Manutenção Preventiva**

1.  **Definição de Custos:**
    *   **CPM:** Custo de uma Manutenção Preventiva (PM) perfeita (que renova o sistema).
    *   **CMR:** Custo de um Reparo Mínimo (MR) realizado após uma falha.
    *   Para este trabalho, foram utilizados os valores: CPM = 1 e CMR = 15 unidades monetárias.

2.  **Função de Custo Esperado por Unidade de Tempo (H(τ)):**
    *   O objetivo é minimizar o custo esperado a longo prazo por unidade de tempo, H(τ), onde τ é o intervalo entre as PMs.
    *   A função de custo para o PLP é:
        H(τ) = (1/τ) \* \[CPM + CMR \* (τ/θ̂)^β̂]

3.  **Cálculo do Intervalo Ótimo de PM (τ\_ótimo):**
    *   O intervalo ótimo τ\_ótimo que minimiza H(τ) é calculado usando a fórmula derivada no artigo de Gilardoni & Colosimo (Equação 5):
        τ\_ótimo = θ̂ \* \[CPM / ((β̂-1) \* CMR)]^(1/β̂)
    *   Este cálculo é válido se β̂ > 1.

4.  **Visualização:**
    *   A função de custo H(τ) é plotada em função de τ, focando na região próxima ao τ\_ótimo calculado, para visualizar o ponto mínimo da curva de custo.

### Implementação em R

O código R para este projeto realiza as seguintes etapas:
1.  Define o dataframe com os dados dos períodos de operação dos transformadores.
2.  Implementa a função de log-verossimilhança para o PLP (abordagem simplificada).
3.  Utiliza `optim()` para estimar os parâmetros β e θ do PLP.
4.  Calcula e plota a MCF paramétrica e a MCF não-paramétrica (Nelson-Aalen).
5.  Calcula o intervalo ótimo de manutenção preventiva (τ\_ótimo) com base nos parâmetros estimados e nos custos fornecidos.
6.  Plota a função de custo esperado por unidade de tempo H(τ).

### Resultados e Conclusões (Exemplo de como preencher)

*   **Parâmetros Estimados do PLP (Letra a):**
    *   β̂ (Forma) ≈ \[Seu valor de beta_hat_plp_a]
    *   θ̂ (Escala) ≈ \[Seu valor de theta_hat_plp_a]
    *   AIC ≈ \[Seu valor de aic_value_plp_a]
*   **Adequação do Modelo (Letra a):**
    *   O valor de β̂ \[maior que/próximo de/menor que] 1 sugere uma taxa de falha \[crescente/constante/decrescente].
    *   A análise visual do gráfico da MCF \[indica um bom ajuste/mostra algumas discrepâncias] entre o modelo PLP paramétrico e a estimativa não-paramétrica.
*   **Política Ótima de Manutenção (Letra b):**
    *   τ\_ótimo ≈ \[Seu valor de tau_otimo em horas] horas (aproximadamente \[Seu valor em dias] dias).
    *   O gráfico da função de custo H(τ) confirma visualmente este ponto ótimo.

Este trabalho demonstra a aplicação de modelos estatísticos para otimizar decisões de manutenção em sistemas reparáveis, visando minimizar custos a longo prazo.