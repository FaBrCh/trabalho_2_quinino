# Carregar pacotes
library(dplyr)
library(ggplot2)
library(survival)

# Criar o dataframe com os dados do artigo
dados_trafo <- data.frame(
  Equi = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 5, 6, 7, 7, 8, 8, 9, 9, 10, 11, 12, 12, 12, 13, 14, 14, 14, 15, 16, 16, 17, 17, 18, 19, 20, 20, 20, 21, 22, 23, 23, 24, 24, 25, 26, 27, 27, 28, 28, 29, 30, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40),
  Failure_Time = c(8839, 17057, 21887, 9280, 16442, 21887, 10445, 13533, 7902, 8414, 13331, 17156, 21887, 16305, 21887, 16802, 21887, 4881, 16625, 7396, 7541, 19590, 2211, 15821, 19746, 19877, 1927, 15813, 21886, 15524, 21886, 21440, 369, 11664, 17031, 21857, 7544, 6039, 2168, 6698, 18840, 21879, 2288, 2499, 10668, 16838, 15550, 21887, 1616, 14041, 20004, 21888, 21888, 21888, 21888, 21888, 21888, 21888, 21888, 21888, 21888),
  Censura = c(1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
)


## 1. Criar dados em formato de intervalo para o Nelson-Aalen
dados_intervalo <- dados_trafo %>%
  group_by(Equi) %>%
  mutate(
    tempo_inicio = lag(Failure_Time, default = 0),
    tempo_fim = Failure_Time,
    evento = Censura
  ) %>%
  ungroup()

## 2. Criar dados em formato de lista para a Máxima Verossimilhança
dados_lista <- dados_intervalo %>%
  group_by(Equi) %>%
  summarise(
    falhas = list(Failure_Time[evento == 1]),
    tempo_final = max(tempo_fim)
  )


## Função para a log-verossimilhança negativa do Power Law Process
# params[1] = beta (forma), params[2] = theta (escala)
neg_log_likelihood_plp <- function(params, dados) {
  beta <- params[1]
  theta <- params[2]
  
  # Parâmetros devem ser positivos
  if (beta <= 0 || theta <= 0) return(Inf)
  
  # sapply itera sobre cada sistema (cada linha de 'dados_lista')
  total_log_likelihood <- sum(sapply(1:nrow(dados), function(i) {
    tempos_falha <- unlist(dados$falhas[[i]])
    tempo_max <- dados$tempo_final[i]
    
    # Termo da verossimilhança para as falhas observadas
    # log(intensidade) = log(dweibull/pweibull)
    log_intensidade <- 0
    if (length(tempos_falha) > 0) {
      log_intensidade <- sum(dweibull(tempos_falha, shape = beta, scale = theta, log = TRUE) - 
                               pweibull(tempos_falha, shape = beta, scale = theta, lower.tail = FALSE, log.p = TRUE))
    }
    
    # Termo da verossimilhança para o período total de observação
    # Intensidade acumulada M(T) = (T / theta)^beta
    intensidade_acumulada <- (tempo_max / theta)^beta
    
    # Log-verossimilhança para este sistema
    return(log_intensidade - intensidade_acumulada)
  }))
  
  # A função optim minimiza, então retornamos o valor negativo
  return(-total_log_likelihood)
}

## Otimização para encontrar os melhores parâmetros
# Usamos "chutes" iniciais razoáveis para os parâmetros
valores_iniciais <- c(2.0, 20000) 
# Usamos os valores iniciais para definir a escala dos parâmetros
otimizacao_corrigida <- optim(
  par = valores_iniciais,
  fn = neg_log_likelihood_plp,
  dados = dados_lista,
  method = "BFGS",
  # AQUI ESTÁ A CORREÇÃO: informamos a escala dos parâmetros
  control = list(parscale = valores_iniciais) 
)

# Extrair e mostrar os parâmetros agora corretos
beta_estimado_correto <- otimizacao_corrigida$par[1]
theta_estimado_correto <- otimizacao_corrigida$par[2]

cat("✅ Parâmetros da Weibull Corrigidos:\n")
cat("   - Beta (Forma):", beta_estimado_correto, "\n")
cat("   - Theta (Escala):", theta_estimado_correto, "\n")



## 1. Calcular a curva não paramétrica (Nelson-Aalen)
surv_obj <- Surv(time = dados_intervalo$tempo_inicio, 
                 time2 = dados_intervalo$tempo_fim, 
                 event = dados_intervalo$evento)
fit_na <- survfit(surv_obj ~ 1)
mcf_na_data <- data.frame(tempo = fit_na$time, mcf = fit_na$cumhaz)

## 2. Definir a função da curva do nosso modelo (Paramétrica)
mcf_parametrica <- function(t, beta, theta) {
  (t / theta)^beta
}

## Gerar o gráfico comparativo com os parâmetros corretos
ggplot(mcf_na_data, aes(x = tempo, y = mcf)) +
  geom_step(aes(color = "Dados (Nelson-Aalen)"), size = 1.2) +
  stat_function(
    fun = mcf_parametrica,
    # Usando os novos parâmetros
    args = list(beta = beta_estimado_correto, theta = theta_estimado_correto),
    aes(color = "Modelo (Weibull)"),
    size = 1
  ) +
  labs(
    title = "Avaliação do Modelo (Versão Corrigida)",
    subtitle = "Com os parâmetros corretos, o ajuste do modelo é excelente.",
    x = "Tempo (horas)",
    y = "Média de Falhas Acumuladas (MCF)",
    color = "Estimador"
  ) +
  scale_color_manual(values = c("Dados (Nelson-Aalen)" = "dodgerblue", "Modelo (Weibull)" = "red")) +
  theme_bw() +
  coord_cartesian(xlim = c(0, 22000))



C_PM <- 1
C_MR <- 15

# Preparando os dados para o gráfico
intervalo_tau <- seq(tau_otimo - 800, tau_otimo + 800, length.out = 500)
dados_custo <- data.frame(
  tau = intervalo_tau,
  custo = (1 / intervalo_tau) * (C_PM + C_MR * (intervalo_tau / theta_estimado_correto)^beta_estimado_correto)
)

# Criar o gráfico com a posição do texto ajustada
ggplot(dados_custo, aes(x = tau, y = custo)) +
  geom_line(color = "blue", size = 1) +
  geom_vline(xintercept = tau_otimo, color = "red", linetype = "dashed", size = 1) +
  annotate("point", x = tau_otimo, y = custo_otimo, color = "red", size = 4) +
  
  # --- INÍCIO DA ALTERAÇÃO ---
  # Diminuímos o valor somado a 'y' para aproximar o texto do ponto
  annotate(
    "text",
    x = tau_otimo,
    y = custo_otimo + 0.000002, # <-- VALOR AJUSTADO AQUI
    label = paste0("T_ótimo ≈ ", round(tau_otimo), "\nH(T_ótimo) ≈ ", format(custo_otimo, scientific = FALSE, digits = 5)),
    hjust = 0.5,
    vjust = 0 # Alinha a base do texto com a coordenada y
  ) +
  # --- FIM DA ALTERAÇÃO ---
  
  labs(
    title = "Função de Custo da Manutenção Preventiva (Focado no Ótimo)",
    x = "Intervalo τ (horas)",
    y = "H(τ) - Custo por Unidade de Tempo"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_line(linetype = 'dotted'),
        panel.grid.minor = element_line(linetype = 'dotted'))

