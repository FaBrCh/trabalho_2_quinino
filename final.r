# --- 1. Criação do Dataframe (Dados Corrigidos - 61 observações) ---
# Estes dados já representam os períodos de observação T_i e o status no final do período.
Equi <- c(1, 1, 1, 2, 2, 2, 3, 3, 4, 5, 6, 7, 7, 8, 8, 9, 9, 10, 11, 12, 
          12, 12, 13, 14, 14, 14, 15, 16, 16, 17, 17, 18, 19, 20, 20, 20, 
          21, 22, 23, 23, 24, 24, 25, 26, 27, 27, 28, 28, 29, 30, 30, 31, 
          32, 33, 34, 35, 36, 37, 38, 39, 40)
Time_T_i <- c(8839, 17057, 21887, 9280, 16442, 21887, 10445, 13533, 7902, 
                 8414, 13331, 17156, 21887, 16305, 21887, 16802, 21887, 4881, 
                 16625, 7396, 7541, 19590, 2211, 15821, 19746, 19877, 1927, 
                 15813, 21886, 15524, 21886, 21440, 369, 11664, 17031, 21857, 
                 7544, 6039, 2168, 6698, 18840, 21879, 2288, 2499, 10668, 
                 16838, 15550, 21887, 1616, 14041, 20004, 21888, 21888, 21888, 
                 21888, 21888, 21888, 21888, 21888, 21888, 21888)
Event_at_T_i <- c(1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 
             0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 
             0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

df_plp_periods <- data.frame(Equi = Equi, T_obs = Time_T_i, Event_final = Event_at_T_i)

#===============================================================================
# LETRA A: Encontrar os parâmetros da "Weibull" (β e θ do PLP) e avaliar adequação
#===============================================================================

# --- 2. Função de Log-Verossimilhança para PLP (Simplificada) ---
# Assume que se Event_final = 1, a falha ocorreu em T_obs (contribui com rho(T_obs)).
# Se Event_final = 0, o período foi censurado/PM em T_obs.
# Todos os períodos contribuem com -M(T_obs).
log_likelihood_PLP_simplified <- function(params, T_observations, final_event_status) {
  beta <- params[1]  # Forma do PLP
  theta <- params[2] # Escala do PLP
  
  if (beta <= 1e-9 || theta <= 1e-9) return(1e10) 
  
  log_lik_contributions <- numeric(length(T_observations))
  
  for (i in 1:length(T_observations)) {
    T_i <- T_observations[i]
    event_i <- final_event_status[i]
    
    log_rho_Ti_contrib <- 0 
    if (event_i == 1) { 
      if (T_i <= 1e-9) { 
          log_rho_Ti_contrib <- -1e10 
      } else {
          log_rho_Ti_contrib <- log(beta) - beta * log(theta) + (beta - 1) * log(T_i)
      }
      if (!is.finite(log_rho_Ti_contrib)) log_rho_Ti_contrib <- -1e10
    }
    
    M_Ti_contrib <- 0 
    if (T_i > 1e-9) {
        M_Ti_contrib <- (T_i / theta)^beta
    }
    if (!is.finite(M_Ti_contrib)) M_Ti_contrib <- 1e10
    
    log_lik_contributions[i] <- log_rho_Ti_contrib - M_Ti_contrib
  }
  
  total_log_lik <- sum(log_lik_contributions)
  
  if (!is.finite(total_log_lik)) return(1e10)
  return(-total_log_lik) 
}

# --- 3. Estimação dos Parâmetros PLP ---
# Usar valores iniciais próximos aos do Excel/Artigo
initial_params_plp <- c(beta = 1.99, theta = 24500) 

mle_fit_plp_simp <- NULL
try({
  mle_fit_plp_simp <- optim(
    par = initial_params_plp, 
    fn = log_likelihood_PLP_simplified, 
    T_observations = df_plp_periods$T_obs, 
    final_event_status = df_plp_periods$Event_final,
    method = "L-BFGS-B", 
    lower = c(0.1, 1000.0), 
    upper = c(5.0, 60000.0), 
    control = list(maxit = 1000, factr = 1e7) 
  )
})

# --- 4. Resultados da Estimação PLP e Gráfico da MCF (Letra a) ---
if (!is.null(mle_fit_plp_simp) && mle_fit_plp_simp$convergence == 0) {
  beta_hat_plp_a <- mle_fit_plp_simp$par[1]
  theta_hat_plp_a <- mle_fit_plp_simp$par[2]
  max_log_lik_plp_a <- -mle_fit_plp_simp$value 
  
  cat("--- LETRA A: Resultados da Estimação PLP (Simplificada) ---\n")
  cat("Beta Estimado (Forma do PLP):", format(beta_hat_plp_a, digits = 5), "\n")
  cat("Theta Estimado (Escala do PLP):", format(theta_hat_plp_a, digits = 5), "\n")
  cat("Log-Verossimilhança Maximizada:", format(max_log_lik_plp_a, digits = 5), "\n")
  
  k_params_plp_a <- 2 
  aic_value_plp_a <- -2 * max_log_lik_plp_a + 2 * k_params_plp_a
  cat("AIC (PLP Simplificado):", format(aic_value_plp_a, digits = 5), "\n\n")

  # --- Cálculo da MCF Não-Paramétrica (Nelson-Aalen) ---
  # Usando a abordagem de contar períodos em risco
  observed_failure_times_agg <- sort(unique(df_plp_periods$T_obs[df_plp_periods$Event_final == 1]))
  mcf_na_values_mod <- numeric(length(observed_failure_times_agg))
  current_mcf_na_mod <- 0
  
  if (length(observed_failure_times_agg) > 0) {
    for (i in 1:length(observed_failure_times_agg)) {
      t_k <- observed_failure_times_agg[i]
      d_k <- sum(df_plp_periods$T_obs == t_k & df_plp_periods$Event_final == 1)
      r_k <- sum(df_plp_periods$T_obs >= t_k) 
      if (r_k > 0) {
        current_mcf_na_mod <- current_mcf_na_mod + (d_k / r_k)
      } else if (d_k > 0) { current_mcf_na_mod <- current_mcf_na_mod + d_k }
      mcf_na_values_mod[i] <- current_mcf_na_mod
    }
  }
  df_mcf_na_plot <- data.frame(Time = observed_failure_times_agg, MCF_NA = mcf_na_values_mod)

  # --- Plotagem do Gráfico da MCF (Letra a) ---
  max_time_plot <- max(df_plp_periods$T_obs, na.rm = TRUE)
  tempos_plot_param <- seq(0, max_time_plot + 1000, length.out = 200)
  # MCF Paramétrica do PLP: M(t) = (t/theta)^beta
  mcf_linha_vermelha_plp <- (tempos_plot_param / theta_hat_plp_a)^beta_hat_plp_a
  
  ylim_max_plot <- 0
  if(nrow(df_mcf_na_plot) > 0 && length(df_mcf_na_plot$MCF_NA) > 0) {
      ylim_max_plot <- max(c(mcf_linha_vermelha_plp, df_mcf_na_plot$MCF_NA), na.rm = TRUE) * 1.1
  } else {
      ylim_max_plot <- max(mcf_linha_vermelha_plp, na.rm = TRUE) * 1.1
  }
  if (!is.finite(ylim_max_plot) || ylim_max_plot <= 0) ylim_max_plot = max(mcf_linha_vermelha_plp, na.rm=TRUE) * 1.2 
  if (!is.finite(ylim_max_plot) || ylim_max_plot <= 0) ylim_max_plot = 1.0 

  plot(tempos_plot_param, mcf_linha_vermelha_plp, type = "l", col = "red", lwd = 2,
       xlab = "Time (hours)", ylab = "MCF", 
       xlim = c(0, max_time_plot + 1000),
       ylim = c(0, ylim_max_plot), 
       main = "Mean Cumulative Function (PLP vs. Nelson-Aalen)")
  grid()
  
  if (nrow(df_mcf_na_plot) > 0 && length(df_mcf_na_plot$Time) == length(df_mcf_na_plot$MCF_NA)) {
    lines(df_mcf_na_plot$Time, df_mcf_na_plot$MCF_NA, type = "s", col = "blue", lwd=1.5) 
    points(df_mcf_na_plot$Time, df_mcf_na_plot$MCF_NA, col = "blue", pch = 19)
  }
  
  legend_text_plp <- paste("PLP Parameters (MLE)\nShape (β):", round(beta_hat_plp_a, 3), 
                       "\nScale (θ):", round(theta_hat_plp_a, 0))
  legend("topleft", legend = legend_text_plp, bty = "n", cex = 0.7, inset = c(0.05, 0.05)) 
  legend("bottomright", legend=c("Parametric PLP", "Non-Parametric N-A"),
         col=c("red", "blue"), lty=1, pch=c(NA,19), lwd=c(2,1.5), cex=0.7)
  
  cat("\n--- Avaliação da Adequação do Modelo PLP (Letra a) ---\n")
  if (beta_hat_plp_a > 1) {
    cat(paste0("O parâmetro de forma do PLP (beta ≈ ", round(beta_hat_plp_a, 3), ") > 1 sugere uma taxa de falha crescente (desgaste).\n"))
    cat("A comparação visual da linha vermelha com os pontos/degraus azuis no gráfico da MCF ajuda a avaliar a adequação.\n")
  } else {
    cat(paste0("O parâmetro de forma do PLP (beta ≈ ", round(beta_hat_plp_a, 3), ") <= 1 sugere taxa de falha constante ou decrescente.\n"))
  }
  
  # Guardar os parâmetros para a letra b)
  beta_estimado_para_b <- beta_hat_plp_a
  theta_estimado_para_b <- theta_hat_plp_a
  
} else {
  cat("Otimização PLP (Simplificada) para Letra A falhou ou não convergiu.\n")
  if (!is.null(mle_fit_plp_simp)) {
    print(mle_fit_plp_simp) 
  }
  beta_estimado_para_b <- NA # Para evitar erro na letra b
  theta_estimado_para_b <- NA
}

#===============================================================================
# LETRA B: Política ótima de manutenção
#===============================================================================
cat("\n\n--- LETRA B: Política Ótima de Manutenção ---\n")

if (!is.na(beta_estimado_para_b) && !is.na(theta_estimado_para_b)) {
  # --- Custos Definidos na Questão b) ---
  CPM <- 1  # Custo da Manutenção Preventiva
  CMR <- 15 # Custo do Reparo Mínimo

  # --- Cálculo do Tempo Ótimo de Manutenção (τ_ótimo) ---
  tau_otimo <- NA
  if (beta_estimado_para_b <= 1) {
    cat("Atenção: O parâmetro beta estimado (", beta_estimado_para_b, ") é menor ou igual a 1.\n")
    cat("A fórmula para o tempo ótimo de PM pode não ser aplicável ou a PM pode não ser benéfica.\n")
  } else {
    termo_custo <- CPM / ((beta_estimado_para_b - 1) * CMR)
    if (termo_custo < 0) {
        cat("Atenção: Termo de custo [CPM / ((β-1) * CMR)] é negativo. Verifique os custos e beta.\n")
    } else {
        tau_otimo <- theta_estimado_para_b * (termo_custo)^(1 / beta_estimado_para_b)
    }
  }

  if (!is.na(tau_otimo)) {
    cat("Parâmetros do PLP Usados (da Letra a):\n")
    cat("  Beta (Forma) Estimado:", format(beta_estimado_para_b, digits = 5), "\n")
    cat("  Theta (Escala) Estimado:", format(theta_estimado_para_b, digits = 5), "\n")
    cat("Custos:\n")
    cat("  Custo de PM (CPM):", CPM, "\n")
    cat("  Custo de MR (CMR):", CMR, "\n")
    cat("Tempo Ótimo de Manutenção Preventiva (τ_ótimo):", format(tau_otimo, digits = 2), "horas\n")
    cat("Tempo Ótimo de Manutenção Preventiva (τ_ótimo):", format(tau_otimo / 24, digits = 2), "dias\n\n")
    
    # --- Plotar a Função de Custo H(τ) para visualizar o mínimo ---
    plot_custo_H_tau_focado <- function(beta, theta, cpm_val, cmr_val, tau_otimo_calc) {
      range_offset <- tau_otimo_calc * 0.10 # Foco de 10% em torno do ótimo
      if (range_offset < 200) range_offset = 200 # Mínimo range para visualização
      
      tau_min_plot <- max(100, tau_otimo_calc - range_offset) 
      tau_max_plot <- tau_otimo_calc + range_offset
      
      taus <- seq(tau_min_plot, tau_max_plot, length.out = 200) 
      H_tau_values <- numeric(length(taus))
      
      for (i in 1:length(taus)) {
        tau_i <- taus[i]
        m_tau_i <- (tau_i / theta)^beta
        H_tau_values[i] <- (1 / tau_i) * (cpm_val + cmr_val * m_tau_i)
      }
      
      valid_indices <- is.finite(H_tau_values)
      min_H_val <- min(H_tau_values[valid_indices], na.rm = TRUE)
      max_H_val <- max(H_tau_values[valid_indices], na.rm = TRUE)
      y_range_padding <- (max_H_val - min_H_val) * 0.1
      if (y_range_padding < 1e-6) y_range_padding = (min_H_val)*0.05 # Mínimo padding
      if (y_range_padding == 0 && min_H_val == 0) y_range_padding = 0.00001


      ylim_plot_focused <- c(max(0, min_H_val - y_range_padding), max_H_val + y_range_padding)
      if(ylim_plot_focused[1] >= ylim_plot_focused[2]) ylim_plot_focused[1] = ylim_plot_focused[2]*0.9


      plot(taus[valid_indices], H_tau_values[valid_indices], type = "l", col = "blue", lwd = 2,
           xlab = "Intervalo de Manutenção Preventiva τ (horas)",
           ylab = "Custo Esperado por Unidade de Tempo H(τ)",
           main = "Função de Custo da Manutenção Preventiva (Focada no Ótimo)",
           ylim = ylim_plot_focused) 
      grid()
      
      if (!is.na(tau_otimo_calc) && tau_otimo_calc > 0 && is.finite(tau_otimo_calc)) {
        m_tau_opt <- (tau_otimo_calc / theta)^beta
        H_tau_opt <- (1 / tau_otimo_calc) * (cpm_val + cmr_val * m_tau_opt)
        if(is.finite(H_tau_opt)){
            abline(v = tau_otimo_calc, col = "red", lty = 2, lwd = 1.5)
            points(tau_otimo_calc, H_tau_opt, col = "red", pch = 19, cex = 1.2)
            text_label <- paste0("τ_ótimo ≈ ", round(tau_otimo_calc,0), "\nH(τ_ótimo) ≈ ", format(H_tau_opt, scientific = FALSE, digits = 5))
            # Ajustar posição do texto para não sobrepor o ponto
            text_x_pos <- if (tau_otimo_calc < mean(range(taus[valid_indices]))) tau_otimo_calc + (0.05 * diff(range(taus[valid_indices]))) else tau_otimo_calc - (0.05 * diff(range(taus[valid_indices])))
            text_y_pos <- H_tau_opt + (0.05 * diff(ylim_plot_focused))
            text_adj <- if (tau_otimo_calc < mean(range(taus[valid_indices]))) c(0,0) else c(1,0) # left or right align

            text(text_x_pos, text_y_pos, labels = text_label, cex = 0.8, col = "red", adj = text_adj) 
        }
      }
    } # Fim da função plot_custo_H_tau_focado
    
    # Plotar a função de custo focada
    plot_custo_H_tau_focado(beta = beta_estimado_para_b, 
                            theta = theta_estimado_para_b, 
                            cpm_val = CPM, 
                            cmr_val = CMR, 
                            tau_otimo_calc = tau_otimo)
                            
  } else {
    cat("Não foi possível calcular o tempo ótimo de manutenção para a letra b, pois os parâmetros da letra a não foram estimados corretamente ou beta <= 1.\n")
  }
} else {
    cat("Parâmetros da letra (a) não foram estimados com sucesso. Não é possível prosseguir para a letra (b).\n")
}