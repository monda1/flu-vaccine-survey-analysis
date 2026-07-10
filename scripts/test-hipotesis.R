library(dplyr)
library(ggplot2)

# =====================================================================
# TEMA 10-11: Contraste de proporciones (health_worker vs resto)
# e Intervalos de confianza por grupo de edad
# =====================================================================

# --------------------------------------------------------------
# 0. Comprobación de las condiciones de aplicación del test
# --------------------------------------------------------------
# El contraste chi-cuadrado (y su equivalente en z) es una aproximación
# asintótica: solo es fiable si, en cada grupo, se espera al menos
# ~5 éxitos y ~5 fracasos (n*p >= 5 y n*(1-p) >= 5). Si esto no se
# cumple en algún grupo, el p-valor del chi-cuadrado no es de fiar
# y conviene usar el test exacto de Fisher en su lugar.

comprobar_condiciones <- function(datos, var_vacuna, nombre_vacuna) {
  df <- datos %>% filter(!is.na(.data[[var_vacuna]]), !is.na(health_worker))
  tabla <- table(df$health_worker, df[[var_vacuna]])
  
  n1 <- sum(tabla["1", ]); n0 <- sum(tabla["0", ])
  p1 <- tabla["1", "1"] / n1; p0 <- tabla["0", "1"] / n0
  
  cond <- c(n1 * p1, n1 * (1 - p1), n0 * p0, n0 * (1 - p0))
  ok <- all(cond >= 5)
  
  cat("--- Condiciones de aplicación (", nombre_vacuna, ") ---\n", sep = "")
  cat("n*p y n*(1-p) por grupo:", round(cond, 1), "\n")
  cat(ifelse(ok,
             "-> Se cumplen (n*p, n*(1-p) >= 5): la aproximación normal/chi-cuadrado es válida.\n\n",
             "-> NO se cumplen todas: usar test exacto de Fisher en vez del chi-cuadrado.\n\n"))
  return(ok)
}

# --------------------------------------------------------------
# 1. Contraste de hipótesis: proporción de vacunados
#    en sanitarios vs no sanitarios
# --------------------------------------------------------------
# H0: p_sanitarios = p_no_sanitarios  (no hay diferencia)
# H1: p_sanitarios != p_no_sanitarios (dos colas)
#
# Aquí tiene sentido plantear
# un contraste UNILATERAL (p_sanitarios > p_no_sanitarios), porque
# hay una razón teórica previa a los datos (mayor exposición al virus,
# acceso más fácil a la vacuna en el propio centro de trabajo, y
# probablemente mayor tasa de recomendación médica -> doctor_recc).
# Usar un test a una cola cuando hay justificación teórica previa
# aumenta la potencia del contraste. Se calculan ambas versiones
# para comparar.

contraste_proporciones <- function(datos, var_vacuna, nombre_vacuna) {
  df <- datos %>%
    filter(!is.na(.data[[var_vacuna]]), !is.na(health_worker))
  
  tabla <- table(df$health_worker, df[[var_vacuna]])
  x_vacunados <- c(tabla["1", "1"], tabla["0", "1"])
  n_total <- c(sum(tabla["1", ]), sum(tabla["0", ]))
  
  cumple <- comprobar_condiciones(datos, var_vacuna, nombre_vacuna)
  
  test_2c <- prop.test(x = x_vacunados, n = n_total, alternative = "two.sided", correct = TRUE)
  test_1c <- prop.test(x = x_vacunados, n = n_total, alternative = "greater", correct = TRUE)
  
  # Alternativa exacta (válida incluso si fallan las condiciones de arriba)
  test_fisher <- fisher.test(matrix(c(x_vacunados[1], n_total[1] - x_vacunados[1],
                                      x_vacunados[2], n_total[2] - x_vacunados[2]),
                                    nrow = 2, byrow = TRUE))
  
  p1 <- x_vacunados[1] / n_total[1]
  p2 <- x_vacunados[2] / n_total[2]
  
  # Tamaño del efecto: Cohen's h para proporciones.
  # A diferencia del p-valor (que depende del tamaño muestral y solo
  # dice "hay diferencia sí/no"), h cuantifica CUÁNTO difieren las
  # proporciones en una escala comparable entre estudios.
  # Referencia : |h| ~ 0.2 pequeño, 0.5 mediano, 0.8 grande.
  h_cohen <- 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))
  
  # Odds ratio: medida de asociación habitual en estudios epidemiológicos,
  # complementaria a la diferencia de proporciones.
  or_val <- (x_vacunados[1] * (n_total[2] - x_vacunados[2])) /
    (x_vacunados[2] * (n_total[1] - x_vacunados[1]))
  
  cat("=== Contraste de proporciones:", nombre_vacuna, "===\n")
  cat("Proporción sanitarios:", round(p1, 4), "(n =", n_total[1], ")\n")
  cat("Proporción no sanitarios:", round(p2, 4), "(n =", n_total[2], ")\n")
  cat("Diferencia observada:", round(p1 - p2, 4), "\n")
  cat("Tamaño del efecto (Cohen's h):", round(h_cohen, 4), "\n")
  cat("Odds ratio:", round(or_val, 4), "\n")
  cat("--- Test dos colas (chi-cuadrado) ---\n")
  cat("Estadístico:", round(test_2c$statistic, 4),
      "| p-valor:", format.pval(test_2c$p.value, digits = 4), "\n")
  cat("IC 95% diferencia de proporciones:",
      round(test_2c$conf.int[1], 4), "a", round(test_2c$conf.int[2], 4), "\n")
  cat("--- Test una cola (H1: p_sanitarios > p_no_sanitarios) ---\n")
  cat("p-valor:", format.pval(test_1c$p.value, digits = 4), "\n")
  cat("--- Test exacto de Fisher (alternativa sin aproximación asintótica) ---\n")
  cat("Odds ratio (Fisher):", round(test_fisher$estimate, 4),
      "| p-valor:", format.pval(test_fisher$p.value, digits = 4), "\n\n")
  
  return(list(chi2 = test_2c, fisher = test_fisher, h_cohen = h_cohen, or = or_val))
}

test_h1n1 <- contraste_proporciones(datos, "h1n1_vaccine", "H1N1")
test_seasonal <- contraste_proporciones(datos, "seasonal_vaccine", "Estacional")

# Se han hecho dos contrastes (H1N1 y estacional) sobre
# la misma comparación (sanitarios vs resto). Si se van a interpretar
# conjuntamente como "confirmación cruzada" de la misma hipótesis,
# el umbral de significación efectivo para no inflar el error de tipo I
# sería alpha/2 = 0.025 (corrección de Bonferroni), no 0.05.

alpha_bonferroni <- 0.05 / 2
cat("Umbral de significación ajustado (Bonferroni, 2 tests):", alpha_bonferroni, "\n\n")

# --------------------------------------------------------------
# 2. Verificación del estadístico
# --------------------------------------------------------------
verificacion_manual <- function(datos, var_vacuna) {
  df <- datos %>% filter(!is.na(.data[[var_vacuna]]), !is.na(health_worker))
  
  p1 <- mean(df[[var_vacuna]][df$health_worker == 1])
  p2 <- mean(df[[var_vacuna]][df$health_worker == 0])
  n1 <- sum(df$health_worker == 1)
  n2 <- sum(df$health_worker == 0)
  
  p_pool <- (sum(df[[var_vacuna]][df$health_worker == 1]) +
               sum(df[[var_vacuna]][df$health_worker == 0])) / (n1 + n2)
  
  ee_pool <- sqrt(p_pool * (1 - p_pool) * (1/n1 + 1/n2))
  z <- (p1 - p2) / ee_pool
  p_valor <- 2 * (1 - pnorm(abs(z)))
  
  cat("Verificación manual -> z =", round(z, 4),
      "| p-valor =", format.pval(p_valor, digits = 4), "\n\n")
}

verificacion_manual(datos, "h1n1_vaccine")
verificacion_manual(datos, "seasonal_vaccine")

# --------------------------------------------------------------
# 3. Intervalos de confianza para la tasa de vacunación
#    desagregada por grupo de edad
# --------------------------------------------------------------
# Se calculan dps versiones del IC:
#  - Wald: la fórmula clásica p_hat +- z*EE (la que ya tenías).
#  - Wilson (score): más precisa cuando p_hat está cerca de 0 o 1,
#    o cuando n no es muy grande, y nunca se sale de [0,1] sin
#    necesidad de recortar manualmente como hace Wald.


ic_por_grupo <- function(datos, var_vacuna, var_grupo, nombre_vacuna) {
  resumen <- datos %>%
    filter(!is.na(.data[[var_vacuna]]), !is.na(.data[[var_grupo]])) %>%
    group_by(grupo = .data[[var_grupo]]) %>%
    summarise(
      n = n(),
      vacunados = sum(.data[[var_vacuna]]),
      p_hat = vacunados / n,
      ee = sqrt(p_hat * (1 - p_hat) / n),
      ic_inf_wald = p_hat - qnorm(0.975) * ee,
      ic_sup_wald = p_hat + qnorm(0.975) * ee,
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      # Intervalo de Wilson calculado manualmente con su fórmula cerrada
      z = qnorm(0.975),
      denom = 1 + z^2 / n,
      centro = (p_hat + z^2 / (2 * n)) / denom,
      margen = (z / denom) * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2)),
      ic_inf_wilson = centro - margen,
      ic_sup_wilson = centro + margen
    ) %>%
    ungroup() %>%
    mutate(
      ic_inf_wald = pmax(ic_inf_wald, 0),
      ic_sup_wald = pmin(ic_sup_wald, 1),
      amplitud_wald = ic_sup_wald - ic_inf_wald,
      amplitud_wilson = ic_sup_wilson - ic_inf_wilson,
      vacuna = nombre_vacuna
    ) %>%
    select(grupo, n, vacunados, p_hat,
           ic_inf_wald, ic_sup_wald, amplitud_wald,
           ic_inf_wilson, ic_sup_wilson, amplitud_wilson, vacuna)
  
  resumen
}

ic_h1n1_edad <- ic_por_grupo(datos, "h1n1_vaccine", "age_group", "H1N1")
ic_seasonal_edad <- ic_por_grupo(datos, "seasonal_vaccine", "age_group", "Estacional")

cat("=== IC 95% tasa de vacunación H1N1 por grupo de edad (Wald vs Wilson) ===\n")
print(ic_h1n1_edad)
cat("\n=== IC 95% tasa de vacunación Estacional por grupo de edad (Wald vs Wilson) ===\n")
print(ic_seasonal_edad)

# Se calculan 5-6 IC por vacuna (uno por grupo de
# edad). Cada IC individual tiene un 95% de confianza, pero la
# probabilidad de que todos acierten a la vez es menor. Si el objetivo
# es comparar grupos de edad entre sí (no solo describir cada uno por
# separado), habría que ajustar el nivel de confianza con Bonferroni:
# 1 - 0.05/k, donde k es el número de grupos comparados.

k_grupos <- n_distinct(ic_h1n1_edad$grupo)
conf_ajustada <- 1 - 0.05 / k_grupos
cat("\nNivel de confianza ajustado (Bonferroni,", k_grupos, "grupos):",
    round(conf_ajustada, 4), "\n\n")

# --------------------------------------------------------------
# 4. Intervalos por edad
#    (se usa Wilson por ser más fiable en los extremos)
# --------------------------------------------------------------
ic_combinado <- bind_rows(ic_h1n1_edad, ic_seasonal_edad)

ggplot(ic_combinado, aes(x = grupo, y = p_hat, color = vacuna)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = ic_inf_wilson, ymax = ic_sup_wilson),
                position = position_dodge(width = 0.5), width = 0.2) +
  labs(title = "Tasa de vacunación por grupo de edad (IC 95% de Wilson)",
       subtitle = "Barras = intervalo de Wilson; ver script para comparación con Wald",
       x = "Grupo de edad", y = "Proporción vacunada", color = "Vacuna") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# --------------------------------------------------------------
# 5. Diferencia sanitarios vs resto,
#    con el tamaño del efecto (Cohen's h) señalado en el subtítulo
# --------------------------------------------------------------
comparacion_hw <- datos %>%
  filter(!is.na(health_worker)) %>%
  mutate(grupo_hw = ifelse(health_worker == 1, "Sanitario", "No sanitario")) %>%
  group_by(grupo_hw) %>%
  summarise(
    p_h1n1 = mean(h1n1_vaccine, na.rm = TRUE),
    p_seasonal = mean(seasonal_vaccine, na.rm = TRUE),
    .groups = "drop"
  )

print(comparacion_hw)

cat("\nTamaño del efecto (Cohen's h):\n",
    " H1N1:", round(test_h1n1$h_cohen, 3),
    "| Estacional:", round(test_seasonal$h_cohen, 3), "\n")