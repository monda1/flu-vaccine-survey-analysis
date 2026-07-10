# ============================================================
# Fenómenos aleatorios, sucesos, axiomas, probabilidad
#          condicionada, independencia, teorema de la probabilidad
#          total y teorema de Bayes
#
# Integra el análisis completo (espacio muestral,
# axiomas, probabilidad condicionada, independencia, teorema de la
# probabilidad total y teorema de Bayes) junto con la interpretación
# de resultados.
#
# Requiere que el objeto `datos` ya exista en el entorno (generado
# por scripts/01_preprocesamiento.R) 
# ============================================================

library(dplyr)

# Paleta reutilizada del script de preprocesamiento (por si se ejecuta suelto)
if (!exists("pal_azul")) pal_azul <- "#3B82F6"
if (!exists("pal_rojo")) pal_rojo <- "#EF4444"

# Aseguramos formato numérico 0/1 de las variables relevantes
datos$h1n1_vaccine     <- as.numeric(as.character(datos$h1n1_vaccine))
datos$seasonal_vaccine <- as.numeric(as.character(datos$seasonal_vaccine))
datos$doctor_recc_h1n1 <- as.numeric(as.character(datos$doctor_recc_h1n1))
datos$health_worker    <- as.numeric(as.character(datos$health_worker))

# ------------------------------------------------------------
# 0. FENÓMENO ALEATORIO 
# ------------------------------------------------------------
# Cada fila de `datos` es el resultado de un fenómeno aleatorio: la
# respuesta de un individuo distinto de la población a la encuesta NHFS
# 2009. No se puede predecir con certeza si una persona concreta se
# vacunará, pero sí se puede estudiar la regularidad estadística del
# fenómeno repitiéndolo sobre muchos individuos (los encuestados), lo
# que permite estimar probabilidades

cat("=== FENÓMENO ALEATORIO ===\n")
cat("Nº de realizaciones del fenómeno (encuestados con dato válido):",
    sum(!is.na(datos$h1n1_vaccine)), "\n")
cat("Resultado posible de cada realización: {No vacunado (0), Vacunado (1)}\n\n")

# ------------------------------------------------------------
# 1. ESPACIO MUESTRAL Y SUCESOS
# ------------------------------------------------------------
# Ω: conjunto de todos los respondentes de la encuesta NHFS 2009
# Suceso A: "el encuestado recibió la vacuna H1N1"       -> h1n1_vaccine == 1
# Suceso B: "el encuestado recibió la vacuna estacional" -> seasonal_vaccine == 1

datos_venn <- datos %>%
  filter(!is.na(h1n1_vaccine), !is.na(seasonal_vaccine))

n_omega   <- nrow(datos_venn)
n_A       <- sum(datos_venn$h1n1_vaccine == 1)
n_B       <- sum(datos_venn$seasonal_vaccine == 1)
n_AiB     <- sum(datos_venn$h1n1_vaccine == 1 & datos_venn$seasonal_vaccine == 1)
n_AuB     <- sum(datos_venn$h1n1_vaccine == 1 | datos_venn$seasonal_vaccine == 1)
n_ninguna <- n_omega - n_AuB

cat("=== ESPACIO MUESTRAL Y SUCESOS ===\n")
cat("Ω  (total de encuestados)                  :", n_omega, "\n")
cat("A  (vacunados H1N1)                        :", n_A, "\n")
cat("B  (vacunados estacional)                  :", n_B, "\n")
cat("A ∩ B (ambas vacunas)                      :", n_AiB, "\n")
cat("A ∪ B (al menos una vacuna)                :", n_AuB, "\n")
cat("Ni A ni B (ninguna vacuna)                 :", n_ninguna, "\n\n")

# Diagrama de Venn
dibujar_circulo <- function(cx, cy, r, ...) {
  theta <- seq(0, 2 * pi, length.out = 200)
  polygon(cx + r * cos(theta), cy + r * sin(theta), ...)
}

plot(NA, xlim = c(-1.3, 2.3), ylim = c(-1.6, 1.5), asp = 1,
     axes = FALSE, xlab = "", ylab = "",
     main = "Diagrama de Venn: vacunación H1N1 (A) y estacional (B)")
dibujar_circulo(0, 0, 1, col = adjustcolor(pal_azul, alpha.f = 0.4), border = pal_azul, lwd = 2)
dibujar_circulo(1, 0, 1, col = adjustcolor(pal_rojo, alpha.f = 0.4), border = pal_rojo, lwd = 2)
text(-0.7, 0.9, "A: H1N1",       col = pal_azul, font = 2)
text(1.7,  0.9, "B: Estacional", col = pal_rojo, font = 2)
text(-0.35, 0, n_A - n_AiB)
text(1.35,  0, n_B - n_AiB)
text(0.5,   0, n_AiB)
text(0.5, -1.4, paste0("Ninguna vacuna (fuera de A y B): ", n_ninguna))

#  Sucesos incompatibles (mutuamente excluyentes)
# A partir de A y B construimos una partición de Ω en 4 sucesos disjuntos
# dos a dos (su intersección es siempre el conjunto vacío)
n_solo_A <- sum(datos_venn$h1n1_vaccine == 1 & datos_venn$seasonal_vaccine == 0)
n_solo_B <- sum(datos_venn$h1n1_vaccine == 0 & datos_venn$seasonal_vaccine == 1)

suma_disjuntos <- n_solo_A + n_solo_B + n_AiB + n_ninguna

cat("=== SUCESOS INCOMPATIBLES (partición de Ω) ===\n")
cat("Solo H1N1 (A ∩ B^c)                        :", n_solo_A, "\n")
cat("Solo estacional (A^c ∩ B)                  :", n_solo_B, "\n")
cat("Ambas vacunas (A ∩ B)                      :", n_AiB, "\n")
cat("Ninguna vacuna (A^c ∩ B^c)                 :", n_ninguna, "\n")
cat("Suma de los 4 sucesos disjuntos            :", suma_disjuntos,
    " (debe coincidir con Ω =", n_omega, ")\n")
cat("Estos 4 sucesos son incompatibles entre sí: la intersección de\n")
cat("cualquier par de ellos es el conjunto vacío (probabilidad 0).\n\n")

# ------------------------------------------------------------
# 2. AXIOMAS Y PROPIEDADES DE LA PROBABILIDAD 
# ------------------------------------------------------------
p_A   <- n_A   / n_omega
p_B   <- n_B   / n_omega
p_AiB <- n_AiB / n_omega

# Regla de la unión: P(A ∪ B) = P(A) + P(B) - P(A ∩ B)
p_AuB_teorico  <- p_A + p_B - p_AiB
p_AuB_empirico <- n_AuB / n_omega

# Regla del complementario: P(A^c) = 1 - P(A)
p_Ac         <- 1 - p_A
p_Ac_directo <- sum(datos_venn$h1n1_vaccine == 0) / n_omega

cat("=== AXIOMAS Y PROPIEDADES DE LA PROBABILIDAD ===\n")
cat("P(A)                                              :", round(p_A, 4), "\n")
cat("P(B)                                              :", round(p_B, 4), "\n")
cat("P(A ∩ B)                                          :", round(p_AiB, 4), "\n")
cat("P(A ∪ B) = P(A)+P(B)-P(A∩B)   [regla de la unión] :", round(p_AuB_teorico, 4), "\n")
cat("P(A ∪ B)   [conteo directo]                       :", round(p_AuB_empirico, 4), "\n")
cat("Diferencia (debe ser ≈ 0)                         :", round(p_AuB_teorico - p_AuB_empirico, 6), "\n\n")

cat("P(A^c) = 1 - P(A)   [complementario]              :", round(p_Ac, 4), "\n")
cat("P(A^c)   [conteo directo]                         :", round(p_Ac_directo, 4), "\n\n")

# --- Monotonía: como A∩B ⊆ A ⊆ A∪B, se cumple P(A∩B) ≤ P(A) ≤ P(A∪B) ---
cat("=== MONOTONÍA DE LA PROBABILIDAD ===\n")
cat("P(A ∩ B)                                          :", round(p_AiB, 4), "\n")
cat("P(A)                                              :", round(p_A, 4), "\n")
cat("P(A ∪ B)                                          :", round(p_AuB_empirico, 4), "\n")
cat("¿Se cumple P(A∩B) ≤ P(A) ≤ P(A∪B)?                :",
    p_AiB <= p_A && p_A <= p_AuB_empirico, "\n")
cat("(A∩B es subconjunto de A, y A es subconjunto de A∪B, por eso sus\n")
cat("probabilidades quedan necesariamente ordenadas así.)\n\n")

# --- Leyes de De Morgan: (A∪B)^c = A^c∩B^c   y   (A∩B)^c = A^c∪B^c ---
n_AuB_c   <- n_omega - n_AuB                                              # (A ∪ B)^c
n_Ac_i_Bc <- sum(datos_venn$h1n1_vaccine == 0 & datos_venn$seasonal_vaccine == 0)  # A^c ∩ B^c
n_AiB_c   <- n_omega - n_AiB                                              # (A ∩ B)^c
n_Ac_u_Bc <- sum(datos_venn$h1n1_vaccine == 0 | datos_venn$seasonal_vaccine == 0)  # A^c ∪ B^c

cat("=== LEYES DE DE MORGAN ===\n")
cat("(A ∪ B)^c                                          :", n_AuB_c, "\n")
cat("A^c ∩ B^c                                          :", n_Ac_i_Bc,
    " (debe coincidir con (A ∪ B)^c)\n")
cat("(A ∩ B)^c                                          :", n_AiB_c, "\n")
cat("A^c ∪ B^c                                          :", n_Ac_u_Bc,
    " (debe coincidir con (A ∩ B)^c)\n\n")

# ------------------------------------------------------------
# 3. PROBABILIDAD CONDICIONADA
# ------------------------------------------------------------
# Suceso C: "el encuestado es trabajador sanitario" -> health_worker == 1
# ¿Ser trabajador sanitario cambia la probabilidad de vacunarse contra H1N1?

datos_cond <- datos %>%
  filter(!is.na(h1n1_vaccine), !is.na(health_worker))

p_C   <- mean(datos_cond$health_worker)
p_AiC <- mean(datos_cond$h1n1_vaccine == 1 & datos_cond$health_worker == 1)

# Definición: P(A|C) = P(A ∩ C) / P(C)
p_A_dado_C <- p_AiC / p_C

p_A_dado_C_directo <- datos_cond %>%
  filter(health_worker == 1) %>%
  summarise(p = mean(h1n1_vaccine)) %>%
  pull(p)

cat("=== PROBABILIDAD CONDICIONADA: P(A|C)  ===\n")
cat("C = ser trabajador sanitario\n")
cat("P(C)                                              :", round(p_C, 4), "\n")
cat("P(A ∩ C)                                          :", round(p_AiC, 4), "\n")
cat("P(A|C) = P(A∩C)/P(C)   [definición]               :", round(p_A_dado_C, 4), "\n")
cat("P(A|C)   [cálculo directo]                        :", round(p_A_dado_C_directo, 4), "\n")
cat("P(A)  marginal, sin condicionar (referencia)      :", round(p_A, 4), "\n")
cat("Interpretación: trabajar en sanidad",
    ifelse(p_A_dado_C > p_A, "aumenta", "disminuye"),
    "la probabilidad de vacunarse frente a H1N1.\n\n")

# ------------------------------------------------------------
# 4. INDEPENDENCIA DE SUCESOS
# ------------------------------------------------------------
vars_behavioral <- c("behavioral_antiviral_meds", "behavioral_avoidance",
                     "behavioral_face_mask", "behavioral_wash_hands",
                     "behavioral_large_gatherings", "behavioral_outside_home",
                     "behavioral_touch_face")

datos[vars_behavioral] <- lapply(datos[vars_behavioral], function(x) as.numeric(as.character(x)))

matriz_cor_behavioral <- cor(datos[vars_behavioral], use = "pairwise.complete.obs")
cat("=== Matriz de correlación entre comportamientos preventivos ===\n")
print(round(matriz_cor_behavioral, 3))
cat("\n")

combinaciones <- combn(vars_behavioral, 2, simplify = FALSE)

# 4a. Test de hipótesis: Chi-cuadrado de independencia
resultados_chi <- data.frame(
  var1 = character(), var2 = character(),
  chi2 = numeric(), p_valor = numeric(), independientes = character(),
  stringsAsFactors = FALSE
)

for (par in combinaciones) {
  tabla <- table(datos[[par[1]]], datos[[par[2]]])
  test <- suppressWarnings(chisq.test(tabla))
  resultados_chi <- rbind(resultados_chi, data.frame(
    var1 = par[1], var2 = par[2],
    chi2 = round(test$statistic, 2),
    p_valor = round(test$p.value, 5),
    independientes = ifelse(test$p.value > 0.05, "Sí (no se rechaza H0)", "No (se rechaza H0)")
  ))
}

cat("=== Tests de independencia (Chi-cuadrado) por pares ===\n")
print(resultados_chi, row.names = FALSE)
cat("\nInterpretación: si p_valor < 0.05, se rechaza la hipótesis de independencia,\n")
cat("es decir, existe evidencia de asociación entre ese par de comportamientos.\n\n")

# 4b. Definición formal de independencia: P(A ∩ B) = P(A) · P(B)
# Complementa el test de hipótesis con la definición del Tema 1.6
resultados_indep_formal <- data.frame(
  var1 = character(), var2 = character(),
  p_interseccion = numeric(), p_producto = numeric(), diferencia = numeric(),
  stringsAsFactors = FALSE
)

for (par in combinaciones) {
  p1    <- mean(datos[[par[1]]], na.rm = TRUE)
  p2    <- mean(datos[[par[2]]], na.rm = TRUE)
  p_int <- mean(datos[[par[1]]] == 1 & datos[[par[2]]] == 1, na.rm = TRUE)
  resultados_indep_formal <- rbind(resultados_indep_formal, data.frame(
    var1 = par[1], var2 = par[2],
    p_interseccion = round(p_int, 4),
    p_producto      = round(p1 * p2, 4),
    diferencia      = round(p_int - p1 * p2, 4)
  ))
}

cat("=== Independencia : P(A∩B) vs P(A)·P(B) ===\n")
print(resultados_indep_formal, row.names = FALSE)
cat("\nSi P(A∩B) ≈ P(A)·P(B) (diferencia cercana a 0), los sucesos son independientes\n")


# 4c. Independencia mutua de tres sucesos (no solo dos a dos)
# Para que tres sucesos D, E, F sean mutuamente independientes no basta con
# que sean independientes dos a dos: además debe cumplirse
# P(D ∩ E ∩ F) = P(D)·P(E)·P(F).
vars_tres <- c("behavioral_wash_hands", "behavioral_touch_face", "behavioral_face_mask")

p_d <- mean(datos[[vars_tres[1]]], na.rm = TRUE)
p_e <- mean(datos[[vars_tres[2]]], na.rm = TRUE)
p_f <- mean(datos[[vars_tres[3]]], na.rm = TRUE)

p_def <- mean(datos[[vars_tres[1]]] == 1 & datos[[vars_tres[2]]] == 1 &
                datos[[vars_tres[3]]] == 1, na.rm = TRUE)

p_producto_tres <- p_d * p_e * p_f

cat("=== INDEPENDENCIA MUTUA DE TRES SUCESOS ===\n")
cat("D =", vars_tres[1], ", E =", vars_tres[2], ", F =", vars_tres[3], "\n")
cat("P(D ∩ E ∩ F)                                      :", round(p_def, 4), "\n")
cat("P(D)·P(E)·P(F)                                    :", round(p_producto_tres, 4), "\n")
cat("Diferencia                                        :", round(p_def - p_producto_tres, 4), "\n")
cat("Los pares de este trío ya se evaluaron dos a dos en la tabla Chi-\n")
cat("cuadrado anterior; para independencia MUTUA hace falta además que se\n")
cat("cumpla esta igualdad conjunta, no solo la de los pares por separado.\n\n")

# ------------------------------------------------------------
# 5. TEOREMA DE LA PROBABILIDAD TOTAL
# ------------------------------------------------------------
# Partición de Ω por grupo de edad: los B_i = {age_group = i} son disjuntos
# y cubren toda la muestra -> partición válida para aplicar el teorema.

datos_particion <- datos %>%
  filter(!is.na(h1n1_vaccine), !is.na(age_group))

tabla_particion <- datos_particion %>%
  group_by(age_group) %>%
  summarise(
    n_grupo     = n(),
    p_Bi        = n_grupo / nrow(datos_particion),
    p_A_dado_Bi = mean(h1n1_vaccine, na.rm = TRUE),
    aporte      = p_Bi * p_A_dado_Bi,
    .groups = "drop"
  )

p_A_total   <- sum(tabla_particion$aporte)
p_A_directo <- mean(datos_particion$h1n1_vaccine)

cat("=== TEOREMA DE LA PROBABILIDAD TOTAL ===\n")
cat("Partición de Ω por age_group (grupos disjuntos que cubren toda la muestra)\n\n")
print(tabla_particion, n = Inf)
cat("\nP(A) = Σ P(A|B_i)·P(B_i)   [prob. total]          :", round(p_A_total, 4), "\n")
cat("P(A)   [cálculo directo, sin partición]           :", round(p_A_directo, 4), "\n")
cat("Diferencia (debe ser ≈ 0)                         :", round(p_A_total - p_A_directo, 6), "\n\n")

# ------------------------------------------------------------
# 6. TEOREMA DE BAYES
# ------------------------------------------------------------
# El denominador P(Recc=1) se calcula aquí mediante el teorema de la
# probabilidad tota, usando la partición {Vacuna=1, Vacuna=0},
# en vez de tomarlo directamente de la media así se conectan ambos
# teoremas 

datos_bayes <- datos %>%
  filter(!is.na(h1n1_vaccine), !is.na(doctor_recc_h1n1))

p_vacuna    <- mean(datos_bayes$h1n1_vaccine)
p_no_vacuna <- 1 - p_vacuna

p_recc_dado_vacuna <- datos_bayes %>%
  filter(h1n1_vaccine == 1) %>%
  summarise(p = mean(doctor_recc_h1n1)) %>%
  pull(p)

p_recc_dado_no_vacuna <- datos_bayes %>%
  filter(h1n1_vaccine == 0) %>%
  summarise(p = mean(doctor_recc_h1n1)) %>%
  pull(p)

# Teorema de la probabilidad total aplicado al denominador de Bayes:
# P(Recc=1) = P(Recc=1|Vacuna=1)·P(Vacuna=1) + P(Recc=1|Vacuna=0)·P(Vacuna=0)
p_recc_total_prob <- p_recc_dado_vacuna * p_vacuna + p_recc_dado_no_vacuna * p_no_vacuna
p_recc_directo     <- mean(datos_bayes$doctor_recc_h1n1)

cat("=== DENOMINADOR DE BAYES VÍA TEOREMA DE LA PROBABILIDAD TOTAL ===\n")
cat("P(Recc=1) [teorema de la probabilidad total]      :", round(p_recc_total_prob, 4), "\n")
cat("P(Recc=1) [cálculo directo]                       :", round(p_recc_directo, 4), "\n")
cat("Diferencia (debe ser ≈ 0)                         :", round(p_recc_total_prob - p_recc_directo, 6), "\n\n")

# Teorema de Bayes:
# P(Vacuna=1 | Recc=1) = P(Recc=1 | Vacuna=1) * P(Vacuna=1) / P(Recc=1)
p_vacuna_dado_recc_bayes <- (p_recc_dado_vacuna * p_vacuna) / p_recc_total_prob

p_vacuna_dado_recc_directo <- datos_bayes %>%
  filter(doctor_recc_h1n1 == 1) %>%
  summarise(p = mean(h1n1_vaccine)) %>%
  pull(p)

# P(Vacuna=1 | Recc=0) -> para completar el panorama
p_vacuna_dado_no_recc <- datos_bayes %>%
  filter(doctor_recc_h1n1 == 0) %>%
  summarise(p = mean(h1n1_vaccine)) %>%
  pull(p)

cat("=== TEOREMA DE BAYES: Recomendación médica y vacunación H1N1 ===\n")
cat("P(Vacuna = 1)                         [marginal / a priori] :", round(p_vacuna, 4), "\n")
cat("P(Recomendación = 1)                  [prob. total]         :", round(p_recc_total_prob, 4), "\n")
cat("P(Recomendación = 1 | Vacuna = 1)     [verosimilitud]       :", round(p_recc_dado_vacuna, 4), "\n")
cat("P(Vacuna = 1 | Recomendación = 1)     [Bayes]               :", round(p_vacuna_dado_recc_bayes, 4), "\n")
cat("P(Vacuna = 1 | Recomendación = 1)     [cálculo directo]     :", round(p_vacuna_dado_recc_directo, 4), "\n")
cat("P(Vacuna = 1 | Recomendación = 0)     [contraste]           :", round(p_vacuna_dado_no_recc, 4), "\n")
cat("Incremento absoluto de probabilidad (posterior - prior)     :",
    round(p_vacuna_dado_recc_bayes - p_vacuna, 4), "\n")
cat("Razón de probabilidades (posterior / prior)                 :",
    round(p_vacuna_dado_recc_bayes / p_vacuna, 4), "\n\n")

# --- Bayes con más de dos hipótesis: partición por age_group ---
# Reutilizamos la partición de la sección 5 (B_i = age_group) y calculamos,
# para cada grupo de edad, la probabilidad a posteriori de pertenecer a
# ese grupo sabiendo que el encuestado se vacunó: P(B_i | A).
tabla_bayes_edad <- tabla_particion %>%
  mutate(
    p_Bi_dado_A = (p_A_dado_Bi * p_Bi) / p_A_total  # Bayes: P(B_i|A)=P(A|B_i)P(B_i)/P(A)
  )

cat("=== TEOREMA DE BAYES CON MÁS DE DOS HIPÓTESIS (partición age_group) ===\n")
print(tabla_bayes_edad %>% select(age_group, p_Bi, p_A_dado_Bi, p_Bi_dado_A), n = Inf)
cat("\nSuma de P(B_i|A) sobre todos los grupos (debe ser ≈ 1)      :",
    round(sum(tabla_bayes_edad$p_Bi_dado_A), 6), "\n")
cat("Interpretación: dado que alguien se ha vacunado, p_Bi_dado_A indica\n")
cat("qué grupo de edad es más probable que sea, revirtiendo la dirección\n")
cat("de la información original P(A|B_i).\n\n")

# ============================================================
# 7. INTERPRETACIÓN ESCRITA DE RESULTADOS (generada dinámicamente)
# ============================================================
# IMPORTANTE: el texto se construye a partir de los objetos ya
# calculados arriba (p_vacuna, p_vacuna_dado_recc_bayes, resultados_chi,
# matriz_cor_behavioral, etc.), en lugar de usar cifras fijas escritas
# a mano. Así el resultado siempre refleja el resultado real de ejecutar
# el script sobre `datos`, evitando afirmaciones que puedan ser falsas
# si los números concretos cambian.

# --- 7.1 Impacto de la recomendación médica (Teorema de Bayes) ---

incremento_pp <- (p_vacuna_dado_recc_bayes - p_vacuna) * 100
razon_prob    <- p_vacuna_dado_recc_bayes / p_vacuna

texto_bayes <- sprintf(
  paste(
    "La probabilidad a priori de vacunarse contra la H1N1 en la muestra es de",
    "aproximadamente %.2f%% (P(V) = %.4f). Al condicionar este suceso a la",
    "existencia de una recomendacion del profesional sanitario (doctor_recc_h1n1 = 1),",
    "la probabilidad a posteriori pasa a %.2f%% (P(V|R) = %.4f), lo que supone un",
    "incremento absoluto de %.1f puntos porcentuales y una razon de probabilidades",
    "(posterior/prior) de %.2f. En ausencia de recomendacion medica, la probabilidad",
    "de vacunacion es de %.2f%% (P(V|R^c) = %.4f). La verosimilitud P(R|V) = %.4f",
    "indica en que medida la recomendacion se asocia con haberse vacunado."
  ),
  p_vacuna * 100, p_vacuna,
  p_vacuna_dado_recc_bayes * 100, p_vacuna_dado_recc_bayes,
  incremento_pp, razon_prob,
  p_vacuna_dado_no_recc * 100, p_vacuna_dado_no_recc,
  p_recc_dado_vacuna
)

cat("=== 7.1 INTERPRETACIÓN: Impacto de la recomendación médica ===\n")
cat(strwrap(texto_bayes, width = 80), sep = "\n")
cat("\n\n")

# --- 7.2 Independencia de los comportamientos preventivos ---

n_total_pares  <- nrow(resultados_chi)
n_pares_dep    <- sum(resultados_chi$independientes == "No (se rechaza H0)")
p_valor_max    <- max(resultados_chi$p_valor)

# Par de variables con mayor correlación (excluyendo la diagonal)
cor_off_diag <- matriz_cor_behavioral
diag(cor_off_diag) <- NA
idx_max   <- which(cor_off_diag == max(cor_off_diag, na.rm = TRUE), arr.ind = TRUE)[1, ]
var_max_1 <- rownames(cor_off_diag)[idx_max[1]]
var_max_2 <- colnames(cor_off_diag)[idx_max[2]]
r_max     <- cor_off_diag[idx_max[1], idx_max[2]]

texto_independencia <- sprintf(
  paste(
    "Los tests Chi-cuadrado de independencia rechazan la hipotesis nula de",
    "independencia en %d de los %d pares de variables behavioral_* evaluados",
    "(p-valor maximo observado entre todos los pares = %.5f). La matriz de",
    "correlacion muestra una estructura de asociacion positiva; el par con mayor",
    "asociacion es %s y %s (r = %.3f). Esto sugiere que los comportamientos",
    "preventivos no son sucesos independientes segun la muestra analizada, algo",
    "a tener en cuenta al elegir modelos que asuman independencia entre variables",
    "(p. ej. Naive Bayes)."
  ),
  n_pares_dep, n_total_pares, p_valor_max,
  var_max_1, var_max_2, r_max
)

cat("=== 7.2 INTERPRETACIÓN: Independencia de comportamientos preventivos ===\n")
cat(strwrap(texto_independencia, width = 80), sep = "\n")
cat("\n")