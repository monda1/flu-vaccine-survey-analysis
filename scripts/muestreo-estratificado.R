# =====================================================================
#  Muestreo estratificado, estimadores y asignación
# Variable de estratificación: income_poverty (nivel de renta respecto
#   al umbral de pobreza del Censo 2008: Below Poverty, <= $75,000
#   Above Poverty, > $75,000)
# Variable objetivo: h1n1_vaccine (Tasa de vacunación H1N1)
# =====================================================================

library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)

# --------------------------------------------------------------
# 1. Preparación del marco muestral efectivo
# Filtramos los valores perdidos en el nivel de renta y la variable objetivo
# --------------------------------------------------------------
df_muestreo <- datos %>%
  filter(!is.na(income_poverty), !is.na(h1n1_vaccine))

n_total <- nrow(df_muestreo)

cat("=== Tamaño del marco muestral efectivo ===\n")
cat("N =", n_total, "observaciones útiles\n\n")

# --------------------------------------------------------------
# 2. Estimador Ingenuo: Muestreo Aleatorio Simple (MAS)
# Ignoramos la estructura por nivel de renta.
# p_mas = Estimación puntual de la proporción
# var_mas = Varianza del estimador MAS
# --------------------------------------------------------------
p_mas <- mean(df_muestreo$h1n1_vaccine)

# Asumimos una población infinita o muy grande respecto a la muestra,
# por lo que el Factor de Corrección por Población Finita (fpc) tiende a 1.
var_mas <- (p_mas * (1 - p_mas)) / (n_total - 1)

cat("=== Estimador Ingenuo (MAS) ===\n")
cat("Tasa global estimada (p_mas):", percent(p_mas, accuracy = 0.01), "\n")
cat("Varianza del estimador (MAS):", formatC(var_mas, format = "e", digits = 4), "\n\n")

# --------------------------------------------------------------
# 3. Estimadores Estratificados por Nivel de Renta (income_poverty)
# A falta de datos censales exactos (N_h), asumimos que el peso de
# cada estrato (W_h) es su peso muestral (n_h / n_total).
# --------------------------------------------------------------
estratos <- df_muestreo %>%
  group_by(income_poverty) %>%
  summarise(
    n_h = n(),                              # Tamaño muestral del estrato
    W_h = n_h / n_total,                    # Peso poblacional (aproximado)
    p_h = mean(h1n1_vaccine),               # Estimador puntual de la proporción
    var_h = (p_h * (1 - p_h)) / (n_h - 1),  # Varianza muestral intra-estrato
    S_h = sqrt(p_h * (1 - p_h)),            # Desviación cuasi-estándar del estrato
    .groups = "drop"
  )

# El estimador global estratificado es la suma ponderada de las proporciones:
p_estratificado <- sum(estratos$W_h * estratos$p_h)

# La varianza del estimador estratificado (asumiendo muestreo independiente por estrato):
# Var(p_st) = Sum(W_h^2 * var_h)
var_estratificada <- sum((estratos$W_h^2) * estratos$var_h)

cat("=== Estimador Estratificado ===\n")
cat("Tasa global (p_st):", percent(p_estratificado, accuracy = 0.01),
    "(Coincide con MAS al usar pesos muestrales)\n")
cat("Varianza del estimador (Estratificado):", formatC(var_estratificada, format = "e", digits = 4), "\n\n")

# --------------------------------------------------------------
# ¿POR QUÉ funciona la estratificación? Descomposición de varianza
# La varianza total de la variable objetivo se descompone en:
#   Varianza TOTAL = Varianza ENTRE estratos + Varianza DENTRO de estratos
# La estratificación solo "compensa" el esfuerzo si la varianza entre
# estratos es una parte sustancial de la varianza total (estratos
# internamente homogéneos, externamente heterogéneos).
#
# La razón de correlación eta^2 mide justamente esto: qué proporción de
# la variabilidad total de h1n1_vaccine es "explicada" por el nivel de renta.
# eta^2 alto  -> income_poverty es una buena variable de estratificación.
# eta^2 bajo  -> estratificar por nivel de renta aporta poca ganancia real.
# --------------------------------------------------------------
media_global <- mean(df_muestreo$h1n1_vaccine)
SST <- sum((df_muestreo$h1n1_vaccine - media_global)^2)

SSB <- estratos %>%
  summarise(SSB = sum(n_h * (p_h - media_global)^2)) %>%
  pull(SSB)

SSW <- SST - SSB
eta2 <- SSB / SST

cat("=== Descomposición de la varianza (ANOVA de un factor) ===\n")
cat("Suma de cuadrados TOTAL      (SST):", round(SST, 2), "\n")
cat("Suma de cuadrados ENTRE      (SSB):", round(SSB, 2), "\n")
cat("Suma de cuadrados DENTRO     (SSW):", round(SSW, 2), "\n")
cat("Razón de correlación (eta^2)      :", round(eta2, 4),
    "-> el", percent(eta2, accuracy = 0.1),
    "de la variabilidad de la vacunación se explica por el nivel de renta\n\n")

# --------------------------------------------------------------
# 4. Evaluación de la Ganancia de Precisión (Efecto de Diseño - Deff)
# --------------------------------------------------------------
ganancia_absoluta <- var_mas - var_estratificada
ganancia_relativa <- (ganancia_absoluta / var_mas) * 100
eficiencia_relativa <- var_mas / var_estratificada # Cuántas veces más varianza tiene MAS

# El Efecto de Diseño (Deff) es el cociente Var(diseño real) / Var(MAS
# con el mismo n). Deff < 1 indica que el diseño estratificado es más
# eficiente que el MAS. Su inverso (1/Deff) se interpreta como el
# "tamaño de muestra efectivo" relativo: n_efectivo = n_total / Deff,
# es decir, cuántas observaciones de un MAS necesitarías para igualar
# la precisión ya obtenida con el diseño estratificado real.
deff <- var_estratificada / var_mas
n_efectivo <- n_total / deff

cat("=== Evaluación de Eficiencia ===\n")
cat("Reducción absoluta de la varianza:", formatC(ganancia_absoluta, format = "e", digits = 4), "\n")
cat("Ganancia porcentual de precisión :", round(ganancia_relativa, 3), "%\n")
cat("Eficiencia Relativa (MAS/Est)    :", round(eficiencia_relativa, 4), "\n")
cat("(Valores > 1 indican que el estimador estratificado es más preciso)\n")
cat("Efecto de Diseño (Deff)          :", round(deff, 4), "\n")
cat("Tamaño de muestra efectivo       :", round(n_efectivo), "de", n_total, "observaciones reales\n")
cat("(Interpretación: con MAS necesitarías", round(n_efectivo),
    "obs. para igualar la precisión que ya da la estratificación con", n_total, ")\n\n")

# --------------------------------------------------------------
# Intervalo de confianza global para el estimador estratificado
# --------------------------------------------------------------
z_critico <- qnorm(0.975)
ic_inferior <- p_estratificado - z_critico * sqrt(var_estratificada)
ic_superior <- p_estratificado + z_critico * sqrt(var_estratificada)

cat("=== Intervalo de Confianza (95%) del estimador estratificado ===\n")
cat("[", percent(ic_inferior, accuracy = 0.01), ",",
    percent(ic_superior, accuracy = 0.01), "]\n\n")

# --------------------------------------------------------------
# 5. Asignación de Observaciones (Allocation)
# Supongamos que queremos diseñar una NUEVA encuesta con n = 10,000
# y queremos comparar cómo se distribuiría la muestra teóricamente.
# --------------------------------------------------------------
n_nuevo_diseno <- 10000
L <- nrow(estratos) # número de estratos

estratos_asignacion <- estratos %>%
  mutate(
    # 5.1 Asignación Igualitaria: mismo n_h en todos los estratos.
    # Sirve de referencia "ingenua" para comparar con las otras dos.
    n_igual = round(n_nuevo_diseno / L),
    
    # 5.2 Asignación Proporcional: n_h proporcional al peso del estrato W_h
    n_proporcional = round(n_nuevo_diseno * W_h),
    
    # 5.3 Asignación de Neyman (Óptima bajo costes iguales): n_h proporcional a W_h * S_h
    numerador_neyman = W_h * S_h,
    n_neyman = round(n_nuevo_diseno * (numerador_neyman / sum(numerador_neyman)))
  ) %>%
  select(income_poverty, W_h, p_h, S_h, n_igual, n_proporcional, n_neyman)

cat("=== Asignación para un nuevo diseño (n = 10,000) ===\n")
print(estratos_asignacion %>% select(-W_h, -S_h) %>% arrange(desc(p_h)))

# --------------------------------------------------------------
# Varianza teórica esperada bajo cada esquema de asignación
# Comparamos, para el mismo n = 10,000, qué varianza tendría el
# estimador global bajo cada regla de reparto. Esto es lo que
# justifica por qué Neyman es "óptima": minimiza Var(p_st) para
# un n total fijo, mientras que la igualitaria suele ser la peor.
# --------------------------------------------------------------
var_teorica <- function(n_h_vector) {
  sum((estratos$W_h^2) * (estratos$S_h^2) / n_h_vector)
}

comparacion_varianzas <- tibble(
  esquema = c("Igualitaria", "Proporcional", "Neyman"),
  varianza_esperada = c(
    var_teorica(estratos_asignacion$n_igual),
    var_teorica(estratos_asignacion$n_proporcional),
    var_teorica(estratos_asignacion$n_neyman)
  )
) %>%
  arrange(varianza_esperada)

cat("\n=== Varianza teórica esperada por esquema de asignación (n = 10,000) ===\n")
print(comparacion_varianzas)
cat("(El esquema con menor varianza esperada es el más eficiente; Neyman debería ganar)\n\n")

# --------------------------------------------------------------
# Asignación de Neyman generalizada con costes distintos por estrato
# En la práctica, encuestar a cada estrato de renta no cuesta lo mismo:
# localizar y entrevistar a hogares de renta baja (más rotación, menos
# contactabilidad) suele ser más costoso que a hogares de renta alta.
# La asignación óptima que minimiza la varianza para un PRESUPUESTO fijo
# (en vez de un n fijo) es:
#   n_h  ∝  (W_h * S_h) / sqrt(c_h)
# En lugar de solo W_h * S_h. Aquí simulamos costes relativos hipotéticos
# (puedes sustituir por costes reales de campo si los tienes).
# --------------------------------------------------------------
set.seed(123)
estratos_coste <- estratos %>%
  mutate(
    c_h = runif(n(), min = 1, max = 3), # coste relativo simulado por estrato
    numerador_neyman_coste = (W_h * S_h) / sqrt(c_h),
    n_neyman_coste = round(n_nuevo_diseno * (numerador_neyman_coste / sum(numerador_neyman_coste)))
  ) %>%
  select(income_poverty, W_h, S_h, c_h, n_neyman_coste)

cat("=== Asignación de Neyman con costes diferenciados por nivel de renta ===\n")
cat("(Costes simulados; sustituir por costes reales de campo si están disponibles)\n")
print(estratos_coste %>% arrange(desc(c_h)))
cat("\n")

# --------------------------------------------------------------
# 6. Visualizaciones
# --------------------------------------------------------------
# Heterogeneidad por nivel de renta (Tasa y Error Estándar)
ggplot(estratos, aes(x = reorder(income_poverty, p_h), y = p_h)) +
  geom_point(color = "darkred", size = 3) +
  geom_errorbar(aes(ymin = p_h - 1.96*sqrt(var_h), ymax = p_h + 1.96*sqrt(var_h)),
                width = 0.2, color = "darkred") +
  geom_hline(yintercept = p_mas, linetype = "dashed", color = "blue") +
  coord_flip() +
  scale_y_continuous(labels = percent) +
  labs(
    title = "Tasa de vacunación H1N1 por Nivel de Renta",
    subtitle = "Intervalos de confianza al 95%. Línea azul = Estimación global.",
    x = "Nivel de Renta (respecto al umbral de pobreza)",
    y = "Proporción de Vacunados (Estimación + IC 95%)"
  ) +
  theme_minimal()

# Asignación Proporcional vs Neyman vs Igualitaria
df_plot_asignacion <- estratos_asignacion %>%
  select(income_poverty, n_igual, n_proporcional, n_neyman) %>%
  pivot_longer(cols = c(n_igual, n_proporcional, n_neyman),
               names_to = "tipo_asignacion", values_to = "n_observaciones")

ggplot(df_plot_asignacion, aes(x = income_poverty, y = n_observaciones, fill = tipo_asignacion)) +
  geom_col(position = "dodge") +
  labs(
    title = "Comparación de Asignación Muestral para n = 10,000",
    subtitle = "Igualitaria vs Proporcional vs Óptima (Neyman)",
    x = "Nivel de Renta",
    y = "Número de observaciones a asignar",
    fill = "Método de Asignación"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Varianza esperada por esquema de asignación
ggplot(comparacion_varianzas, aes(x = reorder(esquema, varianza_esperada), y = varianza_esperada)) +
  geom_col(fill = "steelblue") +
  labs(
    title = "Varianza teórica esperada según el esquema de asignación",
    subtitle = "n = 10,000 | Menor varianza = mayor eficiencia del diseño",
    x = "Esquema de asignación",
    y = "Varianza esperada del estimador global"
  ) +
  theme_minimal()