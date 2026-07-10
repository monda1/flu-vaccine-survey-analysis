library(dplyr)
library(ggplot2)
library(tidyr)

# --------------------------------------------------------------
# 1. Descriptiva univariante de variables categóricas/ordinales
#    sociodemográficas
# --------------------------------------------------------------
describir_categorica <- function(x, nombre) {
  tab <- table(x, useNA = "no")
  prop <- prop.table(tab)
  cat("\n===", nombre, "===\n")
  print(round(prop * 100, 2))
  cat("Moda:", names(tab)[which.max(tab)], "\n")
}

describir_categorica(datos$education, "Nivel educativo")
describir_categorica(datos$income_poverty, "Situación respecto al umbral de pobreza")
describir_categorica(datos$age_group, "Grupo de edad")

# --------------------------------------------------------------
# Tabla de frecuencias completa con acumuladas
# --------------------------------------------------------------
# El resumen de arriba solo da frecuencias relativas simples. Para una
# variable ORDINAL (como las Likert) tiene sentido añadir las
# frecuencias acumuladas: permiten responder preguntas del tipo
# "¿qué porcentaje de la gente está como mucho 'algo preocupada' (<=2)?",
# cosa que la frecuencia simple no puede contestar directamente.

tabla_frecuencias <- function(x, nombre) {
  x <- x[!is.na(x)]
  tab <- table(x)
  df <- data.frame(
    valor = names(tab),
    frecuencia_abs = as.vector(tab)
  )
  df$frecuencia_rel <- round(df$frecuencia_abs / sum(df$frecuencia_abs), 4)
  df$frec_abs_acum <- cumsum(df$frecuencia_abs)
  df$frec_rel_acum <- round(cumsum(df$frecuencia_rel), 4)
  cat("\n=== Tabla de frecuencias (con acumuladas):", nombre, "===\n")
  print(df)
  df
}

tabla_h1n1_concern <- tabla_frecuencias(datos$h1n1_concern, "h1n1_concern")
tabla_opinion_riesgo <- tabla_frecuencias(datos$opinion_h1n1_risk, "opinion_h1n1_risk")

# Polígono de frecuencias acumuladas (representación gráfica propia
# de una distribución acumulada, distinta del histograma de barras)
ggplot(tabla_h1n1_concern, aes(x = as.numeric(valor), y = frec_rel_acum)) +
  geom_line(linewidth = 1, color = "darkred") +
  geom_point(size = 2, color = "darkred") +
  scale_x_continuous(breaks = as.numeric(tabla_h1n1_concern$valor)) +
  labs(title = "Distribución de frecuencias RELATIVAS ACUMULADAS: h1n1_concern",
       x = "Nivel de preocupación", y = "Proporción acumulada") +
  theme_minimal()

# --------------------------------------------------------------
# 2. Variables Likert (ordinales): mediana y moda, NO media
#    : aritméticamente la media asume distancias
#    iguales entre categorías, algo no garantizado en una escala Likert.
#    Usamos mediana (con las categorías como factor ordenado) y moda.
# --------------------------------------------------------------
moda <- function(x) {
  x <- x[!is.na(x)]
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

variables_likert <- c("h1n1_concern", "h1n1_knowledge",
                      "opinion_h1n1_vacc_effective", "opinion_h1n1_risk",
                      "opinion_h1n1_sick_from_vacc", "opinion_seas_vacc_effective",
                      "opinion_seas_risk", "opinion_seas_sick_from_vacc")

resumen_likert <- data.frame(
  variable = variables_likert,
  mediana = sapply(variables_likert, function(v) median(datos[[v]], na.rm = TRUE)),
  moda = sapply(variables_likert, function(v) moda(datos[[v]])),
  # Se incluye la media SOLO a título comparativo, señalando su limitación
  media_referencial = sapply(variables_likert, function(v) round(mean(datos[[v]], na.rm = TRUE), 2))
)
cat("\ Resumen variables Likert (mediana y moda como medidas preferentes)")
print(resumen_likert)

# Visualización como distribución de frecuencias 
datos_likert_long <- datos %>%
  select(all_of(variables_likert)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "valor") %>%
  filter(!is.na(valor))

ggplot(datos_likert_long, aes(x = factor(valor))) +
  geom_bar(fill = "steelblue") +
  facet_wrap(~variable, scales = "free_y") +
  labs(title = "Distribución de frecuencias de variables tipo Likert",
       x = "Categoría de respuesta", y = "Frecuencia") +
  theme_minimal()

# --------------------------------------------------------------
# Media aritmética vs. geométrica vs. armónica
# --------------------------------------------------------------

comparar_medias <- function(x, nombre) {
  x <- x[!is.na(x)] + 1  # +1 para evitar log(0) y división por 0
  media_arit <- mean(x)
  media_geom <- exp(mean(log(x)))
  media_arm <- length(x) / sum(1 / x)
  
  cat("\n=== Comparación de medias:", nombre, "(+1 para evitar ceros) ===\n")
  cat("Media aritmética:", round(media_arit, 4), "\n")
  cat("Media geométrica :", round(media_geom, 4), "\n")
  cat("Media armónica   :", round(media_arm, 4), "\n")
  cat("Cumple AM >= GM >= HM:", media_arit >= media_geom && media_geom >= media_arm, "\n")
  cat("(a mayor dispersión relativa, mayor separación entre las tres medias)\n")
}

comparar_medias(datos$household_adults, "Adultos por hogar")
comparar_medias(datos$household_children, "Menores por hogar")

# --------------------------------------------------------------
# Cuantiles y diagrama de caja
# --------------------------------------------------------------
# La mediana ya se calculó arriba; aquí se completa con el resto de
# cuantiles (cuartiles, deciles) y su representación gráfica habitual,
# el boxplot, que muestra posición Y dispersión a la vez.

cuantiles_resumen <- function(x, nombre) {
  x <- x[!is.na(x)]
  cat("\n=== Cuantiles:", nombre, "===\n")
  cat("Cuartiles (Q1, Q2=mediana, Q3):\n")
  print(round(quantile(x, probs = c(0.25, 0.5, 0.75)), 3))
  cat("Rango intercuartílico (IQR):", round(IQR(x), 3), "\n")
  cat("Deciles:\n")
  print(round(quantile(x, probs = seq(0.1, 0.9, 0.1)), 3))
}

cuantiles_resumen(datos$household_adults, "Adultos por hogar")
cuantiles_resumen(datos$household_children, "Menores por hogar")

datos_hogar_long <- datos %>%
  select(household_adults, household_children) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "valor") %>%
  filter(!is.na(valor))

ggplot(datos_hogar_long, aes(x = variable, y = valor)) +
  geom_boxplot(fill = "lightblue") +
  labs(title = "Diagrama de caja: composición del hogar",
       x = NULL, y = "Número de personas") +
  theme_minimal()

# --------------------------------------------------------------
# Medidas de dispersión
# --------------------------------------------------------------
# dos variables con la misma media pueden ser radicalmente distintas
# (una muy homogénea, otra muy variable) . Se añade el coeficiente de variación (CV) porque,
# a diferencia de la varianza/sd, permite comparar la dispersión
# relativa entre dos variables con escalas o medias distintas.

dispersión_resumen <- function(x, nombre) {
  x <- x[!is.na(x)]
  rango <- max(x) - min(x)
  var_x <- var(x)
  sd_x <- sd(x)
  cv <- sd_x / mean(x) * 100  # en %
  
  cat("\n=== Medidas de dispersión:", nombre, "===\n")
  cat("Rango:", rango, "\n")
  cat("Varianza:", round(var_x, 4), "\n")
  cat("Desviación típica:", round(sd_x, 4), "\n")
  cat("Coeficiente de variación:", round(cv, 2), "%\n")
}

dispersión_resumen(datos$household_adults, "Adultos por hogar")
dispersión_resumen(datos$household_children, "Menores por hogar")

# El CV solo tiene sentido pleno con variables de razón (con cero
# absoluto significativo), que es el caso aquí (0 adultos/menores en el
# hogar es un cero real). No se aplica a las Likert por tratarse de
# escalas ordinales sin cero absoluto interpretable de la misma forma.

# =====================================================================
# Momentos, asimetría y curtosis
# =====================================================================
# dos distribuciones pueden tener la misma
# media y varianza pero forma distinta (una simétrica, otra sesgada;
# una con colas pesadas, otra ligera). El p-valor de un contraste o la
# media por sí solos no revelan esa forma.

momentos_forma <- function(x, nombre) {
  x <- x[!is.na(x)]
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  
  # Asimetría (coeficiente de Fisher, g1): 0 = simétrica;
  # >0 = cola a la derecha (valores altos poco frecuentes);
  # <0 = cola a la izquierda.
  asimetria <- (sum((x - m)^3) / n) / s^3
  
  # Curtosis (exceso sobre la normal, g2): 0 = como la normal;
  # >0 = leptocúrtica (colas más pesadas / pico más agudo que la normal);
  # <0 = platicúrtica (más aplanada que la normal).
  curtosis <- (sum((x - m)^4) / n) / s^4 - 3
  
  cat("\n=== Momentos de forma:", nombre, "===\n")
  cat("Asimetría (g1):", round(asimetria, 4),
      ifelse(asimetria > 0, "-> cola a la derecha", "-> cola a la izquierda"), "\n")
  cat("Curtosis (g2, exceso):", round(curtosis, 4),
      ifelse(curtosis > 0, "-> más apuntada que la normal (leptocúrtica)",
             "-> más aplanada que la normal (platicúrtica)"), "\n")
}

momentos_forma(datos$household_adults, "Adultos por hogar")
momentos_forma(datos$household_children, "Menores por hogar")

# --------------------------------------------------------------
# 3. Estadística bidimensional: tablas de contingencia y chi-cuadrado
#    entre variables sociodemográficas y vacunación
# --------------------------------------------------------------
analizar_asociacion <- function(datos, var_categorica, var_vacuna, nombre_var, nombre_vacuna) {
  df <- datos %>% filter(!is.na(.data[[var_categorica]]), !is.na(.data[[var_vacuna]]))
  
  tabla <- table(df[[var_categorica]], df[[var_vacuna]])
  tabla_prop <- prop.table(tabla, margin = 1)  # % vacunados dentro de cada categoría
  
  test <- chisq.test(tabla)
  
  # V de Cramér como medida de tamaño del efecto (complementa el p-valor)
  n <- sum(tabla)
  k <- min(dim(tabla))
  v_cramer <- sqrt(test$statistic / (n * (k - 1)))
  
  cat("\n=== ", nombre_var, " vs ", nombre_vacuna, " ===\n", sep = "")
  cat("Chi-cuadrado =", round(test$statistic, 2),
      "| gl =", test$parameter,
      "| p-valor =", format.pval(test$p.value, digits = 4), "\n")
  cat("V de Cramér =", round(v_cramer, 4), "\n")
  print(round(tabla_prop * 100, 2))
  
  invisible(list(tabla = tabla, test = test, v_cramer = v_cramer))
}

res_educ_h1n1 <- analizar_asociacion(datos, "education", "h1n1_vaccine", "Educación", "H1N1")
res_pobreza_h1n1 <- analizar_asociacion(datos, "income_poverty", "h1n1_vaccine", "Pobreza", "H1N1")
res_edad_h1n1 <- analizar_asociacion(datos, "age_group", "h1n1_vaccine", "Edad", "H1N1")

res_educ_seas <- analizar_asociacion(datos, "education", "seasonal_vaccine", "Educación", "Estacional")
res_pobreza_seas <- analizar_asociacion(datos, "income_poverty", "seasonal_vaccine", "Pobreza", "Estacional")
res_edad_seas <- analizar_asociacion(datos, "age_group", "seasonal_vaccine", "Edad", "Estacional")

# Tabla resumen de tamaños de efecto (V de Cramér) para comparar
# qué variable sociodemográfica se asocia más fuertemente con vacunarse
resumen_asociaciones <- data.frame(
  variable = rep(c("Educación", "Pobreza", "Edad"), 2),
  vacuna = rep(c("H1N1", "Estacional"), each = 3),
  v_cramer = c(res_educ_h1n1$v_cramer, res_pobreza_h1n1$v_cramer, res_edad_h1n1$v_cramer,
               res_educ_seas$v_cramer, res_pobreza_seas$v_cramer, res_edad_seas$v_cramer),
  p_valor = c(res_educ_h1n1$test$p.value, res_pobreza_h1n1$test$p.value, res_edad_h1n1$test$p.value,
              res_educ_seas$test$p.value, res_pobreza_seas$test$p.value, res_edad_seas$test$p.value)
)
cat("\n=== Resumen de asociaciones (V de Cramér) ===\n")
print(resumen_asociaciones %>% arrange(desc(v_cramer)))

# Visualización de una de las tablas como ejemplo
ggplot(datos %>% filter(!is.na(education), !is.na(h1n1_vaccine)),
       aes(x = education, fill = factor(h1n1_vaccine))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Proporción de vacunación H1N1 por nivel educativo",
       x = "Nivel educativo", y = "Proporción", fill = "Vacunado H1N1") +
  theme_minimal() +
  coord_flip()

# =====================================================================
# Distribuciones bidimensionales: marginales, condicionales,
# covarianza y correlación
# =====================================================================
# Las tablas de contingencia de arriba (educación x vacuna, etc.) son
# distribuciones bidimensionales, y tabla_prop (margin=1) ya es la
# distribución condicional de la vacuna dado el nivel educativo. Lo que
# falta es mostrar explícitamente la distribución marginal y, para las
# dos únicas variables genuinamente numéricas del dataset, la covarianza
# y correlación lineal (que no tienen sentido en variables categóricas
# como educación x vacuna, de ahí que ahí se use chi-cuadrado y no
# correlación).

tabla_bidim <- table(datos$household_adults, datos$household_children)
cat("\n=== Distribución bidimensional: adultos x menores en el hogar ===\n")
print(tabla_bidim)

cat("\nDistribución MARGINAL de adultos (sumando sobre menores):\n")
print(round(prop.table(margin.table(tabla_bidim, margin = 1)), 4))
cat("\nDistribución MARGINAL de menores (sumando sobre adultos):\n")
print(round(prop.table(margin.table(tabla_bidim, margin = 2)), 4))

cat("\nDistribución CONDICIONAL de menores dado 0 adultos en el hogar:\n")
print(round(prop.table(tabla_bidim["0", ]), 4))

# Covarianza y correlación lineal (Pearson)
cov_hogar <- cov(datos$household_adults, datos$household_children, use = "complete.obs")
cor_hogar <- cor(datos$household_adults, datos$household_children, use = "complete.obs")

cat("\nCovarianza (adultos, menores):", round(cov_hogar, 4), "\n")
cat("Correlación de Pearson (r):", round(cor_hogar, 4), "\n")
cat("Interpretación: r cercano a 0 indica que el número de adultos y de",
    "menores en el hogar son prácticamente INDEPENDIENTES entre sí en este dataset",
    "(no que no haya relación de ningún tipo, solo que no es lineal).\n")