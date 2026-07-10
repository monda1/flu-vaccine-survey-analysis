library(dplyr)
library(car)        
library(broom)       
library(ggplot2)

# =====================================================================
# Regresión logística múltiple y multicolinealidad (VIF)
# =====================================================================

vars_modelo <- c("h1n1_vaccine",
                 "h1n1_concern", "h1n1_knowledge",
                 "opinion_h1n1_vacc_effective", "opinion_h1n1_risk", "opinion_h1n1_sick_from_vacc",
                 "opinion_seas_vacc_effective", "opinion_seas_risk", "opinion_seas_sick_from_vacc",
                 "doctor_recc_h1n1", "doctor_recc_seasonal",
                 "chronic_med_condition", "health_worker", "health_insurance",
                 "age_group", "education", "income_poverty", "sex")

datos_modelo <- datos %>%
  select(all_of(vars_modelo)) %>%
  na.omit() %>%
  mutate(
    age_group = factor(age_group),
    education = factor(education),
    income_poverty = factor(income_poverty),
    sex = factor(sex),
    h1n1_vaccine = factor(h1n1_vaccine)
  )

cat("Observaciones tras eliminar NA:", nrow(datos_modelo),
    "de", nrow(datos), "\n\n")

# --------------------------------------------------------------
# 2. Ajuste del modelo de regresión logística múltiple
# --------------------------------------------------------------
modelo_logit <- glm(
  h1n1_vaccine ~ h1n1_concern + h1n1_knowledge +
    opinion_h1n1_vacc_effective + opinion_h1n1_risk + opinion_h1n1_sick_from_vacc +
    opinion_seas_vacc_effective + opinion_seas_risk + opinion_seas_sick_from_vacc +
    doctor_recc_h1n1 + doctor_recc_seasonal +
    chronic_med_condition + health_worker + health_insurance +
    age_group + education + income_poverty + sex,
  data = datos_modelo,
  family = binomial(link = "logit")
)

cat("=== Resumen del modelo ===\n")
print(summary(modelo_logit))

# Odds ratios con IC 95% (más interpretables que los coeficientes en log-odds)
odds_ratios <- tidy(modelo_logit, exponentiate = TRUE, conf.int = TRUE) %>%
  arrange(desc(abs(estimate - 1)))

cat("\n=== Odds Ratios (ordenados por magnitud del efecto) ===\n")
print(odds_ratios, n = Inf)

# --------------------------------------------------------------
# Efectos marginales medios (AME): odds ratio en escala de
# probabilidad
# --------------------------------------------------------------
# Los odds ratios no son comparables directamente entre sí cuando los
# predictores tienen escalas distintas (binario vs. Likert 1-5 vs.
# categórico), y mucha gente los malinterpreta como si fueran "veces
# más probable" cuando en realidad son razón de ODDS, no de
# probabilidades. El Efecto Marginal Promedio (AME) traduce el efecto
# de cada variable a puntos porcentuales de probabilidad de vacunarse,
# manteniendo el resto de variables en sus valores observados
# (enfoque "promedio de efectos individuales", más robusto que fijar
# el resto en la media).

calcular_ame_numerica <- function(modelo, datos, variable, delta = 1) {
  datos_mas <- datos
  datos_mas[[variable]] <- datos_mas[[variable]] + delta
  p_base <- predict(modelo, newdata = datos, type = "response")
  p_mas <- predict(modelo, newdata = datos_mas, type = "response")
  mean(p_mas - p_base, na.rm = TRUE)
}

ame_concern <- calcular_ame_numerica(modelo_logit, datos_modelo, "h1n1_concern")
ame_riesgo_h1n1 <- calcular_ame_numerica(modelo_logit, datos_modelo, "opinion_h1n1_risk")

cat("\n=== Efectos marginales promedio (AME), en puntos de probabilidad ===\n")
cat("Subir 1 punto en h1n1_concern cambia la prob. de vacunarse en:",
    round(ame_concern * 100, 2), "p.p. (en promedio, ceteris paribus)\n")
cat("Subir 1 punto en opinion_h1n1_risk cambia la prob. de vacunarse en:",
    round(ame_riesgo_h1n1 * 100, 2), "p.p.\n")
cat("(A diferencia del odds ratio, esto SÍ es directamente comparable entre",
    "variables y se interpreta en la escala que a la gente le importa: probabilidad real)\n")

# --------------------------------------------------------------
# 3. Bondad de ajuste global
# --------------------------------------------------------------
pseudo_r2 <- 1 - (modelo_logit$deviance / modelo_logit$null.deviance)
cat("\nPseudo R² (McFadden):", round(pseudo_r2, 4), "\n")

# AUC como medida de capacidad predictiva global
pred_prob <- predict(modelo_logit, type = "response")

if (requireNamespace("pROC", quietly = TRUE)) {
  library(pROC)
  roc_obj <- roc(datos_modelo$h1n1_vaccine, pred_prob, quiet = TRUE)
  cat("AUC:", round(auc(roc_obj), 4), "\n")
  
  # La curva ROC completa aporta más que el número solo: muestra el
  # COMPROMISO entre sensibilidad y especificidad en cada punto de
  # corte, no solo un resumen escalar.
  ggroc(roc_obj, color = "darkred", linewidth = 1) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "gray50") +
    labs(title = "Curva ROC: modelo logístico H1N1",
         subtitle = paste0("AUC = ", round(auc(roc_obj), 3)),
         x = "Especificidad", y = "Sensibilidad") +
    theme_minimal()
}

# --------------------------------------------------------------
# Test de Hosmer-Lemeshow: bondad de ajuste específica
# de regresión logística
# --------------------------------------------------------------
# El pseudo R² de McFadden no tiene una interpretación de "% de
# varianza explicada" como el R² lineal, y valores bajos (0.1-0.3) son
# normales en logística y no implican mal ajuste. El test de
# Hosmer-Lemeshow sí evalúa directamente si las probabilidades
# predichas coinciden con las tasas observadas, agrupando a los
# individuos en deciles de probabilidad predicha.
# H0: el modelo ajusta bien (predicciones ~ observado por grupos)
# Un p-valor bajo indica falta de ajuste, no al revés.

if (requireNamespace("ResourceSelection", quietly = TRUE)) {
  library(ResourceSelection)
  hl_test <- hoslem.test(as.numeric(as.character(datos_modelo$h1n1_vaccine)),
                         pred_prob, g = 10)
  cat("\n=== Test de Hosmer-Lemeshow (bondad de ajuste) ===\n")
  cat("Estadístico:", round(hl_test$statistic, 4),
      "| p-valor:", format.pval(hl_test$p.value, digits = 4), "\n")
  cat(ifelse(hl_test$p.value > 0.05,
             "-> No hay evidencia de falta de ajuste (p > 0.05)\n",
             "-> El modelo podría no ajustar bien en algunos rangos de probabilidad (p <= 0.05)\n"))
}

# --------------------------------------------------------------
# Diagnóstico de multicolinealidad: VIF
#    Regla práctica habitual: VIF > 5 indica multicolinealidad
#    preocupante; VIF > 10, problemática. Con GVIF (para factores con
#    más de 2 niveles) se usa GVIF^(1/(2*gl)) como equivalente comparable.
# --------------------------------------------------------------
vif_resultado <- vif(modelo_logit)
print(vif_resultado)

# Normalizamos a un formato comparable independientemente de si son
# variables numéricas (VIF) o factores (GVIF con grados de libertad)
if (is.matrix(vif_resultado)) {
  vif_df <- data.frame(
    variable = rownames(vif_resultado),
    VIF_ajustado = vif_resultado[, "GVIF^(1/(2*Df))"]^2  # elevamos al cuadrado para comparar en escala de VIF
  )
} else {
  vif_df <- data.frame(
    variable = names(vif_resultado),
    VIF_ajustado = vif_resultado
  )
}

vif_df <- vif_df %>% arrange(desc(VIF_ajustado))
cat("\n=== VIF ajustado (comparable), ordenado de mayor a menor ===\n")
print(vif_df)

# Visualización
ggplot(vif_df, aes(x = reorder(variable, VIF_ajustado), y = VIF_ajustado)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "orange") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "Factor de Inflación de la Varianza (VIF) por predictor",
       subtitle = "Líneas de referencia: naranja = 5, rojo = 10",
       x = NULL, y = "VIF ajustado") +
  theme_minimal()

# --------------------------------------------------------------
# Índice de condición: diagnóstico de colinealidad multivariante
# --------------------------------------------------------------
# El VIF detecta colinealidad de cada variable contra el resto tomadas
# en conjunto, pero puede pasar por alto combinaciones lineales que
# involucran a varias variables a la vez sin que ninguna individualmente
# tenga un VIF alto. El índice de condición (Belsley, Kuh y Welsch),
# basado en los autovalores de la matriz de correlación de los
# predictores numéricos, detecta este tipo de colinealidad estructural.
# Regla: índice de condición > 30 sugiere colinealidad grave.

vars_numericas_modelo <- c("h1n1_concern", "h1n1_knowledge",
                           "opinion_h1n1_vacc_effective", "opinion_h1n1_risk", "opinion_h1n1_sick_from_vacc",
                           "opinion_seas_vacc_effective", "opinion_seas_risk", "opinion_seas_sick_from_vacc")

cor_matrix_num <- cor(datos_modelo[vars_numericas_modelo])
autovalores <- eigen(cor_matrix_num)$values
indice_condicion <- sqrt(max(autovalores) / autovalores)

cat("\n=== Índice de condición (colinealidad multivariante) ===\n")
cat("Autovalores:", round(autovalores, 4), "\n")
cat("Índices de condición:", round(indice_condicion, 2), "\n")
cat(ifelse(max(indice_condicion) > 30,
           "-> Colinealidad estructural relevante entre varias variables de opinión a la vez.\n",
           "-> No se detecta colinealidad multivariante grave más allá de lo ya visto en el VIF.\n"))

# --------------------------------------------------------------
# 5. Matriz de correlación entre las variables de opinión
# --------------------------------------------------------------
vars_opinion <- c("opinion_h1n1_vacc_effective", "opinion_h1n1_risk", "opinion_h1n1_sick_from_vacc",
                  "opinion_seas_vacc_effective", "opinion_seas_risk", "opinion_seas_sick_from_vacc")

cor_opinion <- cor(datos_modelo[vars_opinion], use = "complete.obs", method = "spearman")

cat("\n=== Correlación (Spearman, apropiada para variables ordinales) ===\n")
print(round(cor_opinion, 3))

if (requireNamespace("corrplot", quietly = TRUE)) {
  library(corrplot)
  corrplot(cor_opinion, method = "color", type = "upper",
           addCoef.col = "black", tl.col = "black", tl.srt = 45,
           title = "Correlación entre variables de opinión", mar = c(0,0,2,0))
}

# --------------------------------------------------------------
# Interacción: ¿pesa igual la recomendación médica en
# sanitarios que en el resto?
# --------------------------------------------------------------
# Cabría esperar que un
# sanitario, por su propio conocimiento profesional, dependa MENOS de
# que el médico se lo recomiende explícitamente para vacunarse (su
# decisión ya estaría tomada de antemano), mientras que en el resto de
# la población el "empujón" del médico podría ser más determinante.
# Esto es precisamente lo que testa un término de interacción, algo
# que ningún modelo aditivo (como los de arriba) puede capturar.

modelo_interaccion <- glm(
  h1n1_vaccine ~ doctor_recc_h1n1 * health_worker +
    h1n1_concern + chronic_med_condition + health_insurance + age_group,
  data = datos_modelo,
  family = binomial(link = "logit")
)

cat("\n=== Modelo con interacción doctor_recc_h1n1 x health_worker ===\n")
print(summary(modelo_interaccion)$coefficients)

# Test de razón de verosimilitud: ¿la interacción añade algo
# significativo frente al modelo sin ella? (estos dos modelos SÍ están
# anidados de verdad, a diferencia del caso completo vs. reducido de
# más abajo, así que el test es directamente aplicable)
modelo_sin_interaccion <- update(modelo_interaccion, . ~ . - doctor_recc_h1n1:health_worker)
lr_test_interaccion <- anova(modelo_sin_interaccion, modelo_interaccion, test = "Chisq")

cat("\nTest de razón de verosimilitud (¿importa la interacción?):\n")
print(lr_test_interaccion)
cat("Si el p-valor es significativo, el efecto de la recomendación médica",
    "SÍ depende de si la persona es sanitaria o no; un modelo puramente",
    "aditivo estaría ocultando esa diferencia.\n")

# --------------------------------------------------------------
# 6. Modelo alternativo reducido: colapsar variables de opinión
#    redundantes en un índice compuesto, como estrategia para mitigar
#    la multicolinealidad sin perder toda la información que aportan
# --------------------------------------------------------------
datos_modelo <- datos_modelo %>%
  mutate(
    indice_percepcion_riesgo = (opinion_h1n1_risk + opinion_seas_risk) / 2,
    indice_percepcion_eficacia = (opinion_h1n1_vacc_effective + opinion_seas_vacc_effective) / 2
  )

modelo_reducido <- glm(
  h1n1_vaccine ~ h1n1_concern + h1n1_knowledge +
    indice_percepcion_riesgo + indice_percepcion_eficacia +
    doctor_recc_h1n1 + doctor_recc_seasonal +
    chronic_med_condition + health_worker + health_insurance +
    age_group + education + income_poverty + sex,
  data = datos_modelo,
  family = binomial(link = "logit")
)

cat("\n=== Modelo reducido (índices compuestos en vez de opiniones individuales) ===\n")
print(summary(modelo_reducido))
cat("\nVIF del modelo reducido:\n")
print(vif(modelo_reducido))

# Comparación de ambos modelos: ¿se pierde mucha capacidad predictiva
# al simplificar, a cambio de ganar interpretabilidad?
cat("\n=== Comparación de modelos ===\n")
cat("Pseudo R² modelo completo:", round(pseudo_r2, 4), "\n")
pseudo_r2_reducido <- 1 - (modelo_reducido$deviance / modelo_reducido$null.deviance)
cat("Pseudo R² modelo reducido:", round(pseudo_r2_reducido, 4), "\n")

# --------------------------------------------------------------
#  Comparación anidada CORRECTA: ¿aportan algo las opiniones
# como grupo, frente a un modelo sin ellas?
# --------------------------------------------------------------
# Este sí es un test de razón de verosimilitud válido: el modelo sin
# opiniones es un subconjunto exacto de variables del modelo completo
# (no una transformación), así que ambos SÍ están anidados.

modelo_sin_opiniones <- update(
  modelo_logit,
  . ~ . - opinion_h1n1_vacc_effective - opinion_h1n1_risk - opinion_h1n1_sick_from_vacc -
    opinion_seas_vacc_effective - opinion_seas_risk - opinion_seas_sick_from_vacc
)

lr_test_opiniones <- anova(modelo_sin_opiniones, modelo_logit, test = "Chisq")
cat("\n=== Test de razón de verosimilitud: ¿aportan las opiniones como grupo? ===\n")
print(lr_test_opiniones)
cat("H0: las 6 variables de opinión no aportan explicación conjunta adicional.\n",
    "Un p-valor bajo rechaza H0: pese a la colinealidad entre ellas (vista en el VIF",
    "y el índice de condición), como GRUPO sí añaden información relevante al modelo.\n")

# --------------------------------------------------------------
# 7. Regresión Ridge: alternativa a colapsar variables para
# lidiar con la multicolinealidad
# --------------------------------------------------------------
# El índice compuesto (paso 6) reduce colinealidad a costa de perder
# matices entre las variables originales. La regresión Ridge ofrece
# otra vía: en vez de eliminar o combinar variables, penaliza los
# coeficientes (shrinkage) para estabilizarlos cuando hay colinealidad,
# manteniendo todas las variables originales en el modelo.

if (requireNamespace("glmnet", quietly = TRUE)) {
  library(glmnet)
  
  x_ridge <- model.matrix(modelo_logit)[, -1]  # sin el intercepto
  y_ridge <- datos_modelo$h1n1_vaccine
  
  set.seed(123)
  cv_ridge <- cv.glmnet(x_ridge, y_ridge, family = "binomial", alpha = 0)  # alpha=0 -> Ridge
  
  cat("\n=== Regresión Ridge (alpha=0): comparación de coeficientes ===\n")
  cat("Lambda óptimo (validación cruzada):", round(cv_ridge$lambda.min, 5), "\n")
  
  coef_ridge <- as.matrix(coef(cv_ridge, s = "lambda.min"))
  coef_glm <- coef(modelo_logit)
  
  comparacion_coef <- data.frame(
    variable = rownames(coef_ridge),
    coef_glm_original = coef_glm[rownames(coef_ridge)],
    coef_ridge = coef_ridge[, 1]
  ) %>% filter(!is.na(coef_glm_original))
  
  print(comparacion_coef)
  cat("\nInterpretación: Ridge 'encoge' los coeficientes hacia 0 respecto al GLM",
      "original, más agresivamente en las variables MÁS correlacionadas entre sí",
      "(las de opinión). Esto reduce la varianza de las estimaciones a cambio de",
      "un pequeño sesgo -- el compromiso clásico sesgo-varianza frente a la",
      "inestabilidad que la colinealidad introduce en el GLM sin penalizar.\n")
}