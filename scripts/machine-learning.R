library(gbm)
library(pROC)

# ==========================================
# 1. Ingeneria de características
# ==========================================

crear_features_elite <- function(df) {
  # A) Score de Conciencia Sanitaria (Conductas preventivas)
  cols_conducta <- c("behavioral_antiviral_meds", "behavioral_avoidance", 
                     "behavioral_face_mask", "behavioral_wash_hands", 
                     "behavioral_large_gatherings", "behavioral_outside_home", 
                     "behavioral_touch_face")
  df$behavioral_score <- rowSums(df[, cols_conducta], na.rm = TRUE)
  
  # B) INTERACCIÓN MÉDICA (Súper Predictores Cruzados)
  dr_h1n1 <- ifelse(is.na(df$doctor_recc_h1n1), 0, df$doctor_recc_h1n1)
  dr_seas <- ifelse(is.na(df$doctor_recc_seasonal), 0, df$doctor_recc_seasonal)
  df$doctor_recc_any  <- as.factor(ifelse(dr_h1n1 == 1 | dr_seas == 1, 1, 0))
  df$doctor_recc_both <- as.factor(ifelse(dr_h1n1 == 1 & dr_seas == 1, 1, 0))
  
  # C) Índice de Vulnerabilidad Física y Exposición
  chronic <- ifelse(is.na(df$chronic_med_condition), 0, df$chronic_med_condition)
  infant  <- ifelse(is.na(df$child_under_6_months), 0, df$child_under_6_months)
  worker  <- ifelse(is.na(df$health_worker), 0, df$health_worker)
  df$vulnerability_score <- chronic + infant + worker
  
  # D) Índices Psicológicos de Percepción de Riesgo
  eff_h1n1  <- ifelse(is.na(df$opinion_h1n1_vacc_effective), 3, df$opinion_h1n1_vacc_effective)
  risk_h1n1 <- ifelse(is.na(df$opinion_h1n1_risk), 3, df$opinion_h1n1_risk)
  sick_h1n1 <- ifelse(is.na(df$opinion_h1n1_sick_from_vacc), 3, df$opinion_h1n1_sick_from_vacc)
  df$psych_index_h1n1 <- eff_h1n1 + risk_h1n1 - (sick_h1n1 * 0.5)
  
  eff_seas  <- ifelse(is.na(df$opinion_seas_vacc_effective), 3, df$opinion_seas_vacc_effective)
  risk_seas <- ifelse(is.na(df$opinion_seas_risk), 3, df$opinion_seas_risk)
  sick_seas <- ifelse(is.na(df$opinion_seas_sick_from_vacc), 3, df$opinion_seas_sick_from_vacc)
  df$psych_index_seasonal <- eff_seas + risk_seas - (sick_seas * 0.5)
  
  # E) Proxy de Estabilidad (Seguro Médico + Vivienda propia)
  ins <- ifelse(is.na(df$health_insurance) | df$health_insurance == 0, 0, 1)
  own <- ifelse(!is.na(df$rent_or_own) & df$rent_or_own == "Own", 1, 0)
  df$stability_proxy <- ins + own
  
  return(df)
}

datos <- crear_features_elite(datos)
test  <- crear_features_elite(test)

targets      <- c("h1n1_vaccine", "seasonal_vaccine")
ignore_cols  <- c("respondent_id", "en_complete_case")
feature_cols <- setdiff(colnames(datos), c(targets, ignore_cols))

# Alineación estricta de factores
is_char    <- sapply(datos[, feature_cols, drop = FALSE], is.character)
char_names <- feature_cols[is_char]

if (length(char_names) > 0) {
  for (col in char_names) {
    datos[[col]] <- as.factor(datos[[col]])
    if (col %in% colnames(test)) {
      test[[col]] <- factor(test[[col]], levels = levels(datos[[col]]))
    }
  }
}

# ==========================================
# 2. PARTICIÓN DE VALIDACIÓN INDEPENDIENTE (80/20)
# ==========================================
set.seed(42) 
train_idx  <- sample(seq_len(nrow(datos)), size = 0.8 * nrow(datos))
train_data <- datos[train_idx, ]
val_data   <- datos[-train_idx, ]

formula_h1n1     <- as.formula(paste("h1n1_vaccine ~", paste(feature_cols, collapse = " + ")))
formula_seasonal <- as.formula(paste("seasonal_vaccine ~", paste(feature_cols, collapse = " + ")))

# Ajustes robustos anti-overfitting
total_trees <- 1500  
chunk_size  <- 50    
pasos       <- total_trees / chunk_size


# ==========================================
# 3. ENTRENAMIENTO REGULARIZADO: H1N1 VACCINE
# ==========================================
message("\n--- [1/2] Entrenando modelo H1N1 (Anti-Overfitting) ---")
pb_h1n1 <- txtProgressBar(min = 0, max = total_trees, style = 3)

model_h1n1 <- gbm(
  formula = formula_h1n1, data = train_data, distribution = "bernoulli",
  n.trees = chunk_size, 
  interaction.depth = 4,       
  shrinkage = 0.01,            
  bag.fraction = 0.75,         
  n.minobsinnode = 15,         
  train.fraction = 0.80,       
  keep.data = TRUE, verbose = FALSE
)
setTxtProgressBar(pb_h1n1, chunk_size)

for (i in 2:pasos) {
  model_h1n1 <- gbm.more(model_h1n1, n.new.trees = chunk_size)
  setTxtProgressBar(pb_h1n1, i * chunk_size)
}
close(pb_h1n1)

best_trees_h1n1 <- gbm.perf(model_h1n1, method = "test", plot.it = FALSE)
cat(sprintf("-> Árbol óptimo seguro seleccionado para H1N1: %d\n", best_trees_h1n1))


# ==========================================
# 4. ENTRENAMIENTO REGULARIZADO: SEASONAL VACCINE
# ==========================================
message("\n--- [2/2] Entrenando modelo Seasonal (Anti-Overfitting) ---")
pb_seasonal <- txtProgressBar(min = 0, max = total_trees, style = 3)

model_seasonal <- gbm(
  formula = formula_seasonal, data = train_data, distribution = "bernoulli",
  n.trees = chunk_size, 
  interaction.depth = 4, 
  shrinkage = 0.01, 
  bag.fraction = 0.75, 
  n.minobsinnode = 15, 
  train.fraction = 0.80, 
  keep.data = TRUE, verbose = FALSE
)
setTxtProgressBar(pb_seasonal, chunk_size)

for (i in 2:pasos) {
  model_seasonal <- gbm.more(model_seasonal, n.new.trees = chunk_size)
  setTxtProgressBar(pb_seasonal, i * chunk_size)
}
close(pb_seasonal)

best_trees_seasonal <- gbm.perf(model_seasonal, method = "test", plot.it = FALSE)
cat(sprintf("-> Árbol óptimo seguro seleccionado para Seasonal: %d\n", best_trees_seasonal))


preds_h1n1     <- predict(model_h1n1, newdata = val_data, n.trees = best_trees_h1n1, type = "response")
preds_seasonal <- predict(model_seasonal, newdata = val_data, n.trees = best_trees_seasonal, type = "response")

auc_h1n1     <- as.numeric(auc(roc(val_data$h1n1_vaccine, preds_h1n1, quiet = TRUE)))
auc_seasonal <- as.numeric(auc(roc(val_data$seasonal_vaccine, preds_seasonal, quiet = TRUE)))
mean_auc     <- (auc_h1n1 + auc_seasonal) / 2

cat("\n=========================================\n")
cat("      RESULTADO LOCAL  \n")
cat("=========================================\n")
cat(sprintf("H1N1 ROC AUC:         %.4f\n", auc_h1n1))
cat(sprintf("Seasonal ROC AUC:     %.4f\n", auc_seasonal))
cat(sprintf("--> NUEVA ESTIMACIÓN MEDIA LOCAL: %.4f\n", mean_auc))
cat("=========================================\n")

# --- BLOQUE DE INTERPRETABILIDAD  ---
cat("\n=========================================\n")
cat("        INTERPRETABILIDAD      \n")
cat("=========================================\n")

importancia_h1n1     <- summary.gbm(model_h1n1, n.trees = best_trees_h1n1, plotit = FALSE)
importancia_seasonal <- summary.gbm(model_seasonal, n.trees = best_trees_seasonal, plotit = FALSE)

cat("\n[Top Variables Clave para H1N1]:\n")
print(head(importancia_h1n1, 4))

cat("\n[Top Variables Clave para Vacuna Estacional]:\n")
print(head(importancia_seasonal, 4))

# Ventana gráfica: 2 de importancia arriba, 2 de dependencia parcial (PDP) abajo
par(mfrow = c(2, 2))

# 1. Gráfico de Importancia H1N1
summary.gbm(model_h1n1, n.trees = best_trees_h1n1, main = "Importancia: H1N1", las = 2)

# 2. Gráfico de Importancia Estacional
summary.gbm(model_seasonal, n.trees = best_trees_seasonal, main = "Importancia: Estacional", las = 2)

# 3. PDP del impacto marginal de nuestro Índice Psicológico en H1N1
message("\nGenerando gráficos de Dependencia Parcial (PDP)...")
plot.gbm(model_h1n1, i.var = "psych_index_h1n1", n.trees = best_trees_h1n1,
         main = "PDP: Índice Psicológico (H1N1)", xlab = "Índice", ylab = "Efecto Log-Odds")

# 4. PDP del impacto marginal de nuestro Índice Psicológico en Vacuna Estacional
plot.gbm(model_seasonal, i.var = "psych_index_seasonal", n.trees = best_trees_seasonal,
         main = "PDP: Índice Psicológico (Estacional)", xlab = "Índice", ylab = "Efecto Log-Odds")

# Restaurar panel gráfico original
par(mfrow = c(1, 1))


message("\nGenerando predicciones  sobre el set 'test'...")
test_preds_h1n1     <- predict(model_h1n1, newdata = test, n.trees = best_trees_h1n1, type = "response")
test_preds_seasonal <- predict(model_seasonal, newdata = test, n.trees = best_trees_seasonal, type = "response")

submission <- data.frame(
  respondent_id    = test$respondent_id,
  h1n1_vaccine     = test_preds_h1n1,
  seasonal_vaccine = test_preds_seasonal
)

write.csv(submission, "submission.csv", row.names = FALSE)
