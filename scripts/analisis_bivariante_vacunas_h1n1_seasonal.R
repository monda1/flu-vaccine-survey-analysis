library(dplyr)
library(ggplot2)

# ============================================================
#  Distribución conjunta, marginales,
# condicionadas, covarianza, correlación e independencia (chi-cuadrado)
# ============================================================

# --- 1. Tabla de contingencia (distribución conjunta) ---
tabla_conjunta <- table(H1N1 = datos$h1n1_vaccine, Estacional = datos$seasonal_vaccine)
print(tabla_conjunta)

# Distribución conjunta en proporciones (probabilidad conjunta)
prop_conjunta <- prop.table(tabla_conjunta)
cat("\n--- Distribución conjunta (probabilidades) ---\n")
print(round(prop_conjunta, 4))

# --- 2. Distribuciones marginales ---
marginal_h1n1 <- prop.table(table(datos$h1n1_vaccine))
marginal_seasonal <- prop.table(table(datos$seasonal_vaccine))

cat("\n--- Marginal H1N1 ---\n")
print(round(marginal_h1n1, 4))
cat("\n--- Marginal Estacional ---\n")
print(round(marginal_seasonal, 4))

# --- 3. Distribuciones condicionales ---
# P(Estacional | H1N1)
cond_estacional_dado_h1n1 <- prop.table(tabla_conjunta, margin = 1)
cat("\n--- P(Estacional | H1N1) ---\n")
print(round(cond_estacional_dado_h1n1, 4))

# P(H1N1 | Estacional)
cond_h1n1_dado_estacional <- prop.table(tabla_conjunta, margin = 2)
cat("\n--- P(H1N1 | Estacional) ---\n")
print(round(cond_h1n1_dado_estacional, 4))

# --- 4. Covarianza y correlación ---
covarianza <- cov(datos$h1n1_vaccine, datos$seasonal_vaccine, use = "complete.obs")
correlacion <- cor(datos$h1n1_vaccine, datos$seasonal_vaccine, use = "complete.obs")

cat("\nCovarianza:", round(covarianza, 4), "\n")
cat("Correlación (coef. phi):", round(correlacion, 4), "\n")

# --- 5. Test de independencia chi-cuadrado ---
test_chi <- chisq.test(tabla_conjunta)
print(test_chi)

# --- 6. Visualización rápida ---
df_plot <- as.data.frame(prop_conjunta)
ggplot(df_plot, aes(x = H1N1, y = Estacional, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = scales::percent(Freq, accuracy = 0.1))) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Distribución conjunta: Vacuna H1N1 vs Estacional",
       x = "Vacuna H1N1 (0=No, 1=Sí)",
       y = "Vacuna Estacional (0=No, 1=Sí)",
       fill = "Proporción") +
  theme_minimal()


# ============================================================
# Modelo Bernoulli: esperanza, varianza y momentos
# ============================================================

n <- nrow(datos)

# Estimadores puntuales de p (parámetro de cada Bernoulli)
p_h1n1 <- mean(datos$h1n1_vaccine)
p_seasonal <- mean(datos$seasonal_vaccine)

# Esperanza: E[X] = p  (para una Bernoulli)
cat("\n--- Esperanza matemática E[X] (Bernoulli) ---\n")
cat("E[H1N1]      =", round(p_h1n1, 4), "\n")
cat("E[Estacional]=", round(p_seasonal, 4), "\n")

# Varianza: Var(X) = p(1-p)
var_h1n1 <- p_h1n1 * (1 - p_h1n1)
var_seasonal <- p_seasonal * (1 - p_seasonal)
cat("\n--- Varianza Var(X) = p(1-p) ---\n")
cat("Var[H1N1]      =", round(var_h1n1, 4), "\n")
cat("Var[Estacional]=", round(var_seasonal, 4), "\n")

# Momentos: en una Bernoulli, E[X^k] = p para cualquier k >= 1
# (porque 0^k=0 y 1^k=1). Lo comprobamos empíricamente con k=2,3,4.
cat("\n--- Comprobación: E[X^k] = p para cualquier k (propiedad Bernoulli) ---\n")
for (k in 2:4) {
  cat("k =", k, "| E[H1N1^k] =", round(mean(datos$h1n1_vaccine^k), 4),
      "| p =", round(p_h1n1, 4), "\n")
}


# ============================================================
# Independencia: comprobación directa de la definición
# complementa el test chi-cuadrado ya existente
# ============================================================

p_conjunta_11 <- prop_conjunta["1", "1"]     # P(X=1, Y=1) observada
p_producto_11 <- p_h1n1 * p_seasonal         # P(X=1)*P(Y=1) esperada bajo independencia

cat("\n--- Comprobación directa de independencia ---\n")
cat("P(H1N1=1, Estacional=1) observada :", round(p_conjunta_11, 4), "\n")
cat("P(H1N1=1)*P(Estacional=1) esperada:", round(p_producto_11, 4), "\n")
cat("Diferencia:", round(p_conjunta_11 - p_producto_11, 4),
    "-> si fuese ~0 habría independencia; aquí es sustancial.\n")


# ============================================================
# Intervalos de confianza para p:
# Aproximación Normal (TCL) vs Binomial exacta
# ============================================================

z <- qnorm(0.975)

# Aproximación norma
se_h1n1 <- sqrt(p_h1n1 * (1 - p_h1n1) / n)
se_seasonal <- sqrt(p_seasonal * (1 - p_seasonal) / n)

ic_normal_h1n1 <- p_h1n1 + c(-1, 1) * z * se_h1n1
ic_normal_seasonal <- p_seasonal + c(-1, 1) * z * se_seasonal

cat("\n--- IC 95% para p (aproximación Normal a la Binomial) ---\n")
cat("H1N1:      [", round(ic_normal_h1n1[1], 4), ",", round(ic_normal_h1n1[2], 4), "]\n")
cat("Estacional:[", round(ic_normal_seasonal[1], 4), ",", round(ic_normal_seasonal[2], 4), "]\n")

#  Binomial exacta
ic_exacto_h1n1 <- binom.test(sum(datos$h1n1_vaccine), n)$conf.int
ic_exacto_seasonal <- binom.test(sum(datos$seasonal_vaccine), n)$conf.int

cat(" IC 95% para p (Binomial exacta, Clopper-Pearson) ")
cat("H1N1:      [", round(ic_exacto_h1n1[1], 4), ",", round(ic_exacto_h1n1[2], 4), "]\n")
cat("Estacional:[", round(ic_exacto_seasonal[1], 4), ",", round(ic_exacto_seasonal[2], 4), "]\n")


# ============================================================
# Estimación bayesiana de p con distribución Beta
# Beta es el conjugado natural de la Bernoulli/Binomial
# ============================================================

# Prior no informativo: Beta(1,1) = Uniforme(0,1)
alpha_prior <- 1
beta_prior <- 1

# Posterior tras observar los datos: Beta(alpha + éxitos, beta + fracasos)
post_h1n1 <- c(alpha_prior + sum(datos$h1n1_vaccine),
               beta_prior + n - sum(datos$h1n1_vaccine))
post_seasonal <- c(alpha_prior + sum(datos$seasonal_vaccine),
                   beta_prior + n - sum(datos$seasonal_vaccine))

# Media posterior y IC al 95%
media_post_h1n1 <- post_h1n1[1] / sum(post_h1n1)
ic_bayes_h1n1 <- qbeta(c(0.025, 0.975), post_h1n1[1], post_h1n1[2])

media_post_seasonal <- post_seasonal[1] / sum(post_seasonal)
ic_bayes_seasonal <- qbeta(c(0.025, 0.975), post_seasonal[1], post_seasonal[2])

cat("\n--- Estimación bayesiana de p: posterior Beta(a,b) ---\n")
cat("H1N1:       media posterior =", round(media_post_h1n1, 4),
    "| IC creíble 95% = [", round(ic_bayes_h1n1[1], 4), ",", round(ic_bayes_h1n1[2], 4), "]\n")
cat("Estacional: media posterior =", round(media_post_seasonal, 4),
    "| IC creíble 95% = [", round(ic_bayes_seasonal[1], 4), ",", round(ic_bayes_seasonal[2], 4), "]\n")

# Visualización de ambas posteriores
curva_beta <- data.frame(
  p = seq(0, 1, length.out = 1000)
) %>%
  mutate(
    Densidad_H1N1 = dbeta(p, post_h1n1[1], post_h1n1[2]),
    Densidad_Estacional = dbeta(p, post_seasonal[1], post_seasonal[2])
  )

ggplot(curva_beta, aes(x = p)) +
  geom_line(aes(y = Densidad_H1N1, color = "H1N1")) +
  geom_line(aes(y = Densidad_Estacional, color = "Estacional")) +
  labs(title = "Distribución posterior Beta de p (proporción vacunada)",
       x = "p", y = "Densidad", color = "Vacuna") +
  theme_minimal()


# ============================================================
# Cambio de variable: Z = X + Y
#  nº total de vacunas recibidas por persona (0, 1 o 2)
# ============================================================

datos$total_vacunas <- datos$h1n1_vaccine + datos$seasonal_vaccine

dist_Z <- prop.table(table(Z = datos$total_vacunas))
cat("\n--- Distribución de Z = X + Y (nº de vacunas recibidas) ---\n")
print(round(dist_Z, 4))

# Verificación teórica: E[Z] = E[X] + E[Y]  (siempre se cumple, haya o no independencia)
E_Z_teorico <- p_h1n1 + p_seasonal
E_Z_empirico <- mean(datos$total_vacunas)

# Verificación teórica: Var[Z] = Var(X) + Var(Y) + 2*Cov(X,Y)
Var_Z_teorico <- var_h1n1 + var_seasonal + 2 * covarianza
Var_Z_empirico <- var(datos$total_vacunas)

cat("\nE[Z]  teórico =", round(E_Z_teorico, 4), " | empírico =", round(E_Z_empirico, 4), "\n")
cat("Var[Z] teórico =", round(Var_Z_teorico, 4), " | empírico =", round(Var_Z_empirico, 4), "\n")

ggplot(as.data.frame(dist_Z), aes(x = Z, y = Freq)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = scales::percent(Freq, accuracy = 0.1)), vjust = -0.5) +
  labs(title = "Distribución de Z = nº total de vacunas recibidas",
       x = "Z (0, 1 o 2 vacunas)", y = "Proporción") +
  theme_minimal()