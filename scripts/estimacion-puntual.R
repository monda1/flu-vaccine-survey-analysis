library(dplyr)

# =====================================================================
# TEMA 8-9: Estimación puntual por Máxima Verosimilitud (MLE)
# y evaluación de propiedades del estimador
# (Insesgadez, ECM, Consistencia, Eficiencia/Cramér-Rao, Suficiencia,
#  Rao-Blackwell y Método de los Momentos)
# =====================================================================

# --------------------------------------------------------------
# 1. Planteamiento teórico 
# --------------------------------------------------------------
# Sea X_i ~ Bernoulli(p), i = 1,...,n, indicando si el individuo i
# se vacunó (1) o no (0).
#
# Verosimilitud:  L(p) = prod_{i=1}^n p^{x_i} (1-p)^{1-x_i}
# Log-verosimilitud: l(p) = (sum x_i) log(p) + (n - sum x_i) log(1-p)
# Derivando e igualando a 0:
#   dl/dp = (sum x_i)/p - (n - sum x_i)/(1-p) = 0
#   => p_hat_MLE = (1/n) sum x_i = x_barra (proporción muestral)
#
# Propiedades teóricas conocidas del estimador p_hat = X_barra:
#   E[p_hat] = p                  -> INSESGADO
#   Var(p_hat) = p(1-p)/n
#   ECM(p_hat) = Var(p_hat) + Sesgo^2 = p(1-p)/n + 0 = p(1-p)/n

# --------------------------------------------------------------
# 2. Cálculo de la MLE sobre los datos
# --------------------------------------------------------------
calcular_mle <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  p_hat <- mean(x)
  var_p_hat <- p_hat * (1 - p_hat) / n
  ee_p_hat <- sqrt(var_p_hat)
  ic_95 <- p_hat + c(-1, 1) * qnorm(0.975) * ee_p_hat
  
  list(n = n, p_hat = p_hat, var = var_p_hat, ee = ee_p_hat, ic_95 = ic_95)
}

mle_h1n1 <- calcular_mle(datos$h1n1_vaccine)
mle_seasonal <- calcular_mle(datos$seasonal_vaccine)

cat("=== MLE proporción vacunados H1N1 ===\n")
cat("n =", mle_h1n1$n, "\n")
cat("p_hat (MLE) =", round(mle_h1n1$p_hat, 4), "\n")
cat("Var(p_hat) =", round(mle_h1n1$var, 6), "\n")
cat("Error estándar =", round(mle_h1n1$ee, 4), "\n")
cat("IC 95% =", round(mle_h1n1$ic_95, 4), "\n\n")

cat("=== MLE proporción vacunados Estacional ===\n")
cat("n =", mle_seasonal$n, "\n")
cat("p_hat (MLE) =", round(mle_seasonal$p_hat, 4), "\n")
cat("Var(p_hat) =", round(mle_seasonal$var, 6), "\n")
cat("Error estándar =", round(mle_seasonal$ee, 4), "\n")
cat("IC 95% =", round(mle_seasonal$ic_95, 4), "\n\n")


# --------------------------------------------------------------
# 3.  Método de los Momentos (MM,  Máxima Verosimilitud
# --------------------------------------------------------------
# Para Bernoulli, el primer momento poblacional es mu_1 = E[X] = p.
# El método de momentos iguala el momento muestral al poblacional:
#   (1/n) sum x_i = p  =>  p_hat_MM = x_barra
# Es decir, para esta familia MM y MLE coinciden EXACTAMENTE
# (esto no siempre ocurre en otras familias de distribuciones).

p_mm_h1n1 <- mean(datos$h1n1_vaccine)
p_mm_seasonal <- mean(datos$seasonal_vaccine)

cat("--- Método de los Momentos vs MLE ---\n")
cat("H1N1:       p_hat_MM =", round(p_mm_h1n1, 4),
    "| p_hat_MLE =", round(mle_h1n1$p_hat, 4),
    "| Coinciden:", isTRUE(all.equal(p_mm_h1n1, mle_h1n1$p_hat)), "\n")
cat("Estacional: p_hat_MM =", round(p_mm_seasonal, 4),
    "| p_hat_MLE =", round(mle_seasonal$p_hat, 4),
    "| Coinciden:", isTRUE(all.equal(p_mm_seasonal, mle_seasonal$p_hat)), "\n\n")


# --------------------------------------------------------------
# 4. Comrpobacin de la MLE por optimización directa
#    (en vez de usar la fórmula cerrada, maximizamos l(p) numéricamente
#     para comprobar que coincide con x_barra)
# --------------------------------------------------------------
log_verosimilitud <- function(p, x) {
  sum(dbinom(x, size = 1, prob = p, log = TRUE))
}

x_h1n1 <- datos$h1n1_vaccine[!is.na(datos$h1n1_vaccine)]

optim_result <- optimize(
  f = log_verosimilitud,
  interval = c(1e-6, 1 - 1e-6),
  x = x_h1n1,
  maximum = TRUE
)


cat("p_hat (numérico) =", round(optim_result$maximum, 4), "\n")
cat("p_hat (fórmula cerrada) =", round(mean(x_h1n1), 4), "\n")
cat("Coinciden:", isTRUE(all.equal(optim_result$maximum, mean(x_h1n1), tolerance = 1e-4)), "\n\n")

# Gráfico de la log-verosimilitud
curve(sapply(x, log_verosimilitud, x = x_h1n1),
      from = 0.01, to = 0.99, n = 300,
      xlab = "p", ylab = "log-verosimilitud l(p)",
      main = "Log-verosimilitud para p (H1N1)")
abline(v = mean(x_h1n1), col = "red", lty = 2)
legend("bottomright", legend = paste("p_hat =", round(mean(x_h1n1), 4)),
       col = "red", lty = 2, bty = "n")


# --------------------------------------------------------------
# 5. Suficiencia del estadístico T = sum(X_i)
# --------------------------------------------------------------
# Teorema de factorización: L(p) depende de los datos SOLO a través
# de T = sum(x_i), no del orden ni de qué observaciones concretas
# fueron 0 o 1. Lo comprobamos empíricamente: barajamos los datos
# (la suma no cambia) y verificamos que la log-verosimilitud es idéntica.

set.seed(42)
x_barajado <- sample(x_h1n1)  # mismo T = sum(x), distinto orden

ll_original <- log_verosimilitud(mle_h1n1$p_hat, x_h1n1)
ll_barajado <- log_verosimilitud(mle_h1n1$p_hat, x_barajado)

cat(" Comprobación de suficiencia de T = sum(X_i) ")
cat("Suma original  T =", sum(x_h1n1), " | log-verosim. =", round(ll_original, 4), "\n")
cat("Suma barajada  T =", sum(x_barajado), " | log-verosim. =", round(ll_barajado, 4), "\n")
cat("Idénticas (misma T, distinto orden):",
    isTRUE(all.equal(ll_original, ll_barajado)), "\n")
cat("-> Confirma que T = sum(X_i) concentra TODA la información sobre p.\n\n")


# --------------------------------------------------------------
# 6.  Información de Fisher y Cota de Cramér-Rao (Eficiencia)
# --------------------------------------------------------------
# Para una Bernoulli, la información de Fisher de una muestra de
# tamaño n es:  I_n(p) = n / [p(1-p)]
# La Cota de Cramér-Rao (varianza mínima posible para un estimador
# insesgado) es:  CRLB = 1 / I_n(p) = p(1-p) / n
#
# Si Var(p_hat) == CRLB, el estimador es EFICIENTE (óptimo posible).

fisher_info <- function(p_hat, n) n / (p_hat * (1 - p_hat))

I_h1n1 <- fisher_info(mle_h1n1$p_hat, mle_h1n1$n)
CRLB_h1n1 <- 1 / I_h1n1

I_seasonal <- fisher_info(mle_seasonal$p_hat, mle_seasonal$n)
CRLB_seasonal <- 1 / I_seasonal

cat("--- Información de Fisher y Cota de Cramér-Rao ---\n")
cat("H1N1:       I_n(p) =", round(I_h1n1, 2),
    "| CRLB =", round(CRLB_h1n1, 6),
    "| Var(p_hat) =", round(mle_h1n1$var, 6),
    "| Eficiente:", isTRUE(all.equal(CRLB_h1n1, mle_h1n1$var)), "\n")
cat("Estacional: I_n(p) =", round(I_seasonal, 2),
    "| CRLB =", round(CRLB_seasonal, 6),
    "| Var(p_hat) =", round(mle_seasonal$var, 6),
    "| Eficiente:", isTRUE(all.equal(CRLB_seasonal, mle_seasonal$var)), "\n")
cat("-> El MLE ALCANZA la cota de Cramér-Rao: es el estimador insesgado\n")
cat("   de MENOR varianza posible (estimador eficiente).\n\n")


# --------------------------------------------------------------
# 7. Evaluación empírica de insesgadez y ECM vía simulación Monte Carlo
#    Tomamos como "p verdadero" el valor estimado en la muestra completa
#    y simulamos repetidamente muestras de tamaño n para verificar que
#    el estimador se comporta como predice la teoría.
# --------------------------------------------------------------
set.seed(123)

evaluar_estimador_mc <- function(p_verdadero, n, n_simulaciones = 5000) {
  p_hats <- replicate(n_simulaciones, {
    muestra <- rbinom(n, size = 1, prob = p_verdadero)
    mean(muestra)
  })
  
  sesgo <- mean(p_hats) - p_verdadero
  varianza_empirica <- var(p_hats)
  ecm_empirico <- mean((p_hats - p_verdadero)^2)
  
  varianza_teorica <- p_verdadero * (1 - p_verdadero) / n
  ecm_teorico <- varianza_teorica  # porque el sesgo teórico es 0
  
  list(
    p_hats = p_hats,
    sesgo_empirico = sesgo,
    varianza_empirica = varianza_empirica,
    ecm_empirico = ecm_empirico,
    varianza_teorica = varianza_teorica,
    ecm_teorico = ecm_teorico
  )
}

resultado_mc <- evaluar_estimador_mc(p_verdadero = mle_h1n1$p_hat, n = mle_h1n1$n)

cat("=== Evaluación Monte Carlo del estimador (H1N1) ===\n")
cat("Sesgo empírico:", round(resultado_mc$sesgo_empirico, 6), "(teórico: 0)\n")
cat("Varianza empírica:", round(resultado_mc$varianza_empirica, 6),
    "| Varianza teórica (= CRLB):", round(resultado_mc$varianza_teorica, 6), "\n")
cat("ECM empírico:", round(resultado_mc$ecm_empirico, 6),
    "| ECM teórico:", round(resultado_mc$ecm_teorico, 6), "\n\n")

# Histograma de la distribución muestral del estimador
hist(resultado_mc$p_hats, breaks = 50, probability = TRUE,
     main = "Distribución muestral de p_hat (simulación Monte Carlo)",
     xlab = expression(hat(p)), col = "lightblue", border = "white")
abline(v = mle_h1n1$p_hat, col = "red", lwd = 2, lty = 2)
curve(dnorm(x, mean = mle_h1n1$p_hat, sd = sqrt(resultado_mc$varianza_teorica)),
      add = TRUE, col = "darkgreen", lwd = 2)
legend("topright",
       legend = c("p verdadero", "Normal asintótica teórica"),
       col = c("red", "darkgreen"), lty = c(2, 1), lwd = 2, bty = "n")


# --------------------------------------------------------------
# 8.  Consistencia del estimador (convergencia cuando n -> infinito)
# --------------------------------------------------------------
# Un estimador es consistente si, al crecer n, converge en probabilidad
# al valor verdadero p. Lo ilustramos simulando muestras de tamaños
# crecientes y observando cómo se estrecha la distribución de p_hat
# alrededor de p (Ley de los Grandes Números).

p_verdadero <- mle_h1n1$p_hat
n_seq <- c(10, 50, 100, 500, 1000, 5000, mle_h1n1$n)

set.seed(7)
resultados_consistencia <- lapply(n_seq, function(n_i) {
  p_hats_ni <- replicate(2000, mean(rbinom(n_i, 1, p_verdadero)))
  data.frame(n = n_i, p_hat = p_hats_ni)
})
df_consistencia <- bind_rows(resultados_consistencia)
df_consistencia$n <- factor(df_consistencia$n, levels = n_seq)

cat("--- Consistencia: varianza de p_hat según n (debe decrecer como 1/n) ---\n")
tabla_consistencia <- df_consistencia %>%
  group_by(n) %>%
  summarise(media = mean(p_hat), varianza = var(p_hat), .groups = "drop")
print(tabla_consistencia)
cat("-> A medida que n crece, la varianza se reduce y p_hat se concentra\n")
cat("   cada vez más en torno a p =", round(p_verdadero, 4), ": esto ES la consistencia.\n\n")

boxplot(p_hat ~ n, data = df_consistencia,
        main = "Consistencia de p_hat: convergencia al crecer n",
        xlab = "Tamaño muestral (n)", ylab = expression(hat(p)),
        col = "lightblue")
abline(h = p_verdadero, col = "red", lty = 2, lwd = 2)
legend("topright", legend = "p verdadero", col = "red", lty = 2, bty = "n")


# --------------------------------------------------------------
# 9. Teorema de Rao-Blackwell
# --------------------------------------------------------------
# Partimos de un estimador insesgado pero "ingenuo": delta = X_1
# (usar solo la primera observación). Es insesgado, E[X_1] = p,
# pero muy ineficiente: Var(X_1) = p(1-p), MUCHO mayor que p(1-p)/n.
#
# El Teorema de Rao-Blackwell dice que, condicionando en el
# estadístico suficiente T = sum(X_i), obtenemos un estimador
# insesgado de varianza igual o menor:
#   delta_RB = E[X_1 | T = t] = t / n = p_hat
# (resultado clásico: dado T=t, X_1 se comporta como una extracción
#  "aleatoria" entre los t unos y (n-t) ceros, de ahí que su esperanza
#  condicional sea t/n). Es decir: ¡Rao-Blackwell reconstruye la MLE!

delta_ingenuo <- x_h1n1[1]          # estimador crudo: solo la 1ª observación
var_delta_ingenuo <- p_verdadero * (1 - p_verdadero)   # Var(X_1) teórica
delta_rao_blackwell <- sum(x_h1n1) / length(x_h1n1)    # = p_hat

cat("--- Teorema de Rao-Blackwell ---\n")
cat("Estimador ingenuo delta = X_1 =", delta_ingenuo,
    "| Var(delta) teórica =", round(var_delta_ingenuo, 6), "\n")
cat("Rao-Blackwell: E[X_1 | T=t] = t/n =", round(delta_rao_blackwell, 4),
    "| Var(p_hat) =", round(mle_h1n1$var, 6), "\n")
cat("Reducción de varianza: de", round(var_delta_ingenuo, 6), "a",
    round(mle_h1n1$var, 6), "( factor ~", round(var_delta_ingenuo / mle_h1n1$var, 1), "veces menor )\n")
cat("-> El estimador Rao-Blackwellizado COINCIDE con la MLE: esto (junto con\n")
cat("   la suficiencia completa de T) es la prueba de que p_hat es el UMVUE\n")
cat("   (Teorema de Lehmann-Scheffé).\n\n")


# --------------------------------------------------------------
# 10. Tabla resumen final
# --------------------------------------------------------------
resumen <- data.frame(
  Vacuna = c("H1N1", "Estacional"),
  n = c(mle_h1n1$n, mle_seasonal$n),
  p_hat_MLE = c(mle_h1n1$p_hat, mle_seasonal$p_hat),
  p_hat_MM = c(p_mm_h1n1, p_mm_seasonal),
  Var_teorica = c(mle_h1n1$var, mle_seasonal$var),
  CRLB = c(CRLB_h1n1, CRLB_seasonal),
  Es_eficiente = c(isTRUE(all.equal(CRLB_h1n1, mle_h1n1$var)),
                   isTRUE(all.equal(CRLB_seasonal, mle_seasonal$var))),
  EE = c(mle_h1n1$ee, mle_seasonal$ee),
  IC95_inf = c(mle_h1n1$ic_95[1], mle_seasonal$ic_95[1]),
  IC95_sup = c(mle_h1n1$ic_95[2], mle_seasonal$ic_95[2])
)
print(resumen)