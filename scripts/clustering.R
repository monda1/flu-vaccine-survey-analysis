library(conflicted)
library(dplyr)
library(klaR)      
library(ggplot2)
library(tidyr)


conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")

# =====================================================================
# Análisis de conglomerados con k-modes sobre variables
# categóricas/actitudinales
# =====================================================================

# --------------------------------------------------------------
# 1. Selección de variables para el clustering
#    Nos centramos en variables actitudinales (opinión, preocupación,
#    conocimiento, comportamiento), No demográficas, precisamente para
#    poder comprobar después si los clusters resultantes se alinean
#    o no con la demografía (esa es la pregunta de interés)
# --------------------------------------------------------------
vars_cluster <- c("h1n1_concern", "h1n1_knowledge",
                  "opinion_h1n1_vacc_effective", "opinion_h1n1_risk", "opinion_h1n1_sick_from_vacc",
                  "opinion_seas_vacc_effective", "opinion_seas_risk", "opinion_seas_sick_from_vacc",
                  "behavioral_antiviral_meds", "behavioral_avoidance", "behavioral_face_mask",
                  "behavioral_wash_hands", "behavioral_large_gatherings",
                  "behavioral_outside_home", "behavioral_touch_face",
                  "doctor_recc_h1n1", "doctor_recc_seasonal")

datos_cluster <- datos %>%
  dplyr::select(all_of(vars_cluster)) %>%
  na.omit() %>%
  mutate(across(everything(), as.factor)) %>%
  as.data.frame()   # <- FIX: kmodes() no funciona bien con tibbles

cat("Observaciones para el clustering:", nrow(datos_cluster),
    "de", nrow(datos), "\n\n")

# --------------------------------------------------------------
# 2. Selección del número óptimo de clusters (k)
#    k-modes no tiene un equivalente directo al codo de inercia de
#    k-means, pero el coste (suma de disimilitudes intra-cluster)
#    sirve de forma análoga
# --------------------------------------------------------------
set.seed(123)

# --------------------------------------------------------------
# múltiples arranques aleatorios por valor de k
# --------------------------------------------------------------
# La función original recibía un parámetro n_starts pero nunca lo usaba
# dentro del cuerpo: solo se ejecutaba kmodes() una vez por k. Esto es
# un problema real: k-modes (igual que k-means) es
# sensible al punto de partida aleatorio y puede converger a un óptimo
# local distinto cada vez. Sin varios arranques, el "codo" que se
# dibuja depende en parte de qué arranque tuvo suerte, no solo de k.
# Ahora sí se repite n_starts veces por cada k y nos quedamos con el
# mejor (menor coste intra-cluster).

evaluar_k <- function(k, datos, n_starts = 5) {
  costes_intentos <- sapply(1:n_starts, function(i) {
    modelo <- kmodes(datos, modes = k, iter.max = 20, weighted = FALSE, fast = TRUE)
    sum(modelo$withindiff)
  })
  min(costes_intentos)
}

k_candidatos <- 2:8
costes <- sapply(k_candidatos, function(k) evaluar_k(k, datos_cluster))

df_codo <- data.frame(k = k_candidatos, coste = costes)

ggplot(df_codo, aes(x = k, y = coste)) +
  geom_line(color = "steelblue") +
  geom_point(size = 3, color = "steelblue") +
  scale_x_continuous(breaks = k_candidatos) +
  labs(title = "Método del codo para k-modes (mejor de 5 arranques por k)",
       x = "Número de clusters (k)", y = "Coste (disimilitud intra-cluster)") +
  theme_minimal()

print(df_codo)

# --------------------------------------------------------------
# Validación con silueta media
# --------------------------------------------------------------
# La reducción decrece de
# forma monótona y suave" (sin codo marcado)
# El coeficiente de silueta mide, para cada punto, cuánto más parecido
# es a su propio cluster que al cluster vecino más cercano
# (rango [-1, 1]; valores más altos = clusters mejor separados).
# Para poder calcularlo con variables categóricas se usa la distancia
# de Gower (que generaliza la disimilitud de coincidencia simple).
#
# calcular una matriz de distancias completa con
# ~23.000 filas implicaría del orden de 270 millones de pares, muy
# costoso en memoria. Se usa una m.a  para este
# diagnóstico suficiente para comparar k's entre sí, no para
# reasignar el clustering final.

if (requireNamespace("cluster", quietly = TRUE)) {
  library(cluster)
  
  set.seed(456)
  muestra_idx <- sample(1:nrow(datos_cluster), size = min(2000, nrow(datos_cluster)))
  muestra_silueta <- datos_cluster[muestra_idx, ]
  dist_gower_muestra <- daisy(muestra_silueta, metric = "gower")
  
  evaluar_silueta_k <- function(k) {
    modelo <- kmodes(muestra_silueta, modes = k, iter.max = 20, weighted = FALSE)
    sil <- silhouette(modelo$cluster, dist_gower_muestra)
    mean(sil[, "sil_width"])
  }
  
  siluetas_k <- sapply(k_candidatos, evaluar_silueta_k)
  df_silueta <- data.frame(k = k_candidatos, silueta_media = siluetas_k)
  
  cat("\n=== Silueta media por k (sobre muestra de", nrow(muestra_silueta), "obs.) ===\n")
  print(df_silueta)
  cat("Regla orientativa: >0.5 estructura razonable, 0.25-0.5 débil pero real,",
      "<0.25 estructura poco clara (habitual con variables actitudinales/Likert,",
      "donde no se espera una separación tan nítida como en datos numéricos continuos).\n\n")
  
  ggplot(df_silueta, aes(x = k, y = silueta_media)) +
    geom_line(color = "darkorange") +
    geom_point(size = 3, color = "darkorange") +
    scale_x_continuous(breaks = k_candidatos) +
    labs(title = "Silueta media por número de clusters",
         subtitle = "Complementa al codo: mide separación real, no solo coste decreciente",
         x = "k", y = "Silueta media (Gower, sobre muestra)") +
    theme_minimal()
}

# --------------------------------------------------------------
# 3. Ajuste del modelo final con el k elegido
# --------------------------------------------------------------
k_elegido <- 4

modelo_kmodes <- kmodes(datos_cluster, modes = k_elegido, iter.max = 50, weighted = FALSE)

cat("=== Tamaño de cada cluster ===\n")
print(table(modelo_kmodes$cluster))

cat("\n=== Modas de cada cluster (perfil central) ===\n")
print(modelo_kmodes$modes)

# el número asignado a cada cluster (1, 2, 3, 4) es arbitrario y puede
# cambiar de una ejecución a otra aunque la ESTRUCTURA encontrada sea
# la misma . Cualquier informe que diga "el
# Conglomerado 3 es el más vacunado" debe verificar esa correspondencia
# cada vez que se re-ejecute el modelo, no asumir que el número persiste.

# --------------------------------------------------------------
# 4. Añadimos la asignación de cluster a los datos originales
#    (usando el mismo subconjunto de filas sin NA que se usó en el
#    clustering, para poder cruzarlo después con vacunación y demografía)
# --------------------------------------------------------------
indices_completos <- datos %>%
  dplyr::select(all_of(vars_cluster)) %>%
  complete.cases() %>%
  which()

datos_con_cluster <- datos[indices_completos, ] %>%
  mutate(cluster = factor(modelo_kmodes$cluster))

# --------------------------------------------------------------
# 5. Perfilado de los clusters: distribución de cada variable
#    actitudinal dentro de cada cluster
# --------------------------------------------------------------
perfil_clusters <- datos_con_cluster %>%
  dplyr::select(cluster, all_of(vars_cluster)) %>%
  pivot_longer(-cluster, names_to = "variable", values_to = "valor") %>%
  group_by(cluster, variable, valor) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster, variable) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Visualización: perfil de un subconjunto de variables clave por cluster
vars_clave <- c("h1n1_concern", "opinion_h1n1_risk", "opinion_seas_risk",
                "behavioral_wash_hands", "doctor_recc_h1n1")

ggplot(perfil_clusters %>% filter(variable %in% vars_clave),
       aes(x = valor, y = prop, fill = cluster)) +
  geom_col(position = "dodge") +
  facet_wrap(~variable, scales = "free_x") +
  labs(title = "Perfil actitudinal de cada cluster (variables seleccionadas)",
       x = "Valor de la variable", y = "Proporción dentro del cluster", fill = "Cluster") +
  theme_minimal()

# --------------------------------------------------------------
# 6. ¿los clusters se alinean con la tasa de
#    vacunación observada?
# --------------------------------------------------------------
tasa_vacunacion_cluster <- datos_con_cluster %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    tasa_h1n1 = mean(h1n1_vaccine, na.rm = TRUE),
    tasa_seasonal = mean(seasonal_vaccine, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== Tasa de vacunación por cluster actitudinal ===\n")
print(tasa_vacunacion_cluster)

ggplot(tasa_vacunacion_cluster %>%
         pivot_longer(cols = c(tasa_h1n1, tasa_seasonal),
                      names_to = "vacuna", values_to = "tasa"),
       aes(x = cluster, y = tasa, fill = vacuna)) +
  geom_col(position = "dodge") +
  labs(title = "Tasa de vacunación según cluster actitudinal",
       x = "Cluster", y = "Proporción vacunada", fill = "Vacuna") +
  theme_minimal()

# --------------------------------------------------------------
#  ¿los clusters se corresponden con
#    categorías demográficas obvias, o cortan "en diagonal"?
#    Cruzamos cluster con variables demográficas mediante chi-cuadrado
#    y V de Cramér, y visualmente con tablas de proporciones
# --------------------------------------------------------------
vars_demograficas <- c("age_group", "education", "income_poverty",
                       "sex", "race", "employment_status")

evaluar_asociacion_cluster <- function(datos, var_demo) {
  df <- datos %>% filter(!is.na(.data[[var_demo]]))
  tabla <- table(df$cluster, df[[var_demo]])
  test <- chisq.test(tabla)
  n <- sum(tabla)
  k_dim <- min(dim(tabla))
  v_cramer <- sqrt(test$statistic / (n * (k_dim - 1)))
  data.frame(variable = var_demo,
             chi2 = round(test$statistic, 2),
             p_valor = format.pval(test$p.value, digits = 4),
             v_cramer = round(v_cramer, 4))
}

resumen_demografia <- bind_rows(
  lapply(vars_demograficas, function(v) evaluar_asociacion_cluster(datos_con_cluster, v))
)

cat("\n=== ¿Los clusters se explican por variables demográficas? ===\n")
cat("(V de Cramér bajo pese a p-valor significativo -> los clusters NO son",
    "un simple proxy de esa variable demográfica; con N > 20.000 hasta",
    "diferencias marginales resultan 'significativas', por eso el tamaño",
    "del efecto manda aquí más que el p-valor)\n\n")
print(resumen_demografia %>% arrange(desc(v_cramer)))

# Visualización de la composición demográfica de cada cluster
ggplot(datos_con_cluster %>% filter(!is.na(age_group)),
       aes(x = cluster, fill = age_group)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Composición por grupo de edad dentro de cada cluster",
       subtitle = "Si los clusters fueran solo un proxy de la edad, cada barra sería casi monocolor",
       x = "Cluster", y = "Proporción", fill = "Grupo de edad") +
  theme_minimal()

ggplot(datos_con_cluster %>% filter(!is.na(education)),
       aes(x = cluster, fill = education)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Composición por nivel educativo dentro de cada cluster",
       x = "Cluster", y = "Proporción", fill = "Educación") +
  theme_minimal()

# =====================================================================
# 8. ¿Cambia algo si respetamos el ORDEN de las variables Likert?
#    k-modes con distancia de coincidencia simple trata "concern=1 vs 2"
#    exactamente igual que "concern=1 vs 3": para k-modes, todo lo que
#    no coincide es un desacuerdo binario, sin noción de "más cerca" o
#    "más lejos". Eso tira información real (el orden de la escala).
#    PAM (Partitioning Around Medoids) sobre una distancia de Gower que
#    SÍ declara las variables Likert como ordinales ("ordratio") es una
#    alternativa más fiel a la naturaleza de los datos.
# =====================================================================
if (requireNamespace("cluster", quietly = TRUE)) {
  
  vars_ordinales <- c("h1n1_concern", "h1n1_knowledge",
                      "opinion_h1n1_vacc_effective", "opinion_h1n1_risk", "opinion_h1n1_sick_from_vacc",
                      "opinion_seas_vacc_effective", "opinion_seas_risk", "opinion_seas_sick_from_vacc")
  
  datos_gower <- datos_cluster
  datos_gower[vars_ordinales] <- lapply(datos_gower[vars_ordinales],
                                        function(x) factor(x, ordered = TRUE))
  
  # Misma muestra que en la validación por silueta, por coherencia y
  # coste computacional (Gower en las ~23.000 filas completas es
  # posible pero lento; para un chequeo de robustez metodológica una
  # muestra representativa es suficiente).
  muestra_gower <- datos_gower[muestra_idx, ]
  
  dist_gower_ordinal <- daisy(muestra_gower, metric = "gower",
                              type = list(ordratio = vars_ordinales))
  
  modelo_pam <- pam(dist_gower_ordinal, k = k_elegido, diss = TRUE)
  
  cat("\n=== Comparación k-modes vs. PAM+Gower (misma muestra, k =", k_elegido, ") ===\n")
  cat("Silueta media k-modes (calculada antes):",
      round(mean(silhouette(kmodes(muestra_silueta, modes = k_elegido)$cluster,
                            dist_gower_muestra)[, "sil_width"]), 4), "\n")
  cat("Silueta media PAM+Gower (respetando el orden Likert):",
      round(mean(silhouette(modelo_pam)[, "sil_width"]), 4), "\n")

}