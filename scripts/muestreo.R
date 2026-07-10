library(conflicted)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")

# =====================================================================
# IV. MUESTREO
# Población objetivo vs muestreada, marco de muestreo,
#         sesgos de selección y de medición
#
# La NHFS es una encuesta telefónica RDD de doble marco (fijo + móvil),
# dirigida a todas las personas de EE. UU., realizada oct-2009/jun-2010.
# Aquí NO tenemos acceso al censo real ni a los pesos muestrales
# (no viene columna de "weight" en el fichero), así que el análisis se
# centra en lo que SÍ podemos comprobar con los datos disponibles:
#   1) Comparación con puntos de referencia poblacionales aproximados
#      (Census Bureau 2010) para detectar posible sesgo de cobertura
#   2) El sesgo de "complete case" que introducimos nosotros al
#      hacer na.omit() antes de clusterizar (Tema 26)
#   3) Si la falta de respuesta en variables actitudinales/clínicas
#      está relacionada con la demografía (missingness no aleatorio)
#   4) Indicios indirectos de sesgo de medición (autodeclaración)
#      comparando recomendación médica vs vacunación autodeclarada
#   5) Unidad de muestreo vs unidad de observación: el
#      mecanismo de selección DENTRO del hogar como fuente de sesgo
#      no corregida al no haber pesos muestrales
#   6) Taxonomía: clasificar cada hallazgo anterior como error
#      muestral o ajeno al muestreo, y por qué el muestreo seguía
#      siendo preferible a un censo en este contexto concreto
#   7) Estimador insesgado y medida de error: demostración por
#      simulación, y por qué el error "de manual" (SRS) probablemente
#      subestima el error real de este diseño
# =====================================================================

# --------------------------------------------------------------
# 0. Vector de variables usado en el clustering
# --------------------------------------------------------------
vars_cluster <- c("h1n1_concern", "h1n1_knowledge",
                  "opinion_h1n1_vacc_effective", "opinion_h1n1_risk", "opinion_h1n1_sick_from_vacc",
                  "opinion_seas_vacc_effective", "opinion_seas_risk", "opinion_seas_sick_from_vacc",
                  "behavioral_antiviral_meds", "behavioral_avoidance", "behavioral_face_mask",
                  "behavioral_wash_hands", "behavioral_large_gatherings",
                  "behavioral_outside_home", "behavioral_touch_face",
                  "doctor_recc_h1n1", "doctor_recc_seasonal")

vars_demograficas <- c("age_group", "education", "income_poverty",
                       "sex", "race", "employment_status")

# --------------------------------------------------------------
# 1. SESGO DE COBERTURA: comparación con puntos de referencia
#    poblacionales aproximados (Census Bureau, 2010)
# --------------------------------------------------------------
ref_edad_censo <- c(
  "18 - 34 Years" = 0.30,
  "35 - 44 Years" = 0.18,
  "45 - 54 Years" = 0.19,
  "55 - 64 Years" = 0.16,
  "65+ Years"     = 0.17
)

dist_muestra_edad <- datos %>%
  filter(!is.na(age_group)) %>%
  count(age_group) %>%
  mutate(prop_muestra = n / sum(n))


print(levels(factor(datos$age_group)))

cat("\n=== Comparación distribución muestral vs referencia censal aproximada ===\n")
print(dist_muestra_edad)

ggplot(dist_muestra_edad, aes(x = age_group, y = prop_muestra)) +
  geom_col(fill = "steelblue") +
  scale_y_continuous(labels = percent) +
  labs(title = "Distribución por edad en la muestra NHFS",
       subtitle = "Contrastar visualmente con la estructura de edad de EE. UU. (Census 2010)",
       x = "Grupo de edad", y = "Proporción en la muestra") +
  theme_minimal()

# --------------------------------------------------------------
# 2. SESGO 
# --------------------------------------------------------------
indices_completos <- datos %>%
  dplyr::select(all_of(vars_cluster)) %>%
  complete.cases() %>%
  which()

datos$en_complete_case <- FALSE
datos$en_complete_case[indices_completos] <- TRUE

cat("\n=== Filas retenidas tras na.omit() para el clustering ===\n")
cat(length(indices_completos), "de", nrow(datos),
    sprintf("(%.1f%%)\n", 100 * length(indices_completos) / nrow(datos)))

comparar_demografia_attrition <- function(var_demo) {
  datos %>%
    filter(!is.na(.data[[var_demo]])) %>%
    group_by(en_complete_case, .data[[var_demo]]) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    mutate(variable = var_demo) %>%
    rename(valor = 2)
}

comparacion_attrition <- bind_rows(lapply(vars_demograficas, comparar_demografia_attrition))

test_attrition <- function(var_demo) {
  df <- datos %>% filter(!is.na(.data[[var_demo]]))
  tabla <- table(df$en_complete_case, df[[var_demo]])
  test <- chisq.test(tabla)
  data.frame(variable = var_demo,
             chi2 = round(unname(test$statistic), 2),
             p_valor = format.pval(test$p.value, digits = 4))
}

cat("\n=== ¿La retención en el subconjunto complete-case depende de la demografía? ===\n")
resultado_attrition <- bind_rows(lapply(vars_demograficas, test_attrition))
print(resultado_attrition)

ggplot(comparacion_attrition %>% filter(variable == "income_poverty"),
       aes(x = valor, y = prop, fill = en_complete_case)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = percent) +
  labs(title = "Renta: muestra completa vs subconjunto complete-case",
       subtitle = "Si las barras difieren mucho, el na.omit() no es neutral respecto a la renta",
       x = "Nivel de renta/pobreza", y = "Proporción", fill = "¿En complete-case?") +
  theme_minimal()

# --------------------------------------------------------------
# 3. ¿LA FALTA DE RESPUESTA EN VARIABLES ACTITUDINALES ES ALEATORIA?
# --------------------------------------------------------------
tabla_missingness <- sapply(vars_cluster, function(v) mean(is.na(datos[[v]])))
tabla_missingness <- data.frame(variable = names(tabla_missingness),
                                pct_perdido = round(100 * tabla_missingness, 2)) %>%
  arrange(desc(pct_perdido))

cat("\n=== % de valores perdidos por variable actitudinal/clínica ===\n")
print(tabla_missingness)

test_missingness_demografia <- function(var_objetivo, var_demo) {
  df <- datos %>%
    filter(!is.na(.data[[var_demo]])) %>%
    mutate(perdido = is.na(.data[[var_objetivo]]))
  tabla <- table(df$perdido, df[[var_demo]])
  if (nrow(tabla) < 2) {
    return(data.frame(variable_objetivo = var_objetivo, variable_demo = var_demo,
                      chi2 = NA, p_valor = NA))
  }
  test <- suppressWarnings(chisq.test(tabla))
  data.frame(variable_objetivo = var_objetivo, variable_demo = var_demo,
             chi2 = round(unname(test$statistic), 2),
             p_valor = format.pval(test$p.value, digits = 4))
}

resumen_missingness <- bind_rows(
  lapply(vars_cluster, function(vo) {
    bind_rows(lapply(vars_demograficas, function(vd) test_missingness_demografia(vo, vd)))
  })
)

cat("\n=== Asociación entre missingness y demografía (p < 0.05 sugiere MAR/MNAR) ===\n")
print(resumen_missingness %>% arrange(p_valor) %>% head(15))

# --------------------------------------------------------------
# 4. INDICIO INDIRECTO DE SESGO DE MEDICIÓN
# --------------------------------------------------------------
coherencia_h1n1 <- datos %>%
  filter(!is.na(doctor_recc_h1n1), !is.na(h1n1_vaccine)) %>%
  group_by(doctor_recc_h1n1) %>%
  summarise(n = n(),
            tasa_vacunacion = mean(h1n1_vaccine),
            .groups = "drop")

cat("\n=== Coherencia recomendación médica vs vacunación H1N1 autodeclarada ===\n")
print(coherencia_h1n1)

coherencia_concern <- datos %>%
  filter(!is.na(h1n1_concern), !is.na(h1n1_vaccine)) %>%
  group_by(h1n1_concern) %>%
  summarise(n = n(),
            tasa_vacunacion = mean(h1n1_vaccine),
            .groups = "drop")

cat("\n=== Tasa de vacunación H1N1 autodeclarada según nivel de preocupación ===\n")
print(coherencia_concern)

ggplot(coherencia_h1n1, aes(x = factor(doctor_recc_h1n1), y = tasa_vacunacion)) +
  geom_col(fill = "darkorange") +
  scale_y_continuous(labels = percent) +
  labs(title = "Tasa de vacunación H1N1 autodeclarada según recomendación médica",
       x = "¿Recomendación médica de vacunarse?", y = "Tasa de vacunación autodeclarada") +
  theme_minimal()

# =====================================================================
# 5. Unidad de MUESTREO vs unidad de OBSERVACIÓN
# =====================================================================
# En un RDD como este, la unidad de MUESTREO es el número de teléfono
# (aproximadamente, el hogar), mientras que la unidad de OBSERVACIÓN
# que acaba en cada fila del dataset es el INDIVIDUO seleccionado
# dentro de ese hogar para responder. Esta distinción no es solo
# terminológica: si dentro de cada hogar se selecciona a UN adulto al
# azar (procedimiento habitual en este tipo de encuestas) SIN ponderar
# después por el tamaño del hogar, las personas que viven en hogares
# GRANDES quedan sistemáticamente infrarrepresentadas por persona,
# porque cada una de ellas "compite" con más convivientes por la
# única plaza de respondente de su hogar.
#
# No tenemos el marco real de hogares de EE. UU. en 2009 para probarlo
# directamente sobre estos datos, así que se ilustra el MECANISMO con
# una simulación mínima y autocontenida  no es una prueba de que
# exista sesgo aquí, es la demostración de por qué el mecanismo, de
# existir, produciría sesgo si no se pondera.

set.seed(789)
poblacion_hogares <- data.frame(
  hogar_id = 1:10000,
  tamano_hogar = sample(1:4, 10000, replace = TRUE, prob = c(0.30, 0.35, 0.20, 0.15))
)

# Universo real de personas (cada persona de cada hogar es un individuo)
poblacion_personas <- poblacion_hogares %>%
  tidyr::uncount(tamano_hogar) %>%
  mutate(persona_id = row_number())

# Muestreo tal como lo hace la encuesta: 1 hogar = 1 selección,
# después 1 persona al azar DENTRO de ese hogar, SIN ponderar por tamaño
hogares_muestreados <- poblacion_hogares[sample(1:nrow(poblacion_hogares), 2000), ]

cat("\n=== Simulación: sesgo por selección dentro del hogar sin ponderar ===\n")
cat("Distribución REAL de tamaños de hogar en la población simulada:\n")
print(round(prop.table(table(poblacion_hogares$tamano_hogar)), 3))

cat("\nDistribución de tamaños de hogar de las PERSONAS seleccionadas",
    "(1 por hogar muestreado, sin ponderar):\n")
print(round(prop.table(table(hogares_muestreados$tamano_hogar)), 3))

cat("\nDistribución de tamaños de hogar en la población de PERSONAS",
    "(lo que se querría representar, sumando sobre todos los individuos):\n")
print(round(prop.table(table(poblacion_personas$tamano_hogar)), 3))

cat("\nConclusión de la simulación: la muestra de personas seleccionadas",
    "así reproduce la distribución de tamaños de hogar (unidad de muestreo),",
    "NO la distribución real de personas por tamaño de hogar (unidad de",
    "observación de interés) -- las personas de hogares grandes quedan",
    "infrarrepresentadas por individuo salvo que se pondere por 1/tamaño_hogar.\n")

cat("\nEn el dataset real: sin columna de pesos muestrales, no podemos",
    "corregir este mecanismo aunque sí disponemos de household_adults y",
    "household_children, que podrían usarse para reconstruir un peso",
    "aproximado (1 / (household_adults + household_children + 1)) si se",
    "quisiera aproximar una corrección post-hoc.\n")

# --------------------------------------------------------------
# 5bis. Comprobación descriptiva en los datos reales (sin benchmark)
# --------------------------------------------------------------
# Esto es puramente descriptivo: muestra CÓMO es la distribución
# muestral de adultos por hogar, pero sin un marco poblacional real
# de EE. UU. para 2009 no podemos concluir si hay o no sesgo -- solo
# dejar constancia de la limitación, como en el punto 1 con la edad.
dist_hogar_muestra <- datos %>%
  filter(!is.na(household_adults)) %>%
  count(household_adults) %>%
  mutate(prop = n / sum(n))

cat("\n=== Distribución muestral de adultos por hogar (sin benchmark externo) ===\n")
print(dist_hogar_muestra)

# =====================================================================
# 6. Taxonomía: error muestral vs. ajeno al muestreo,
#    y por qué seguía siendo preferible el muestreo frente al censo
# =====================================================================
# Todo lo detectado en los puntos 1-5 puede clasificarse dentro de la
# distinción clásica: ERROR MUESTRAL (variabilidad inevitable por
# observar solo una parte de la población, se reduce aumentando n) vs.
# ERRORES AJENOS AL MUESTREO (no se arreglan con más observaciones:
# de cobertura, de no respuesta, de medición, de procesamiento).
# Aumentar el tamaño muestral no corrige NINGUNO de los que siguen:

taxonomia_errores <- data.frame(
  hallazgo = c(
    "Sesgo de cobertura (edad vs Census aprox.)",
    "Sesgo de complete-case (na.omit antes de clusterizar)",
    "Missingness no aleatorio en variables actitudinales",
    "Posible sobrerreporte por deseabilidad social",
    "Selección dentro del hogar sin ponderar (simulación)",
    "Variabilidad de p_hat de un cluster a otro (IC del prop.test)"
  ),
  tipo_error = c(
    "Ajeno al muestreo (cobertura)",
    "Ajeno al muestreo (no respuesta / selección)",
    "Ajeno al muestreo (no respuesta)",
    "Ajeno al muestreo (medición)",
    "Ajeno al muestreo (selección, corregible con pesos)",
    "MUESTRAL"
  ),
  se_corrige_con_mas_n = c("No", "No", "No", "No", "No", "Sí")
)

cat("\n=== Taxonomía de los errores detectados en este análisis ===\n")
print(taxonomia_errores)
cat("\nNota clave: de los 6 problemas detectados en todo el análisis, 5 son",
    "AJENOS al muestreo. Esto es importante porque recoger una muestra más",
    "grande NO los solucionaría -- solo reduciría el único error genuinamente",
    "muestral (el último de la tabla). Un error común es pensar que 'más datos'",
    "arregla cualquier problema de un estudio; aquí se ve que no es así.\n")

# Ventajas del muestreo frente al censo, aplicadas a este caso concreto
# (no como definición abstracta, sino con el argumento que de verdad
# aplica a una encuesta de gripe activa):
cat("\n=== ¿Por qué muestreo y no censo, en este caso concreto? ===\n")
cat("- Inmediatez: la NHFS se hizo DURANTE una pandemia activa (oct 2009-jun 2010).",
    "Un censo completo de EE. UU. tarda años en planificarse y ejecutarse; la",
    "utilidad de datos sobre una pandemia llegados años tarde sería nula para",
    "orientar la respuesta sanitaria en tiempo real.\n")
cat("- Coste: encuestar por RDD a una muestra de decenas de miles de personas",
    "es una fracción del coste de un censo poblacional completo.\n")
cat("- El PRECIO a pagar por esa rapidez y ese coste menor es precisamente el",
    "error muestral (variabilidad de p_hat) Y la exposición a errores ajenos",
    "al muestreo como los de la tabla de arriba, que un censo bien ejecutado",
    "no eliminaría del todo pero sí podría mitigar en parte (ej. cobertura).\n")

# =====================================================================
# 7. Muestreo probabilístico, estimador insesgado
#    y medidas de error: demostración por simulación
# =====================================================================
# Hasta ahora, en scripts anteriores, se calculaban IC y errores
# estándar de p_hat sin conectar explícitamente esos cálculos con el
# concepto teórico de "estimador insesgado". Aquí se demuestra
# empíricamente qué significa eso: si se repite el muestreo muchas
# veces, la media de las estimaciones converge al valor poblacional
# (insesgadez), y la DISPERSIÓN de esas estimaciones es precisamente
# el error muestral que la fórmula analítica intenta predecir.
#
# Tratamos aquí el dataset COMPLETO como si fuera "la población"
# (simplificación pedagógica: en la realidad esto ya es una muestra de
# la población de EE. UU., pero sirve para aislar el concepto de
# estimador insesgado sin depender de un censo real inexistente).

datos_h1n1_completos <- datos %>% filter(!is.na(h1n1_vaccine))
p_poblacional <- mean(datos_h1n1_completos$h1n1_vaccine)
N_poblacion <- nrow(datos_h1n1_completos)

cat("\n=== Simulación de muestreo repetido (tratando el dataset como 'población') ===\n")
cat("Proporción poblacional real de vacunados H1N1 (p):", round(p_poblacional, 4), "\n")

n_muestra <- 200
n_repeticiones <- 1000

set.seed(2024)
estimaciones <- replicate(n_repeticiones, {
  muestra <- sample(datos_h1n1_completos$h1n1_vaccine, size = n_muestra, replace = FALSE)
  mean(muestra)
})

media_estimaciones <- mean(estimaciones)
se_empirico <- sd(estimaciones)
se_analitico <- sqrt(p_poblacional * (1 - p_poblacional) / n_muestra)

cat("\nMedia de", n_repeticiones, "estimaciones muestrales (n =", n_muestra, "cada una):",
    round(media_estimaciones, 4), "\n")
cat("(Se acerca al valor poblacional real ", round(p_poblacional, 4),
    " -> ilustra que p_hat es un ESTIMADOR INSESGADO: no acierta en cada",
    " muestra individual, pero no tiene un sesgo sistemático hacia ningún lado)\n", sep = "")
cat("Error estándar EMPÍRICO (sd de las 1000 estimaciones):", round(se_empirico, 4), "\n")
cat("Error estándar ANALÍTICO (fórmula sqrt(p(1-p)/n)):", round(se_analitico, 4), "\n")
cat("Ambos deberían ser muy parecidos: la fórmula analítica no es magia,",
    "predice exactamente la dispersión que se observaría si de verdad se",
    "repitiera el muestreo muchas veces.\n")

ggplot(data.frame(estimacion = estimaciones), aes(x = estimacion)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  geom_vline(xintercept = p_poblacional, color = "darkred", linewidth = 1, linetype = "dashed") +
  labs(title = "Distribución muestral de p_hat (1000 muestras simuladas de n=200)",
       subtitle = paste0("Línea roja = valor poblacional real (p = ", round(p_poblacional, 3), ")"),
       x = "Proporción estimada en cada muestra", y = "Frecuencia") +
  theme_minimal()

# --------------------------------------------------------------
# Efecto de diseño (deff): por qué el SE "de manual" (SRS)
# probablemente SUBESTIMA el error real de esta encuesta concreta
# --------------------------------------------------------------
# La simulación de arriba asume muestreo aleatorio simple (SRS), que
# es el supuesto detrás de la fórmula sqrt(p(1-p)/n). Pero la NHFS es
# un diseño MÁS COMPLEJO: RDD de doble marco, con selección dentro del
# hogar (visto en el punto 5). Los diseños complejos suelen tener un
# EFECTO DE DISEÑO (deff) > 1: la varianza real de las estimaciones es
# mayor que la que predice la fórmula de SRS, típicamente porque
# personas del mismo hogar/región no son observaciones del todo
# independientes entre sí.
#
# No podemos calcular el deff exacto de esta encuesta porque no hay
# columna de hogar/cluster ni de pesos en el fichero (limitación real,
# no un cálculo omitido por pereza). Lo que SÍ se puede decir con
# honestidad metodológica es la dirección del sesgo: el error estándar
# "de manual" calculado en scripts anteriores (prop.test, IC de Wilson)
# es, casi con toda seguridad, un LÍMITE INFERIOR del error real,
# no una estimación exacta. Los IC publicados en los temas de
# proporciones deberían leerse con esa reserva.
cat("\n=== Aviso sobre el efecto de diseño (deff) ===\n")
cat("Los errores estándar e IC calculados en scripts anteriores (prop.test,",
    "Wilson) asumen muestreo aleatorio simple. La NHFS es un diseño más",
    "complejo (RDD + selección dentro del hogar), así que su error real",
    "probablemente sea ALGO MAYOR que el reportado -- sin pesos ni variable",
    "de cluster/hogar en el fichero, no se puede cuantificar el deff exacto,",
    "pero sí advertir que los IC de los temas anteriores son optimistas,",
    "no que sean incorrectos.\n")