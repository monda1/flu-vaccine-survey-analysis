## ============================================================================
## Modelo predictivo de vacunación contra influenza H1N1 y estacional
## NHFS 2009 - Script R Preprocesamiento
## ============================================================================

# Setup
library(tidyverse)
library(tidymodels)
library(recipes)
library(embed)        
library(knitr)
library(kableExtra)

training_set_features <- read_csv("data/raw/caracteristicas_datos_entrenamiento.csv")
training_set_labels   <- read_csv("data/raw/etiquetas_datos_entrenamiento.csv")
test <- read_csv("data/raw/test.csv")

# Unión
datos <- training_set_features |>
  left_join(training_set_labels, by = "respondent_id")

# Paleta
pal_azul <- "#3B82F6"
pal_rojo <- "#EF4444"
pal_gris <- "#6B7280"

opts_chunk$set(fig.width = 8, fig.height = 5, out.width = "92%")


# Estructura
tribble(
  ~Conjunto,                   ~Filas,                         ~Columnas,
  "training_set_features",     nrow(training_set_features),    ncol(training_set_features),
  "training_set_labels",       nrow(training_set_labels),      ncol(training_set_labels),
  "Datos unidos (análisis)",   nrow(datos),                    ncol(datos)
) |>
  kable(format.args = list(big.mark = ",")) |>
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE)


# Prevalencia
datos |>
  summarise(
    `N`                        = n(),
    `Vacunados H1N1 (n)`       = sum(h1n1_vaccine,     na.rm = TRUE),
    `Vacunados H1N1 (%)`       = mean(h1n1_vaccine,    na.rm = TRUE) * 100,
    `Vacunados estacional (n)` = sum(seasonal_vaccine, na.rm = TRUE),
    `Vacunados estacional (%)` = mean(seasonal_vaccine,na.rm = TRUE) * 100
  ) |>
  kable(digits = 1, format.args = list(big.mark = ",")) |>
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE)


# Barras respuesta
datos |>
  select(h1n1_vaccine, seasonal_vaccine) |>
  pivot_longer(everything(), names_to = "vacuna", values_to = "valor") |>
  mutate(
    vacuna = recode(vacuna,
                    h1n1_vaccine     = "H1N1",
                    seasonal_vaccine = "Estacional"),
    valor = factor(valor, labels = c("No vacunado", "Vacunado"))
  ) |>
  count(vacuna, valor) |>
  group_by(vacuna) |>
  mutate(pct = n / sum(n)) |>
  ggplot(aes(x = vacuna, y = pct, fill = valor)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = scales::percent(pct, accuracy = 0.1)),
            position = position_stack(vjust = 0.5), color = "white", size = 3.5) +
  scale_fill_manual(values = c("No vacunado" = pal_gris, "Vacunado" = pal_azul)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = NULL, fill = NULL,
       title = "Proporción de vacunación por etiqueta") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")


# Tipo variables
training_set_features |>
  select(-respondent_id) |>
  summarise(across(everything(), class)) |>
  pivot_longer(everything(), names_to = "Variable", values_to = "Tipo R") |>
  mutate(
    Bloque = case_when(
      Variable %in% c("h1n1_concern", "h1n1_knowledge",
                      "opinion_h1n1_vacc_effective", "opinion_h1n1_risk",
                      "opinion_h1n1_sick_from_vacc", "opinion_seas_vacc_effective",
                      "opinion_seas_risk", "opinion_seas_sick_from_vacc") ~ "1 Actitudes",
      str_starts(Variable, "behavioral_")                                 ~ "2 Comportamientos",
      Variable %in% c("doctor_recc_h1n1", "doctor_recc_seasonal",
                      "chronic_med_condition", "child_under_6_months",
                      "health_worker", "health_insurance")                ~ "3 Clínico-laboral",
      TRUE                                                                ~ "4 Sociodemográfico"
    )
  ) |>
  arrange(Bloque, Variable) |>
  kable() |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, font_size = 12) |>
  collapse_rows(columns = 3, valign = "top")


# Vacios calculo
missing_tbl <- training_set_features |>
  select(-respondent_id) |>
  summarise(across(everything(), ~ sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  mutate(
    pct_missing   = n_missing / nrow(training_set_features),
    tiene_missing = n_missing > 0
  ) |>
  arrange(desc(pct_missing))

n_vars_con_missing <- sum(missing_tbl$tiene_missing)
n_vars_sin_missing <- sum(!missing_tbl$tiene_missing)


# Vacios grafico
missing_tbl |>
  filter(tiene_missing) |>
  mutate(
    severidad = case_when(
      pct_missing >= 0.30 ~ "Alta (≥ 30 %)",
      pct_missing >= 0.10 ~ "Moderada (10–30 %)",
      TRUE                ~ "Baja (< 10 %)"
    ),
    severidad = factor(severidad,
                       levels = c("Alta (≥ 30 %)", "Moderada (10–30 %)", "Baja (< 10 %)"))
  ) |>
  ggplot(aes(x = reorder(variable, pct_missing), y = pct_missing, fill = severidad)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = scales::percent(pct_missing, accuracy = 0.1)),
            hjust = -0.08, size = 2.8) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.60),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c(
    "Alta (≥ 30 %)"      = pal_rojo,
    "Moderada (10–30 %)" = "#F59E0B",
    "Baja (< 10 %)"      = pal_azul)) +
  labs(x = NULL, y = "% de valores faltantes", fill = "Severidad",
       title    = "Missingness por variable predictora",
       subtitle = "Umbral de severidad: rojo ≥ 30 %, ámbar 10–30 %, azul < 10 %") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom")


# Vacios tabla
missing_tbl |>
  filter(tiene_missing) |>
  select(Variable = variable, `N faltantes` = n_missing, `% faltantes` = pct_missing) |>
  mutate(`% faltantes` = scales::percent(`% faltantes`, accuracy = 0.1)) |>
  kable(format.args = list(big.mark = ",")) |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, font_size = 12)


# Vacios patron
training_set_features |>
  select(-respondent_id) |>
  mutate(n_na = rowSums(is.na(across(everything())))) |>
  count(n_na) |>
  mutate(pct = n / sum(n)) |>
  ggplot(aes(x = factor(n_na), y = pct)) +
  geom_col(fill = pal_azul, alpha = 0.85) +
  geom_text(aes(label = scales::percent(pct, accuracy = 0.1)),
            vjust = -0.4, size = 3) +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Número de variables con NA en la misma observación",
       y = "% de observaciones",
       title = "Distribución de la carga de missingness por individuo") +
  theme_minimal(base_size = 11)


# Categoricas
training_set_features |>
  select(-respondent_id) |>
  select(where(is.character)) |>
  summarise(across(everything(), n_distinct)) |>
  pivot_longer(everything(), names_to = "Variable", values_to = "Niveles únicos") |>
  arrange(desc(`Niveles únicos`)) |>
  kable() |>
  kable_styling(bootstrap_options = c("striped", "hover"),
                full_width = FALSE, font_size = 12)


# Numericas resumen
training_set_features |>
  select(-respondent_id) |>
  select(where(is.numeric)) |>
  pivot_longer(everything(), names_to = "variable", values_to = "valor") |>
  group_by(variable) |>
  summarise(
    Min     = min(valor,  na.rm = TRUE),
    P25     = quantile(valor, 0.25, na.rm = TRUE),
    Mediana = median(valor, na.rm = TRUE),
    P75     = quantile(valor, 0.75, na.rm = TRUE),
    Max     = max(valor,  na.rm = TRUE),
    `N NA`  = sum(is.na(valor))
  ) |>
  kable(digits = 2, format.args = list(big.mark = ",")) |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, font_size = 11)


# Desbalance
datos |>
  summarise(
    `Ratio H1N1 (mayoritaria:minoritaria)`      =
      max(table(h1n1_vaccine))     / min(table(h1n1_vaccine)),
    `Ratio estacional (mayoritaria:minoritaria)` =
      max(table(seasonal_vaccine)) / min(table(seasonal_vaccine))
  ) |>
  kable(digits = 2) |>
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE)


# Duplicados
n_ids_unicos <- n_distinct(training_set_features$respondent_id)
n_filas      <- nrow(training_set_features)

cat("Filas totales:    ", n_filas,                "\n")
cat("IDs únicos:       ", n_ids_unicos,            "\n")
cat("Filas duplicadas: ", n_filas - n_ids_unicos,  "\n")




info_tipos_pre <- tibble(
  variable  = names(datos),
  clase     = sapply(datos, function(x) class(x)[1]),
  ordenado  = sapply(datos, is.ordered),
  n_niveles = sapply(datos, function(x)
    if (is.factor(x)) nlevels(x) else NA_integer_),
  niveles   = sapply(datos, function(x)
    if (is.factor(x)) paste(levels(x), collapse = " | ") else "-"),
  n_na      = sapply(datos, function(x) sum(is.na(x))),
  pct_na    = round(sapply(datos, function(x) mean(is.na(x))) * 100, 1)
)


cat("\n=== Todas las variables con NA antes de imputar ===\n")
info_tipos_pre |>
  filter(n_na > 0) |>
  arrange(desc(pct_na)) |>
  select(variable, clase, n_na, pct_na) |>
  kable(col.names = c("Variable", "Clase", "N NA", "% NA")) |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, font_size = 11) |>
  print()


# Conversion factor
vars_character <- c(
  "employment_occupation", "employment_industry",
  "income_poverty", "rent_or_own", "employment_status",
  "education", "marital_status",
  "sex", "race", "census_msa", "hhs_geo_region"
)

# Conversión
datos <- datos |>
  mutate(across(all_of(vars_character), as.factor))

# Tabla de niveles por variable
map_dfr(vars_character, function(v) {
  niveles <- levels(datos[[v]])
  tibble(
    Variable    = v,
    `N niveles` = length(niveles),
    Niveles     = paste(niveles, collapse = " · ")
  )
}) |>
  kable() |>
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = TRUE, font_size = 11
  ) |>
  column_spec(3, width = "55em")


# Niveles grafico
datos |>
  select(all_of(vars_character)) |>
  pivot_longer(everything(), names_to = "variable", values_to = "nivel") |>
  mutate(nivel = fct_explicit_na(nivel, na_level = "(NA)")) |>
  count(variable, nivel) |>
  group_by(variable) |>
  mutate(pct = n / sum(n)) |>
  ggplot(aes(x = reorder(nivel, pct), y = pct,
             fill = nivel == "(NA)")) +
  geom_col(alpha = 0.85, show.legend = FALSE) +
  geom_text(aes(label = scales::percent(pct, accuracy = 0.1)),
            hjust = -0.08, size = 2.4) +
  coord_flip() +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("FALSE" = pal_azul, "TRUE" = pal_rojo)) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.14))) +
  labs(x = NULL, y = NULL,
       title    = "Distribución de niveles en variables categóricas",
       subtitle = "En rojo: proporción de valores faltantes (NA)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"))


# Eliminacion hash
tribble(
  ~Variable,                ~`N niveles`, ~`% NA`,  ~Motivo,
  "employment_industry",    21L,          "49.9 %", "Niveles hash + missingness > 50 %",
  "employment_occupation",  23L,          "50.4 %", "Niveles hash + missingness > 50 %",
  "hhs_geo_region",         10L,          "0.0 %",  "Niveles hash (geografía no interpretable)"
) |>
  kable() |>
  kable_styling(bootstrap_options = c("striped", "hover"),
                full_width = FALSE)

# Eliminación del data frame principal
datos <- datos |>
  select(-employment_industry, -employment_occupation, -hhs_geo_region)

# Actualizar vector de variables categóricas activas
# Incluimos hhs_geo_region en el setdiff para que no aparezca en futuros análisis
vars_character <- setdiff(vars_character, 
                          c("employment_industry", "employment_occupation", "hhs_geo_region"))

cat(sprintf("Variables categóricas restantes: %d\n", length(vars_character)))
cat(paste(" ·", vars_character, collapse = "\n"), "\n")


# Actualizacion factor
vars_character <- c(
  "income_poverty", "rent_or_own", "employment_status",
  "education", "marital_status",
  "sex", "race", "census_msa"
)

# Conversión
datos <- datos |>
  mutate(across(all_of(vars_character), as.factor))

# Tabla de niveles por variable
map_dfr(vars_character, function(v) {
  niveles <- levels(datos[[v]])
  tibble(
    Variable    = v,
    `N niveles` = length(niveles),
    Niveles     = paste(niveles, collapse = " · ")
  )
}) |>
  kable() |>
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = TRUE, font_size = 11
  ) |>
  column_spec(3, width = "55em")


# Niveles grafico actualizado
datos |>
  select(all_of(vars_character)) |>
  pivot_longer(everything(), names_to = "variable", values_to = "nivel") |>
  mutate(nivel = fct_explicit_na(nivel, na_level = "(NA)")) |>
  count(variable, nivel) |>
  group_by(variable) |>
  mutate(pct = n / sum(n)) |>
  ggplot(aes(x = reorder(nivel, pct), y = pct,
             fill = nivel == "(NA)")) +
  geom_col(alpha = 0.85, show.legend = FALSE) +
  geom_text(aes(label = scales::percent(pct, accuracy = 0.1)),
            hjust = -0.08, size = 2.4) +
  coord_flip() +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("FALSE" = pal_azul, "TRUE" = pal_rojo)) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.14))) +
  labs(x = NULL, y = NULL,
       title    = "Distribución de niveles en variables categóricas",
       subtitle = "En rojo: proporción de valores faltantes (NA) · tras eliminar variables hash") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"))


# Imputa numericas
tribble(
  ~Variable,                     ~`% NA`,  ~Mecanismo,  ~Método,
  "health_insurance",            "46.0 %", "MAR",       "Moda (binaria 0/1)",
  "doctor_recc_h1n1",            "8.1 %",  "MAR",       "Mediana",
  "doctor_recc_seasonal",        "8.1 %",  "MAR",       "Mediana",
  "chronic_med_condition",       "3.6 %",  "MCAR",      "Mediana",
  "child_under_6_months",        "3.1 %",  "MCAR",      "Mediana",
  "health_worker",               "3.0 %",  "MCAR",      "Mediana",
  "opinion_seas_sick_from_vacc", "2.0 %",  "MCAR",      "Mediana",
  "opinion_seas_risk",           "1.9 %",  "MCAR",      "Mediana",
  "opinion_seas_vacc_effective", "1.7 %",  "MCAR",      "Mediana",
  "opinion_h1n1_vacc_effective", "1.5 %",  "MCAR",      "Mediana",
  "opinion_h1n1_risk",           "1.5 %",  "MCAR",      "Mediana",
  "opinion_h1n1_sick_from_vacc", "1.5 %",  "MCAR",      "Mediana",
  "household_adults",            "0.9 %",  "MCAR",      "Mediana",
  "household_children",          "0.9 %",  "MCAR",      "Mediana",
  "behavioral_avoidance",        "0.8 %",  "MCAR",      "Mediana",
  "behavioral_touch_face",       "0.5 %",  "MCAR",      "Mediana",
  "h1n1_knowledge",              "0.4 %",  "MCAR",      "Mediana",
  "h1n1_concern",                "0.3 %",  "MCAR",      "Mediana",
  "behavioral_antiviral_meds",   "0.3 %",  "MCAR",      "Mediana",
  "behavioral_large_gatherings", "0.3 %",  "MCAR",      "Mediana",
  "behavioral_outside_home",     "0.3 %",  "MCAR",      "Mediana",
  "behavioral_wash_hands",       "0.2 %",  "MCAR",      "Mediana",
  "behavioral_face_mask",        "0.1 %",  "MCAR",      "Mediana"
) |>
  kable() |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, font_size = 11)


# Imputa numericas viz
vars_plot_num <- c("doctor_recc_h1n1", "health_insurance",
                   "opinion_seas_risk", "chronic_med_condition",
                   "household_adults",  "h1n1_concern")

datos_num_imp <- datos |>
  mutate(across(
    where(is.numeric) & !all_of(c("respondent_id",
                                  "h1n1_vaccine", "seasonal_vaccine")),
    ~ if_else(is.na(.), median(., na.rm = TRUE), .)
  )) |>
  # health_insurance: moda (binaria)
  mutate(health_insurance = if_else(
    is.na(health_insurance),
    as.numeric(names(which.max(table(health_insurance)))),
    health_insurance
  ))

bind_rows(
  datos         |> select(all_of(vars_plot_num)) |> mutate(etapa = "Antes"),
  datos_num_imp |> select(all_of(vars_plot_num)) |> mutate(etapa = "Después")
) |>
  pivot_longer(-etapa, names_to = "variable", values_to = "valor") |>
  drop_na(valor) |>
  ggplot(aes(x = factor(valor), fill = etapa)) +
  geom_bar(position = "dodge", alpha = 0.85) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("Antes" = pal_gris, "Después" = pal_azul)) +
  labs(x = NULL, y = "Frecuencia", fill = NULL,
       title    = "Efecto de la imputación por mediana/moda en variables numéricas",
       subtitle = "La distribución post-imputación debe aproximar la pre-imputación") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text      = element_text(face = "bold"),
        axis.text.x     = element_text(angle = 30, hjust = 1))


# Imputa factors
tribble(
  ~Variable,           ~`% NA`,  ~Mecanismo, ~Método,
  "income_poverty",    "16.6 %", "MAR",      "Hot Deck KNN (k = 5)",
  "rent_or_own",       "7.6 %",  "MCAR",     "Moda",
  "employment_status", "5.5 %",  "MCAR",     "Moda",
  "education",         "5.3 %",  "MCAR",     "Moda",
  "marital_status",    "5.3 %",  "MCAR",     "Moda"
) |>
  kable() |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, font_size = 11)


# Imputa factors viz
vars_factor_imp <- c("income_poverty", "rent_or_own",
                     "employment_status", "education", "marital_status")

moda <- function(x) {
  ux <- na.omit(unique(x))
  ux[which.max(tabulate(match(x, ux)))]
}

datos_fct_imp <- datos |>
  mutate(across(all_of(vars_factor_imp),
                ~ if_else(is.na(.), moda(.), .)))

bind_rows(
  datos         |> select(all_of(vars_factor_imp)) |>
    mutate(across(everything(), as.character), etapa = "Antes"),
  datos_fct_imp |> select(all_of(vars_factor_imp)) |>
    mutate(across(everything(), as.character), etapa = "Después")
) |>
  pivot_longer(-etapa, names_to = "variable", values_to = "valor") |>
  drop_na(valor) |>
  count(etapa, variable, valor) |>
  group_by(etapa, variable) |>
  mutate(pct = n / sum(n)) |>
  ggplot(aes(x = reorder(valor, pct), y = pct, fill = etapa)) +
  geom_col(position = "dodge", alpha = 0.85) +
  geom_text(aes(label = scales::percent(pct, accuracy = 0.1)),
            position = position_dodge(width = 0.9),
            hjust = -0.08, size = 2.3) +
  coord_flip() +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Antes" = pal_gris, "Después" = pal_azul)) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = NULL, fill = NULL,
       title    = "Efecto de la imputación en variables categóricas",
       subtitle = "Azul = post-imputación · Gris = pre-imputación (sin NA)") +
  theme_minimal(base_size = 10) +
  theme(legend.position    = "bottom",
        strip.text         = element_text(face = "bold"),
        panel.grid.major.y = element_blank())


# Verificacion nulos

moda <- function(x) {
  ux <- na.omit(unique(x))
  ux[which.max(tabulate(match(x, ux)))]
}

datos_imputado <- datos |>
  # Numéricas: mediana (todas excepto id y respuestas)
  mutate(across(
    where(is.numeric) &
      !all_of(c("respondent_id", "h1n1_vaccine", "seasonal_vaccine")),
    ~ if_else(is.na(.), median(., na.rm = TRUE), .)
  )) |>
  # health_insurance: moda (binaria)
  mutate(health_insurance = if_else(
    is.na(health_insurance),
    as.numeric(names(which.max(table(health_insurance)))),
    health_insurance
  )) |>
  # Factores: moda
  mutate(across(
    all_of(c("rent_or_own", "employment_status",
             "education", "marital_status", "income_poverty")),
    ~ if_else(is.na(.), moda(.), .)
  ))

# Conteo de NA por variable
na_post <- datos_imputado |>
  select(-respondent_id) |>
  summarise(across(everything(), ~ sum(is.na(.)))) |>
  pivot_longer(everything(),
               names_to  = "Variable",
               values_to = "NA residuales") |>
  filter(`NA residuales` > 0)

# Resultado
if (nrow(na_post) == 0) {
  cat("Ningún valor faltante residual. El dataset está listo para el modelado.\n")
} else {
  cat("Variables con NA residual:\n")
  na_post |>
    kable() |>
    kable_styling(bootstrap_options = c("striped", "hover"),
                  full_width = FALSE) |>
    print()
}

# Tabla resumen final
tibble(
  Métrica = c(
    "Observaciones",
    "Variables eliminadas (hash)",
    "Variables en el dataset final",
    "NA totales antes de imputar",
    "NA totales después de imputar"
  ),
  Valor = c(
    nrow(datos_imputado),
    3L,
    ncol(datos_imputado |> select(-respondent_id)),
    sum(is.na(datos |> select(-respondent_id))),
    sum(is.na(datos_imputado |> select(-respondent_id)))
  )
) |>
  kable() |>
  kable_styling(bootstrap_options = c("striped", "hover"),
                full_width = FALSE)



# Variables que tenían NA originalmente
vars_con_na <- datos |>
  select(-respondent_id) |>
  summarise(across(everything(), ~ sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_na") |>
  filter(n_na > 0) |>
  pull(variable)

bind_rows(
  datos |>
    select(all_of(vars_con_na)) |>
    summarise(across(everything(), ~ sum(is.na(.)))) |>
    pivot_longer(everything(),
                 names_to  = "variable",
                 values_to = "n_na") |>
    mutate(etapa = "Antes"),
  datos_imputado |>
    select(any_of(vars_con_na)) |>
    summarise(across(everything(), ~ sum(is.na(.)))) |>
    pivot_longer(everything(),
                 names_to  = "variable",
                 values_to = "n_na") |>
    mutate(etapa = "Después")
) |>
  ggplot(aes(x = reorder(variable, n_na), y = n_na, fill = etapa)) +
  geom_col(position = "dodge", alpha = 0.85) +
  geom_text(aes(label = n_na),
            position = position_dodge(width = 0.9),
            hjust = -0.1, size = 2.8) +
  coord_flip() +
  scale_fill_manual(values = c("Antes" = pal_rojo, "Después" = pal_azul)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "N valores faltantes", fill = NULL,
       title    = "Missingness antes y después de la imputación",
       subtitle = "Azul = post-imputación (debe ser 0 en todas las variables)") +
  theme_minimal(base_size = 10) +
  theme(legend.position    = "bottom",
        panel.grid.major.y = element_blank(),
        strip.text         = element_text(face = "bold"))


#  diccionario actitudes
tribble(
  ~Variable,                     ~Tipo,     ~Escala, ~Descripción,
  "h1n1_concern",                "Ordinal", "0–3",
  "Nivel de preocupación del encuestado por la enfermedad H1N1. 0 = nada preocupado · 3 = muy preocupado",
  "h1n1_knowledge",              "Ordinal", "0–2",
  "Nivel de conocimiento autopercibido sobre el H1N1. 0 = ninguno · 2 = mucho",
  "opinion_h1n1_vacc_effective", "Ordinal", "1–5",
  "Opinión sobre la eficacia de la vacuna H1N1. 1 = nada efectiva · 5 = muy efectiva",
  "opinion_h1n1_risk",           "Ordinal", "1–5",
  "Opinión sobre el riesgo de enfermar de H1N1 sin vacuna. 1 = muy bajo · 5 = muy alto",
  "opinion_h1n1_sick_from_vacc", "Ordinal", "1–5",
  "Preocupación por enfermar a causa de la vacuna H1N1. 1 = nada · 5 = mucho",
  "opinion_seas_vacc_effective", "Ordinal", "1–5",
  "Opinión sobre la eficacia de la vacuna estacional. 1 = nada efectiva · 5 = muy efectiva",
  "opinion_seas_risk",           "Ordinal", "1–5",
  "Opinión sobre el riesgo de enfermar de gripe estacional sin vacuna. 1 = muy bajo · 5 = muy alto",
  "opinion_seas_sick_from_vacc", "Ordinal", "1–5",
  "Preocupación por enfermar a causa de la vacuna estacional. 1 = nada · 5 = mucho"
) |>
  kable() |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = TRUE, font_size = 11) |>
  column_spec(4, width = "30em")


# Diccionario comportamientos
tribble(
  ~Variable,                     ~Escala, ~Descripción,
  "behavioral_antiviral_meds",   "0 / 1",
  "Toma de medicamentos antivirales como medida preventiva frente al H1N1",
  "behavioral_avoidance",        "0 / 1",
  "Evitar el contacto con personas con síntomas de gripe",
  "behavioral_face_mask",        "0 / 1",
  "Uso de mascarilla facial",
  "behavioral_wash_hands",       "0 / 1",
  "Lavado frecuente de manos o uso de desinfectante",
  "behavioral_large_gatherings", "0 / 1",
  "Reducción de la asistencia a reuniones de gran aforo",
  "behavioral_outside_home",     "0 / 1",
  "Reducción del tiempo fuera del hogar",
  "behavioral_touch_face",       "0 / 1",
  "Evitar tocarse la cara con las manos"
) |>
  kable() |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = TRUE, font_size = 11) |>
  column_spec(3, width = "35em")


# Diccionario clinico
tribble(
  ~Variable,               ~Escala, ~Descripción,
  "doctor_recc_h1n1",      "0 / 1",
  "El médico recomendó explícitamente al encuestado vacunarse contra el H1N1",
  "doctor_recc_seasonal",  "0 / 1",
  "El médico recomendó explícitamente vacunarse contra la gripe estacional",
  "chronic_med_condition", "0 / 1",
  "El encuestado padece alguna enfermedad crónica (diabetes, asma, cardiopatía, etc.)",
  "child_under_6_months",  "0 / 1",
  "Convivencia con un menor de 6 meses en el hogar",
  "health_worker",         "0 / 1",
  "El encuestado trabaja en el sector sanitario",
  "health_insurance",      "0 / 1",
  "El encuestado dispone de seguro médico privado o público"
) |>
  kable() |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = TRUE, font_size = 11) |>
  column_spec(3, width = "35em")


# Diccionario socio
tribble(
  ~Variable,           ~Tipo,       ~Categorías, ~Descripción,
  "age_group",         "Nominal",   "5 grupos",
  "Grupo de edad: 18–34 · 35–44 · 45–54 · 55–64 · 65+",
  "education",         "Ordinal",   "4 niveles",
  "Nivel educativo: < 12 Years · 12 Years · Some College · College Graduate",
  "race",              "Nominal",   "4 categorías",
  "Raza/etnia autodeclarada: White · Black · Hispanic · Other or Multiple",
  "sex",               "Binaria",   "Female / Male",
  "Sexo del encuestado",
  "income_poverty",    "Ordinal",   "3 niveles",
  "Ingresos del hogar respecto al umbral de pobreza: Below Poverty · <= $75,000 Above Poverty · > $75,000",
  "marital_status",    "Binaria",   "2 categorías",
  "Estado civil: Married · Not Married",
  "rent_or_own",       "Binaria",   "2 categorías",
  "Régimen de tenencia de la vivienda: Own · Rent",
  "employment_status", "Nominal",   "3 categorías",
  "Situación laboral: Employed · Unemployed · Not in Labor Force",
  "census_msa",        "Nominal",   "3 categorías",
  "Tipo de área censal (MSA): Principle City · Not Principle City · Non-MSA",
  "household_adults",  "Discreta",  "0–3+",
  "Número de adultos (≥ 18 años) en el hogar, excluido el encuestado",
  "household_children","Discreta",  "0–3+",
  "Número de menores de 18 años en el hogar"
) |>
  kable() |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = TRUE, font_size = 11) |>
  column_spec(4, width = "28em")


# Diccionario respuesta
tribble(
  ~Variable,          ~Tipo,    ~Escala, ~Descripción,
  "h1n1_vaccine",     "Binaria","0 / 1",
  "El encuestado recibió la vacuna contra el H1N1. 0 = no · 1 = sí",
  "seasonal_vaccine", "Binaria","0 / 1",
  "El encuestado recibió la vacuna contra la gripe estacional. 0 = no · 1 = sí"
) |>
  kable() |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = TRUE, font_size = 11) |>
  column_spec(4, width = "35em")
