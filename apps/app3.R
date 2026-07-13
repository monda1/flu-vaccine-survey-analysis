#
# Shiny App: Muestreo — población objetivo vs muestreada, sesgos de
# cobertura/selección/medición, unidad de muestreo vs observación,
# taxonomía de errores y estimador insesgado
# Datos: National 2009 H1N1 Flu Survey (DrivenData - Flu Shot Learning)
#
# Ejecuta con el botón 'Run App'.
#

library(shiny)
library(dplyr)
library(tidyr)
library(readr)
library(bslib)
library(bsicons)
library(ggplot2)
library(scales)

# --------------------------------------------------------------
# 1. Carga de datos (ruta robusta)
# --------------------------------------------------------------
encontrar_ruta_datos <- function(archivo) {
  candidatos <- c(
    file.path("data", "raw", archivo),
    file.path("..", "data", "raw", archivo)
  )
  ruta <- candidatos[file.exists(candidatos)]
  if (length(ruta) == 0) {
    stop(
      "No se encontró '", archivo, "'. Colócalo en 'data/raw/' en la raíz ",
      "de tu proyecto.\nDirectorio de trabajo actual: ", getwd()
    )
  }
  ruta[1]
}

features <- read_csv(encontrar_ruta_datos("training_set_features.csv"), show_col_types = FALSE)
labels   <- read_csv(encontrar_ruta_datos("training_set_labels.csv"),   show_col_types = FALSE)
datos <- features %>% left_join(labels, by = "respondent_id")

vars_cluster <- c("h1n1_concern", "h1n1_knowledge",
                  "opinion_h1n1_vacc_effective", "opinion_h1n1_risk", "opinion_h1n1_sick_from_vacc",
                  "opinion_seas_vacc_effective", "opinion_seas_risk", "opinion_seas_sick_from_vacc",
                  "behavioral_antiviral_meds", "behavioral_avoidance", "behavioral_face_mask",
                  "behavioral_wash_hands", "behavioral_large_gatherings",
                  "behavioral_outside_home", "behavioral_touch_face",
                  "doctor_recc_h1n1", "doctor_recc_seasonal")
vars_cluster <- intersect(vars_cluster, names(datos))

vars_demograficas <- intersect(
  c("age_group", "education", "income_poverty", "sex", "race", "employment_status"),
  names(datos)
)

indices_completos <- datos %>%
  dplyr::select(all_of(vars_cluster)) %>%
  complete.cases() %>%
  which()
datos$en_complete_case <- FALSE
datos$en_complete_case[indices_completos] <- TRUE

ref_edad_censo <- c(
  "18 - 34 Years" = 0.30, "35 - 44 Years" = 0.18, "45 - 54 Years" = 0.19,
  "55 - 64 Years" = 0.16, "65+ Years" = 0.17
)

taxonomia_errores <- data.frame(
  Hallazgo = c(
    "Sesgo de cobertura (edad vs Census aprox.)",
    "Sesgo de complete-case (na.omit antes de clusterizar)",
    "Missingness no aleatorio en variables actitudinales",
    "Posible sobrerreporte por deseabilidad social",
    "Selección dentro del hogar sin ponderar (simulación)",
    "Variabilidad de p̂ entre muestras (IC del prop.test)"
  ),
  `Tipo de error` = c(
    "Ajeno al muestreo (cobertura)",
    "Ajeno al muestreo (no respuesta / selección)",
    "Ajeno al muestreo (no respuesta)",
    "Ajeno al muestreo (medición)",
    "Ajeno al muestreo (selección, corregible con pesos)",
    "MUESTRAL"
  ),
  `¿Se corrige con más n?` = c("No", "No", "No", "No", "No", "Sí"),
  check.names = FALSE
)

# --------------------------------------------------------------
# 2. UI
# --------------------------------------------------------------
ui <- page_sidebar(
  title = "Muestreo · Sesgos, unidad muestral y estimador insesgado",
  theme = bs_theme(
    version = 5,
    bootswatch = "minty",
    base_font = font_google("Inter"),
    heading_font = font_google("Inter"),
    primary = "#2C6E63"
  ),
  
  sidebar = sidebar(
    width = 300,
    selectInput(
      "var_demo", "Variable demográfica a examinar:",
      choices = vars_demograficas,
      selected = if ("income_poverty" %in% vars_demograficas) "income_poverty" else vars_demograficas[1]
    ),
    selectInput(
      "var_objetivo", "Variable actitudinal/clínica (missingness):",
      choices = vars_cluster,
      selected = vars_cluster[1]
    ),
    hr(),
    sliderInput("n_hogares_sim", "Simulación hogar: nº de hogares",
                min = 2000, max = 20000, value = 10000, step = 1000),
    sliderInput("n_hogares_muestreados", "Simulación hogar: hogares muestreados",
                min = 200, max = 5000, value = 2000, step = 200),
    hr(),
    sliderInput("n_muestra_est", "Estimador: tamaño de cada muestra (n)",
                min = 50, max = 1000, value = 200, step = 50),
    sliderInput("n_repeticiones_est", "Estimador: nº de repeticiones",
                min = 200, max = 5000, value = 1000, step = 200),
    hr(),
    selectInput(
      "vacuna_estrat", "Muestreo estratificado — vacuna:",
      choices = c("H1N1" = "h1n1_vaccine", "Estacional" = "seasonal_vaccine"),
      selected = "h1n1_vaccine"
    ),
    selectInput(
      "var_estrato", "Muestreo estratificado — variable de estratos:",
      choices = vars_demograficas,
      selected = if ("income_poverty" %in% vars_demograficas) "income_poverty" else vars_demograficas[1]
    ),
    sliderInput("n_nuevo_diseno", "Tamaño del nuevo diseño (n)",
                min = 1000, max = 50000, value = 10000, step = 1000),
    hr(),
    p(class = "text-muted small",
      "La NHFS es una encuesta telefónica RDD sin columna de pesos muestrales
       ni de hogar/cluster en este fichero, por lo que el análisis se limita
       a lo verificable directamente con los datos disponibles.")
  ),
  
  navset_card_underline(
    nav_panel(
      "Cobertura poblacional",
      card_body(
        plotOutput("plot_cobertura", height = "360px"),
        tableOutput("tabla_cobertura"),
        p(class = "text-muted small mt-2",
          "Comparación de la distribución muestral por edad frente a una referencia
           aproximada del Census Bureau (2010). Diferencias grandes sugieren sesgo
           de cobertura (quién tiene más probabilidad de responder a una RDD).")
      )
    ),
    nav_panel(
      "Sesgo complete-case",
      card_body(
        plotOutput("plot_attrition", height = "320px"),
        tableOutput("tabla_attrition_test"),
        p(class = "text-muted small mt-2",
          sprintf("Se retienen %d de %d filas (%.1f%%) tras exigir datos completos en
           las %d variables usadas para clustering. La tabla contrasta, por
           chi-cuadrado, si esa retención depende de cada variable demográfica.",
                  length(indices_completos), nrow(datos),
                  100 * length(indices_completos) / nrow(datos), length(vars_cluster)))
      )
    ),
    nav_panel(
      "Missingness no aleatorio",
      card_body(
        plotOutput("plot_missingness", height = "320px"),
        tableOutput("tabla_missingness_demo"),
        p(class = "text-muted small mt-2",
          "Arriba: % de valores perdidos por variable actitudinal/clínica. Abajo:
           test χ² de si la falta de respuesta en la variable objetivo elegida está
           asociada a cada variable demográfica (p < 0.05 sugiere no aleatoriedad).")
      )
    ),
    nav_panel(
      "Sesgo de medición",
      card_body(
        plotOutput("plot_medicion", height = "320px"),
        tableOutput("tabla_medicion"),
        p(class = "text-muted small mt-2",
          "Indicio indirecto (no prueba) de sesgo de medición por autodeclaración:
           coherencia entre recomendación médica y vacunación declarada, y gradiente
           de vacunación según nivel de preocupación por la H1N1.")
      )
    ),
    nav_panel(
      "Unidad de muestreo vs observación",
      card_body(
        plotOutput("plot_hogares", height = "340px"),
        tableOutput("tabla_hogares"),
        p(class = "text-muted small mt-2",
          "Simulación autocontenida (no usa los datos reales): si se selecciona 1
           persona al azar por hogar sin ponderar por tamaño, las personas de
           hogares grandes quedan infrarrepresentadas frente a la población real
           de individuos.")
      )
    ),
    nav_panel(
      "Estimador insesgado",
      card_body(
        plotOutput("plot_estimador", height = "320px"),
        tableOutput("tabla_estimador"),
        hr(),
        h6("Taxonomía de errores detectados"),
        tableOutput("tabla_taxonomia"),
        p(class = "text-muted small mt-2",
          "El SE analítico √(p(1-p)/n) asume muestreo aleatorio simple. La NHFS
           tiene un diseño más complejo (RDD + selección dentro del hogar), así que
           su error real probablemente sea algo mayor (efecto de diseño, deff > 1)
           — sin variable de cluster/peso no se puede cuantificar exactamente.")
      )
    ),
    nav_panel(
      "Muestreo estratificado",
      card_body(
        layout_columns(
          col_widths = c(3, 3, 3, 3),
          value_box(title = "p̂ MAS", value = textOutput("p_mas_val"),
                    showcase = bs_icon("dice-5"), theme = "secondary"),
          value_box(title = "p̂ estratificado", value = textOutput("p_st_val"),
                    showcase = bs_icon("layers"), theme = "primary"),
          value_box(title = "η² (estrato explica)", value = textOutput("eta2_val"),
                    showcase = bs_icon("pie-chart"), theme = "success"),
          value_box(title = "Deff", value = textOutput("deff_val"),
                    showcase = bs_icon("speedometer2"), theme = "warning")
        ),
        plotOutput("plot_estratos", height = "320px"),
        tableOutput("tabla_estratos"),
        hr(),
        h6("Comparación de esquemas de asignación"),
        plotOutput("plot_asignacion", height = "320px"),
        tableOutput("tabla_varianza_asignacion"),
        p(class = "text-muted small mt-2",
          "MAS ignora los estratos; el estimador estratificado pondera cada estrato
           por su peso muestral (W_h = n_h/n). η² indica qué parte de la variabilidad
           de la vacunación explica la variable de estratificación elegida. Deff < 1
           significa que el diseño estratificado es más eficiente que un MAS del mismo
           tamaño; el tamaño de muestra efectivo (tabla) traduce esa ganancia a
           'observaciones de MAS equivalentes'. Abajo: para un nuevo diseño con el n
           elegido en el sidebar, se compara la varianza esperada bajo asignación
           igualitaria, proporcional y óptima (Neyman) — Neyman debería minimizarla.")
      )
    )
  )
)

# --------------------------------------------------------------
# 3. Server
# --------------------------------------------------------------
server <- function(input, output, session) {
  
  # ---- 1. Cobertura poblacional ----
  dist_edad <- reactive({
    datos %>%
      filter(!is.na(age_group)) %>%
      count(age_group) %>%
      mutate(prop_muestra = n / sum(n)) %>%
      mutate(ref_censo = ref_edad_censo[as.character(age_group)])
  })
  
  output$plot_cobertura <- renderPlot({
    df <- dist_edad() %>%
      select(age_group, Muestra = prop_muestra, `Referencia Census` = ref_censo) %>%
      pivot_longer(-age_group, names_to = "fuente", values_to = "prop")
    
    ggplot(df, aes(x = age_group, y = prop, fill = fuente)) +
      geom_col(position = "dodge", width = 0.6) +
      scale_y_continuous(labels = percent) +
      scale_fill_manual(values = c("Muestra" = "#2C6E63", "Referencia Census" = "#D9534F")) +
      labs(x = "Grupo de edad", y = "Proporción", fill = NULL) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
  })
  
  output$tabla_cobertura <- renderTable({
    dist_edad() %>%
      mutate(prop_muestra = round(prop_muestra, 4),
             ref_censo = round(ref_censo, 4),
             diferencia = round(prop_muestra - ref_censo, 4)) %>%
      select(`Grupo de edad` = age_group, n, `Prop. muestra` = prop_muestra,
             `Ref. Census` = ref_censo, Diferencia = diferencia)
  })
  
  # ---- 2. Sesgo complete-case ----
  test_attrition <- function(var_demo) {
    df <- datos %>% filter(!is.na(.data[[var_demo]]))
    tabla <- table(df$en_complete_case, df[[var_demo]])
    test <- suppressWarnings(chisq.test(tabla))
    data.frame(Variable = var_demo,
               `Chi-cuadrado` = round(unname(test$statistic), 2),
               `p-valor` = format.pval(test$p.value, digits = 4),
               `¿Depende de esta variable?` = ifelse(test$p.value < 0.05, "Sí", "No"),
               check.names = FALSE)
  }
  
  output$tabla_attrition_test <- renderTable({
    bind_rows(lapply(vars_demograficas, test_attrition))
  })
  
  output$plot_attrition <- renderPlot({
    var_demo <- input$var_demo
    df <- datos %>%
      filter(!is.na(.data[[var_demo]])) %>%
      group_by(en_complete_case, valor = .data[[var_demo]]) %>%
      summarise(n = n(), .groups = "drop_last") %>%
      mutate(prop = n / sum(n)) %>%
      ungroup()
    
    ggplot(df, aes(x = valor, y = prop, fill = en_complete_case)) +
      geom_col(position = "dodge", width = 0.6) +
      scale_y_continuous(labels = percent) +
      scale_fill_manual(values = c(`TRUE` = "#2C6E63", `FALSE` = "#D9534F"),
                        labels = c(`TRUE` = "En complete-case", `FALSE` = "Excluido")) +
      labs(x = var_demo, y = "Proporción dentro de cada grupo", fill = NULL) +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            panel.grid.minor = element_blank())
  })
  
  # ---- 3. Missingness ----
  output$plot_missingness <- renderPlot({
    tabla <- sapply(vars_cluster, function(v) mean(is.na(datos[[v]])))
    df <- data.frame(variable = names(tabla), pct = 100 * tabla) %>% arrange(desc(pct))
    ggplot(df, aes(x = reorder(variable, pct), y = pct)) +
      geom_col(fill = "#2C6E63", width = 0.6) +
      coord_flip() +
      labs(x = NULL, y = "% de valores perdidos") +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank())
  })
  
  output$tabla_missingness_demo <- renderTable({
    var_obj <- input$var_objetivo
    resultados <- lapply(vars_demograficas, function(vd) {
      df <- datos %>%
        filter(!is.na(.data[[vd]])) %>%
        mutate(perdido = is.na(.data[[var_obj]]))
      tabla <- table(df$perdido, df[[vd]])
      if (nrow(tabla) < 2) {
        return(data.frame(`Variable demográfica` = vd, `Chi-cuadrado` = NA,
                          `p-valor` = "—", `¿No aleatorio?` = "n/a", check.names = FALSE))
      }
      test <- suppressWarnings(chisq.test(tabla))
      data.frame(`Variable demográfica` = vd,
                 `Chi-cuadrado` = round(unname(test$statistic), 2),
                 `p-valor` = format.pval(test$p.value, digits = 4),
                 `¿No aleatorio?` = ifelse(test$p.value < 0.05, "Sí (MAR/MNAR)", "No evidencia"),
                 check.names = FALSE)
    })
    bind_rows(resultados)
  })
  
  # ---- 4. Sesgo de medición ----
  output$plot_medicion <- renderPlot({
    coherencia <- datos %>%
      filter(!is.na(doctor_recc_h1n1), !is.na(h1n1_vaccine)) %>%
      group_by(doctor_recc_h1n1) %>%
      summarise(tasa_vacunacion = mean(h1n1_vaccine), .groups = "drop")
    
    ggplot(coherencia, aes(x = factor(doctor_recc_h1n1), y = tasa_vacunacion)) +
      geom_col(fill = "#2C6E63", width = 0.5) +
      scale_y_continuous(labels = percent) +
      labs(x = "¿Recomendación médica de vacunarse?", y = "Tasa de vacunación autodeclarada") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
  })
  
  output$tabla_medicion <- renderTable({
    datos %>%
      filter(!is.na(h1n1_concern), !is.na(h1n1_vaccine)) %>%
      group_by(`Nivel de preocupación` = h1n1_concern) %>%
      summarise(n = n(), `Tasa de vacunación` = round(mean(h1n1_vaccine), 4), .groups = "drop")
  })
  
  # ---- 5. Unidad de muestreo vs observación ----
  hogares_sim <- reactive({
    set.seed(789)
    n_hogares <- input$n_hogares_sim
    n_muestreados <- min(input$n_hogares_muestreados, n_hogares)
    
    poblacion_hogares <- data.frame(
      hogar_id = 1:n_hogares,
      tamano_hogar = sample(1:4, n_hogares, replace = TRUE, prob = c(0.30, 0.35, 0.20, 0.15))
    )
    poblacion_personas <- poblacion_hogares %>%
      tidyr::uncount(tamano_hogar) %>%
      mutate(persona_id = row_number())
    
    hogares_muestreados <- poblacion_hogares[sample(1:nrow(poblacion_hogares), n_muestreados), ]
    
    list(
      hogares = round(prop.table(table(poblacion_hogares$tamano_hogar)), 4),
      seleccionadas = round(prop.table(table(hogares_muestreados$tamano_hogar)), 4),
      personas = round(prop.table(table(poblacion_personas$tamano_hogar)), 4)
    )
  })
  
  output$plot_hogares <- renderPlot({
    s <- hogares_sim()
    df <- bind_rows(
      data.frame(tamano = names(s$hogares), prop = as.numeric(s$hogares), fuente = "Hogares (población)"),
      data.frame(tamano = names(s$seleccionadas), prop = as.numeric(s$seleccionadas), fuente = "Personas seleccionadas (sin ponderar)"),
      data.frame(tamano = names(s$personas), prop = as.numeric(s$personas), fuente = "Personas (población real)")
    )
    ggplot(df, aes(x = tamano, y = prop, fill = fuente)) +
      geom_col(position = "dodge", width = 0.7) +
      scale_y_continuous(labels = percent) +
      scale_fill_manual(values = c("Hogares (población)" = "#8FC1B5",
                                   "Personas seleccionadas (sin ponderar)" = "#2C6E63",
                                   "Personas (población real)" = "#D9534F")) +
      labs(x = "Tamaño del hogar", y = "Proporción", fill = NULL) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom", panel.grid.minor = element_blank())
  })
  
  output$tabla_hogares <- renderTable({
    s <- hogares_sim()
    data.frame(
      `Tamaño hogar` = names(s$hogares),
      `Hogares (pob.)` = as.numeric(s$hogares),
      `Personas seleccionadas` = as.numeric(s$seleccionadas),
      `Personas (pob. real)` = as.numeric(s$personas),
      check.names = FALSE
    )
  })
  
  # ---- 6. Estimador insesgado ----
  estimador_sim <- reactive({
    datos_h1n1 <- datos %>% filter(!is.na(h1n1_vaccine))
    p_poblacional <- mean(datos_h1n1$h1n1_vaccine)
    n_muestra <- input$n_muestra_est
    n_rep <- input$n_repeticiones_est
    
    set.seed(2024)
    estimaciones <- replicate(n_rep, {
      muestra <- sample(datos_h1n1$h1n1_vaccine, size = n_muestra, replace = FALSE)
      mean(muestra)
    })
    
    list(
      estimaciones = estimaciones,
      p_poblacional = p_poblacional,
      media = mean(estimaciones),
      se_empirico = sd(estimaciones),
      se_analitico = sqrt(p_poblacional * (1 - p_poblacional) / n_muestra)
    )
  })
  
  output$plot_estimador <- renderPlot({
    r <- estimador_sim()
    ggplot(data.frame(estimacion = r$estimaciones), aes(x = estimacion)) +
      geom_histogram(bins = 30, fill = "#8FC1B5", color = "white") +
      geom_vline(xintercept = r$p_poblacional, color = "#D9534F", linewidth = 1, linetype = "dashed") +
      labs(x = "Proporción estimada en cada muestra", y = "Frecuencia") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
  })
  
  output$tabla_estimador <- renderTable({
    r <- estimador_sim()
    data.frame(
      Métrica = c("p poblacional (dataset completo)", "Media de las estimaciones",
                  "Error estándar EMPÍRICO (sd)", "Error estándar ANALÍTICO √(p(1-p)/n)"),
      Valor = sprintf("%.4f", c(r$p_poblacional, r$media, r$se_empirico, r$se_analitico))
    )
  })
  
  output$tabla_taxonomia <- renderTable({
    taxonomia_errores
  })
  
  # ---- 7. Muestreo estratificado ----
  estratificado_reactivo <- reactive({
    var_y <- input$vacuna_estrat
    var_h <- input$var_estrato
    
    df <- datos %>% filter(!is.na(.data[[var_h]]), !is.na(.data[[var_y]]))
    n_total <- nrow(df)
    
    p_mas <- mean(df[[var_y]])
    var_mas <- (p_mas * (1 - p_mas)) / (n_total - 1)
    
    estratos <- df %>%
      group_by(estrato = .data[[var_h]]) %>%
      summarise(
        n_h = n(),
        W_h = n_h / n_total,
        p_h = mean(.data[[var_y]]),
        var_h = (p_h * (1 - p_h)) / pmax(n_h - 1, 1),
        S_h = sqrt(p_h * (1 - p_h)),
        .groups = "drop"
      )
    
    p_estratificado <- sum(estratos$W_h * estratos$p_h)
    var_estratificada <- sum((estratos$W_h^2) * estratos$var_h)
    
    media_global <- mean(df[[var_y]])
    SST <- sum((df[[var_y]] - media_global)^2)
    SSB <- sum(estratos$n_h * (estratos$p_h - media_global)^2)
    eta2 <- SSB / SST
    
    deff <- var_estratificada / var_mas
    n_efectivo <- n_total / deff
    
    z <- qnorm(0.975)
    ic_inf <- p_estratificado - z * sqrt(var_estratificada)
    ic_sup <- p_estratificado + z * sqrt(var_estratificada)
    
    # Asignación para un nuevo diseño
    n_nuevo <- input$n_nuevo_diseno
    L <- nrow(estratos)
    estratos_asig <- estratos %>%
      mutate(
        n_igual = round(n_nuevo / L),
        n_proporcional = round(n_nuevo * W_h),
        numerador_neyman = W_h * S_h,
        n_neyman = round(n_nuevo * (numerador_neyman / sum(numerador_neyman)))
      )
    
    var_teorica <- function(n_h_vector) sum((estratos$W_h^2) * (estratos$S_h^2) / pmax(n_h_vector, 1))
    comparacion_varianzas <- data.frame(
      Esquema = c("Igualitaria", "Proporcional", "Neyman"),
      `Varianza esperada` = c(
        var_teorica(estratos_asig$n_igual),
        var_teorica(estratos_asig$n_proporcional),
        var_teorica(estratos_asig$n_neyman)
      ),
      check.names = FALSE
    ) %>% arrange(`Varianza esperada`)
    
    list(
      n_total = n_total, p_mas = p_mas, var_mas = var_mas,
      estratos = estratos, p_estratificado = p_estratificado,
      var_estratificada = var_estratificada, eta2 = eta2, deff = deff,
      n_efectivo = n_efectivo, ic_inf = ic_inf, ic_sup = ic_sup,
      estratos_asig = estratos_asig, comparacion_varianzas = comparacion_varianzas
    )
  })
  
  output$p_mas_val <- renderText(sprintf("%.3f", estratificado_reactivo()$p_mas))
  output$p_st_val  <- renderText(sprintf("%.3f", estratificado_reactivo()$p_estratificado))
  output$eta2_val  <- renderText(sprintf("%.3f", estratificado_reactivo()$eta2))
  output$deff_val  <- renderText(sprintf("%.3f", estratificado_reactivo()$deff))
  
  output$plot_estratos <- renderPlot({
    r <- estratificado_reactivo()
    ggplot(r$estratos, aes(x = reorder(estrato, p_h), y = p_h)) +
      geom_point(color = "#D9534F", size = 3) +
      geom_errorbar(aes(ymin = p_h - 1.96 * sqrt(var_h), ymax = p_h + 1.96 * sqrt(var_h)),
                    width = 0.15, color = "#D9534F") +
      geom_hline(yintercept = r$p_mas, linetype = "dashed", color = "#2C6E63") +
      coord_flip() +
      scale_y_continuous(labels = percent) +
      labs(x = NULL, y = "Proporción vacunados (IC 95%) · línea = estimación global MAS") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
  })
  
  output$tabla_estratos <- renderTable({
    r <- estratificado_reactivo()
    r$estratos %>%
      mutate(across(c(W_h, p_h, var_h, S_h), ~ round(.x, 4))) %>%
      rename(Estrato = estrato) %>%
      bind_rows(data.frame(
        Estrato = "GLOBAL ESTRATIFICADO", n_h = r$n_total, W_h = NA,
        p_h = round(r$p_estratificado, 4), var_h = round(r$var_estratificada, 6), S_h = NA
      ))
  })
  
  output$plot_asignacion <- renderPlot({
    r <- estratificado_reactivo()
    df <- r$estratos_asig %>%
      select(estrato, n_igual, n_proporcional, n_neyman) %>%
      pivot_longer(-estrato, names_to = "tipo", values_to = "n")
    
    ggplot(df, aes(x = estrato, y = n, fill = tipo)) +
      geom_col(position = "dodge", width = 0.65) +
      scale_fill_manual(values = c(n_igual = "#8FC1B5", n_proporcional = "#2C6E63", n_neyman = "#D9534F"),
                        labels = c(n_igual = "Igualitaria", n_proporcional = "Proporcional", n_neyman = "Neyman")) +
      labs(x = NULL, y = "Observaciones asignadas", fill = NULL) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            panel.grid.minor = element_blank())
  })
  
  output$tabla_varianza_asignacion <- renderTable({
    estratificado_reactivo()$comparacion_varianzas
  }, digits = 8)
}

# --------------------------------------------------------------
# 4. Run
# --------------------------------------------------------------
shinyApp(ui = ui, server = server)