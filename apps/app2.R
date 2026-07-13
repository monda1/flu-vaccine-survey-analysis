#
# Shiny App: Contraste de proporciones (health_worker vs resto)
# e Intervalos de confianza por grupo
# Datos: National 2009 H1N1 Flu Survey (DrivenData - Flu Shot Learning)
#
# Ejecuta con el botón 'Run App'.
#

library(shiny)
library(dplyr)
library(readr)
library(bslib)
library(bsicons)
library(ggplot2)

# --------------------------------------------------------------
# 1. Carga y preparación de datos (ruta robusta, igual que app.R)
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

# --------------------------------------------------------------
# 2. Funciones estadísticas (de tu script original)
# --------------------------------------------------------------
comprobar_condiciones <- function(datos, var_vacuna) {
  df <- datos %>% filter(!is.na(.data[[var_vacuna]]), !is.na(health_worker))
  tabla <- table(df$health_worker, df[[var_vacuna]])
  
  n1 <- sum(tabla["1", ]); n0 <- sum(tabla["0", ])
  p1 <- tabla["1", "1"] / n1; p0 <- tabla["0", "1"] / n0
  
  cond <- c(n1 * p1, n1 * (1 - p1), n0 * p0, n0 * (1 - p0))
  ok <- all(cond >= 5)
  
  list(ok = ok, valores = cond)
}

contraste_proporciones <- function(datos, var_vacuna) {
  df <- datos %>% filter(!is.na(.data[[var_vacuna]]), !is.na(health_worker))
  
  tabla <- table(df$health_worker, df[[var_vacuna]])
  x_vacunados <- c(tabla["1", "1"], tabla["0", "1"])
  n_total <- c(sum(tabla["1", ]), sum(tabla["0", ]))
  
  test_2c <- prop.test(x = x_vacunados, n = n_total, alternative = "two.sided", correct = TRUE)
  test_1c <- prop.test(x = x_vacunados, n = n_total, alternative = "greater", correct = TRUE)
  test_fisher <- fisher.test(matrix(c(x_vacunados[1], n_total[1] - x_vacunados[1],
                                      x_vacunados[2], n_total[2] - x_vacunados[2]),
                                    nrow = 2, byrow = TRUE))
  
  p1 <- x_vacunados[1] / n_total[1]
  p2 <- x_vacunados[2] / n_total[2]
  
  h_cohen <- 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))
  or_val <- (x_vacunados[1] * (n_total[2] - x_vacunados[2])) /
    (x_vacunados[2] * (n_total[1] - x_vacunados[1]))
  
  list(
    p1 = p1, p2 = p2, n1 = n_total[1], n2 = n_total[2],
    diff = p1 - p2, h_cohen = h_cohen, or = or_val,
    chi2_stat = unname(test_2c$statistic), chi2_p = test_2c$p.value,
    chi2_ic = test_2c$conf.int, p_1c = test_1c$p.value,
    fisher_or = unname(test_fisher$estimate), fisher_p = test_fisher$p.value,
    condiciones = comprobar_condiciones(datos, var_vacuna)
  )
}

ic_por_grupo <- function(datos, var_vacuna, var_grupo) {
  datos %>%
    filter(!is.na(.data[[var_vacuna]]), !is.na(.data[[var_grupo]])) %>%
    group_by(grupo = as.factor(.data[[var_grupo]])) %>%
    summarise(
      n = n(),
      vacunados = sum(.data[[var_vacuna]]),
      p_hat = vacunados / n,
      ee = sqrt(p_hat * (1 - p_hat) / n),
      ic_inf_wald = pmax(p_hat - qnorm(0.975) * ee, 0),
      ic_sup_wald = pmin(p_hat + qnorm(0.975) * ee, 1),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      z = qnorm(0.975),
      denom = 1 + z^2 / n,
      centro = (p_hat + z^2 / (2 * n)) / denom,
      margen = (z / denom) * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2)),
      ic_inf_wilson = centro - margen,
      ic_sup_wilson = centro + margen
    ) %>%
    ungroup() %>%
    mutate(amplitud_wald = ic_sup_wald - ic_inf_wald,
           amplitud_wilson = ic_sup_wilson - ic_inf_wilson) %>%
    select(grupo, n, vacunados, p_hat,
           ic_inf_wald, ic_sup_wald, amplitud_wald,
           ic_inf_wilson, ic_sup_wilson, amplitud_wilson)
}

vars_grupo_candidatas <- c("age_group", "education", "race", "sex", "income_poverty",
                           "marital_status", "rent_or_own", "employment_status", "census_msa")
vars_grupo <- intersect(vars_grupo_candidatas, names(datos))

# --------------------------------------------------------------
# 3. UI
# --------------------------------------------------------------
ui <- page_sidebar(
  title = "Contraste de proporciones · Sanitarios vs resto",
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
      "vacuna", "Vacuna",
      choices = c("H1N1" = "h1n1_vaccine", "Estacional" = "seasonal_vaccine"),
      selected = "h1n1_vaccine"
    ),
    hr(),
    selectInput(
      "var_grupo", "Desagregar intervalos por:",
      choices = vars_grupo,
      selected = if ("age_group" %in% vars_grupo) "age_group" else vars_grupo[1]
    ),
    radioButtons(
      "tipo_ic", "Tipo de intervalo",
      choices = c("Wilson (recomendado)" = "wilson", "Wald (clásico)" = "wald"),
      selected = "wilson"
    ),
    hr(),
    downloadButton("descargar_ic", "Descargar IC por grupo (CSV)", class = "btn-sm w-100"),
    hr(),
    p(class = "text-muted small",
      "Contrasta si la proporción de vacunados difiere entre sanitarios y
       no sanitarios (chi-cuadrado, test exacto de Fisher, tamaño del
       efecto) y calcula intervalos de confianza de la tasa de vacunación
       por grupo.")
  ),
  
  layout_columns(
    col_widths = c(3, 3, 3, 3),
    value_box(
      title = "p̂ sanitarios",
      value = textOutput("p1_val"),
      showcase = bs_icon("heart-pulse"),
      theme = "primary"
    ),
    value_box(
      title = "p̂ no sanitarios",
      value = textOutput("p2_val"),
      showcase = bs_icon("people"),
      theme = "secondary"
    ),
    value_box(
      title = "Diferencia",
      value = textOutput("diff_val"),
      showcase = bs_icon("arrow-left-right"),
      theme = "success"
    ),
    value_box(
      title = "p-valor (χ²)",
      value = textOutput("pval_val"),
      showcase = bs_icon("check2-circle"),
      theme = "warning"
    )
  ),
  
  navset_card_underline(
    nav_panel(
      "Contraste sanitarios vs resto",
      card_body(
        plotOutput("plot_hw", height = "320px"),
        tableOutput("tabla_tests"),
        p(class = "text-muted small mt-2",
          "H0: misma proporción de vacunados en ambos grupos. Se reportan el test
           chi-cuadrado (dos colas), su versión unilateral (sanitarios > resto,
           justificada por mayor exposición y acceso a la vacuna) y el test exacto
           de Fisher. Umbral de Bonferroni para 2 vacunas comparadas: α = 0.025.")
      )
    ),
    nav_panel(
      "Condiciones del test",
      card_body(
        tableOutput("tabla_condiciones"),
        p(class = "text-muted small mt-2",
          "La aproximación chi-cuadrado/normal requiere n·p ≥ 5 y n·(1-p) ≥ 5 en
           cada grupo. Si no se cumple, el test exacto de Fisher (pestaña anterior)
           es la referencia fiable.")
      )
    ),
    nav_panel(
      "IC por grupo",
      card_body(
        plotOutput("plot_ic_grupo", height = "380px"),
        tableOutput("tabla_ic_grupo"),
        p(class = "text-muted small mt-2",
          "Comparar varios grupos a la vez infla el error de tipo I; para
           comparaciones múltiples, el nivel de confianza ajustado por Bonferroni
           se muestra debajo de la tabla.")
      )
    )
  )
)

# --------------------------------------------------------------
# 4. Server
# --------------------------------------------------------------
server <- function(input, output, session) {
  
  test_reactivo <- reactive({
    contraste_proporciones(datos, input$vacuna)
  })
  
  ic_reactivo <- reactive({
    ic_por_grupo(datos, input$vacuna, input$var_grupo)
  })
  
  # ---- Value boxes ----
  output$p1_val   <- renderText(sprintf("%.3f", test_reactivo()$p1))
  output$p2_val   <- renderText(sprintf("%.3f", test_reactivo()$p2))
  output$diff_val <- renderText(sprintf("%+.3f", test_reactivo()$diff))
  output$pval_val <- renderText(format.pval(test_reactivo()$chi2_p, digits = 3))
  
  # ---- Gráfico sanitarios vs resto ----
  output$plot_hw <- renderPlot({
    df <- datos %>%
      filter(!is.na(health_worker), !is.na(.data[[input$vacuna]])) %>%
      mutate(grupo_hw = ifelse(health_worker == 1, "Sanitario", "No sanitario")) %>%
      group_by(grupo_hw) %>%
      summarise(n = n(), p_hat = mean(.data[[input$vacuna]]),
                ee = sqrt(p_hat * (1 - p_hat) / n), .groups = "drop") %>%
      mutate(ic_inf = pmax(p_hat - qnorm(0.975) * ee, 0),
             ic_sup = pmin(p_hat + qnorm(0.975) * ee, 1))
    
    ggplot(df, aes(x = grupo_hw, y = p_hat, fill = grupo_hw)) +
      geom_col(width = 0.55, show.legend = FALSE) +
      geom_errorbar(aes(ymin = ic_inf, ymax = ic_sup), width = 0.15, color = "#2C6E63") +
      scale_fill_manual(values = c("Sanitario" = "#2C6E63", "No sanitario" = "#8FC1B5")) +
      labs(x = NULL, y = "Proporción vacunados (p̂)") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank())
  })
  
  output$tabla_tests <- renderTable({
    t <- test_reactivo()
    data.frame(
      Métrica = c("n sanitarios", "n no sanitarios", "Tamaño del efecto (Cohen's h)",
                  "Odds ratio", "Estadístico χ²", "p-valor χ² (dos colas)",
                  "IC 95% diferencia (χ²)", "p-valor (una cola, sanitarios > resto)",
                  "Odds ratio (Fisher)", "p-valor (Fisher, exacto)"),
      Valor = c(
        as.character(t$n1), as.character(t$n2),
        sprintf("%.4f", t$h_cohen), sprintf("%.4f", t$or),
        sprintf("%.4f", t$chi2_stat), format.pval(t$chi2_p, digits = 4),
        sprintf("[%.4f, %.4f]", t$chi2_ic[1], t$chi2_ic[2]),
        format.pval(t$p_1c, digits = 4),
        sprintf("%.4f", t$fisher_or), format.pval(t$fisher_p, digits = 4)
      )
    )
  })
  
  # ---- Condiciones de aplicación ----
  output$tabla_condiciones <- renderTable({
    t <- test_reactivo()
    data.frame(
      Grupo = c("Sanitarios (éxitos)", "Sanitarios (fracasos)",
                "No sanitarios (éxitos)", "No sanitarios (fracasos)"),
      `n·p o n·(1-p)` = round(t$condiciones$valores, 1),
      `¿>= 5?` = t$condiciones$valores >= 5,
      check.names = FALSE
    )
  })
  
  # ---- IC por grupo ----
  output$plot_ic_grupo <- renderPlot({
    df <- ic_reactivo()
    if (input$tipo_ic == "wilson") {
      df <- df %>% rename(ic_inf = ic_inf_wilson, ic_sup = ic_sup_wilson)
    } else {
      df <- df %>% rename(ic_inf = ic_inf_wald, ic_sup = ic_sup_wald)
    }
    
    ggplot(df, aes(x = reorder(grupo, p_hat), y = p_hat)) +
      geom_point(color = "#2C6E63", size = 3) +
      geom_errorbar(aes(ymin = ic_inf, ymax = ic_sup), width = 0.2, color = "#2C6E63") +
      coord_flip() +
      labs(x = NULL, y = "Proporción vacunados (p̂)") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
  })
  
  output$tabla_ic_grupo <- renderTable({
    ic_reactivo() %>%
      mutate(across(where(is.numeric), ~ round(.x, 4)))
  })
  
  output$descargar_ic <- downloadHandler(
    filename = function() paste0("ic_", input$var_grupo, "_", input$vacuna, ".csv"),
    content = function(file) {
      write.csv(ic_reactivo(), file, row.names = FALSE)
    }
  )
}

# --------------------------------------------------------------
# 5. Run
# --------------------------------------------------------------
shinyApp(ui = ui, server = server)