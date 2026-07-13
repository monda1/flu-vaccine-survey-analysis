#
# Shiny App: Estimación por Máxima Verosimilitud de la proporción de vacunados
# Datos: National 2009 H1N1 Flu Survey (DrivenData - Flu Shot Learning)
#
# Ejecuta con el botón 'Run App'.
#

library(shiny)
library(dplyr)
library(readr)
library(bslib)
library(ggplot2)

# --------------------------------------------------------------
# 1. Carga y preparación de datos
# --------------------------------------------------------------
# Busca data/raw/ tanto si app.R está en la raíz del proyecto como si
# está dentro de una subcarpeta (p. ej. flu-shot-learning/app/app.R).
encontrar_ruta_datos <- function(archivo) {
  candidatos <- c(
    file.path("data", "raw", archivo),        # app.R en la raíz del proyecto
    file.path("..", "data", "raw", archivo)    # app.R en una subcarpeta (p. ej. /app)
  )
  ruta <- candidatos[file.exists(candidatos)]
  if (length(ruta) == 0) {
    stop(
      "No se encontró '", archivo, "'. Colócalo en 'data/raw/' en la raíz ",
      "de tu proyecto, o ajusta 'encontrar_ruta_datos()' con la ruta correcta.\n",
      "Directorio de trabajo actual: ", getwd()
    )
  }
  ruta[1]
}

features <- read_csv(encontrar_ruta_datos("training_set_features.csv"), show_col_types = FALSE)
labels   <- read_csv(encontrar_ruta_datos("training_set_labels.csv"),   show_col_types = FALSE)

datos <- features %>%
  left_join(labels, by = "respondent_id")

# --------------------------------------------------------------
# 2. Funciones estadísticas (reutilizadas de tu script original)
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

log_verosimilitud <- function(p, x) sum(dbinom(x, size = 1, prob = p, log = TRUE))

fisher_info <- function(p_hat, n) n / (p_hat * (1 - p_hat))

# Variables categóricas/binarias candidatas para explorar tasas de vacunación
vars_candidatas <- c(
  "h1n1_concern", "h1n1_knowledge", "doctor_recc_h1n1", "doctor_recc_seasonal",
  "chronic_med_condition", "child_under_6_months", "health_worker",
  "health_insurance", "age_group", "education", "race", "sex",
  "income_poverty", "marital_status", "rent_or_own", "employment_status",
  "census_msa"
)
vars_disponibles <- intersect(vars_candidatas, names(datos))

# --------------------------------------------------------------
# 3. UI
# --------------------------------------------------------------
ui <- page_sidebar(
  title = "Estimación de Máxima Verosimilitud · Vacunación contra la gripe",
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
    sliderInput(
      "n_sim", "Simulaciones Monte Carlo",
      min = 500, max = 10000, value = 3000, step = 500
    ),
    sliderInput(
      "n_muestra_consistencia", "Tamaños de muestra a comparar (máx.)",
      min = 100, max = 5000, value = 2000, step = 100
    ),
    hr(),
    selectInput(
      "variable_explorar", "Explorar tasa de vacunación por:",
      choices = vars_disponibles,
      selected = if ("age_group" %in% vars_disponibles) "age_group" else vars_disponibles[1]
    ),
    hr(),
    downloadButton("descargar_resumen", "Descargar resumen (CSV)", class = "btn-sm w-100"),
    hr(),
    p(class = "text-muted small",
      "La app estima p = P(vacunado) por Máxima Verosimilitud y
       comprueba sus propiedades: insesgadez, eficiencia (Cramér-Rao),
       suficiencia y consistencia.")
  ),
  
  layout_columns(
    col_widths = c(4, 4, 4),
    value_box(
      title = "p̂ (MLE)",
      value = textOutput("p_hat_val"),
      showcase = bsicons::bs_icon("bullseye"),
      theme = "primary"
    ),
    value_box(
      title = "Error estándar",
      value = textOutput("ee_val"),
      showcase = bsicons::bs_icon("rulers"),
      theme = "secondary"
    ),
    value_box(
      title = "IC 95%",
      value = textOutput("ic_val"),
      showcase = bsicons::bs_icon("arrows-angle-expand"),
      theme = "success"
    )
  ),
  
  navset_card_underline(
    nav_panel(
      "Log-verosimilitud",
      card_body(
        plotOutput("plot_loglik", height = "380px"),
        p(class = "text-muted small mt-2",
          "El máximo de l(p) se alcanza exactamente en p̂ = x̄, la proporción muestral.")
      )
    ),
    nav_panel(
      "Monte Carlo (insesgadez y ECM)",
      card_body(
        plotOutput("plot_mc", height = "380px"),
        tableOutput("tabla_mc")
      )
    ),
    nav_panel(
      "Consistencia",
      card_body(
        plotOutput("plot_consistencia", height = "380px"),
        p(class = "text-muted small mt-2",
          "A medida que n crece, la varianza de p̂ se reduce y se concentra en torno a p.")
      )
    ),
    nav_panel(
      "Eficiencia (Cramér-Rao)",
      card_body(
        tableOutput("tabla_eficiencia"),
        p(class = "text-muted small mt-2",
          "Si Var(p̂) coincide con la Cota de Cramér-Rao, el estimador es eficiente
           (varianza mínima posible entre estimadores insesgados).")
      )
    ),
    nav_panel(
      "Explorar por grupo",
      card_body(
        plotOutput("plot_grupo", height = "380px"),
        p(class = "text-muted small mt-2",
          "Proporción estimada de vacunados (p̂ por grupo) con intervalo de confianza al 95%,
           calculada de forma independiente para cada categoría.")
      )
    )
  )
)

# --------------------------------------------------------------
# 4. Server
# --------------------------------------------------------------
server <- function(input, output, session) {
  
  mle_reactivo <- reactive({
    calcular_mle(datos[[input$vacuna]])
  })
  
  x_reactivo <- reactive({
    v <- datos[[input$vacuna]]
    v[!is.na(v)]
  })
  
  # ---- Value boxes ----
  output$p_hat_val <- renderText(sprintf("%.3f", mle_reactivo()$p_hat))
  output$ee_val     <- renderText(sprintf("%.4f", mle_reactivo()$ee))
  output$ic_val      <- renderText({
    ic <- mle_reactivo()$ic_95
    sprintf("[%.3f, %.3f]", ic[1], ic[2])
  })
  
  # ---- Log-verosimilitud ----
  output$plot_loglik <- renderPlot({
    mle <- mle_reactivo()
    x <- x_reactivo()
    p_seq <- seq(0.01, 0.99, length.out = 300)
    ll <- sapply(p_seq, log_verosimilitud, x = x)
    df <- data.frame(p = p_seq, ll = ll)
    
    ggplot(df, aes(p, ll)) +
      geom_line(color = "#2C6E63", linewidth = 1) +
      geom_vline(xintercept = mle$p_hat, color = "#D9534F", linetype = "dashed") +
      annotate("text", x = mle$p_hat, y = min(ll),
               label = sprintf("  p̂ = %.3f", mle$p_hat),
               hjust = 0, color = "#D9534F") +
      labs(x = "p", y = "log-verosimilitud  l(p)") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
  })
  
  # ---- Monte Carlo ----
  mc_reactivo <- reactive({
    mle <- mle_reactivo()
    set.seed(123)
    p_hats <- replicate(input$n_sim, {
      muestra <- rbinom(mle$n, size = 1, prob = mle$p_hat)
      mean(muestra)
    })
    list(p_hats = p_hats, mle = mle)
  })
  
  output$plot_mc <- renderPlot({
    res <- mc_reactivo()
    df <- data.frame(p_hat = res$p_hats)
    
    ggplot(df, aes(p_hat)) +
      geom_histogram(aes(y = after_stat(density)), bins = 45,
                     fill = "#8FC1B5", color = "white") +
      stat_function(fun = dnorm,
                    args = list(mean = res$mle$p_hat, sd = sqrt(res$mle$var)),
                    color = "#D9534F", linewidth = 1) +
      geom_vline(xintercept = res$mle$p_hat, color = "#D9534F", linetype = "dashed") +
      labs(x = expression(hat(p)), y = "densidad") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
  })
  
  output$tabla_mc <- renderTable({
    res <- mc_reactivo()
    sesgo <- mean(res$p_hats) - res$mle$p_hat
    ecm <- mean((res$p_hats - res$mle$p_hat)^2)
    data.frame(
      Métrica = c("Sesgo empírico", "Varianza empírica", "ECM empírico",
                  "Varianza teórica", "ECM teórico"),
      Valor = c(sesgo, var(res$p_hats), ecm, res$mle$var, res$mle$var)
    )
  }, digits = 6)
  
  # ---- Consistencia ----
  output$plot_consistencia <- renderPlot({
    mle <- mle_reactivo()
    n_max <- input$n_muestra_consistencia
    n_seq <- unique(pmin(c(10, 50, 100, 500, 1000, n_max), mle$n))
    n_seq <- sort(n_seq)
    
    set.seed(7)
    df <- bind_rows(lapply(n_seq, function(n_i) {
      p_hats <- replicate(1000, mean(rbinom(n_i, 1, mle$p_hat)))
      data.frame(n = factor(n_i, levels = n_seq), p_hat = p_hats)
    }))
    
    ggplot(df, aes(n, p_hat)) +
      geom_boxplot(fill = "#8FC1B5", outlier.alpha = 0.3) +
      geom_hline(yintercept = mle$p_hat, color = "#D9534F", linetype = "dashed") +
      labs(x = "Tamaño muestral (n)", y = expression(hat(p))) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
  })
  
  # ---- Eficiencia ----
  output$tabla_eficiencia <- renderTable({
    mle <- mle_reactivo()
    I_n <- fisher_info(mle$p_hat, mle$n)
    crlb <- 1 / I_n
    data.frame(
      Métrica = c("Información de Fisher I_n(p)", "Cota de Cramér-Rao",
                  "Var(p̂) observada", "¿Es eficiente?"),
      Valor = c(as.character(round(I_n, 2)),
                as.character(round(crlb, 6)),
                as.character(round(mle$var, 6)),
                as.character(isTRUE(all.equal(crlb, mle$var))))
    )
  })
  # ---- Explorar por grupo ----
  output$plot_grupo <- renderPlot({
    var <- input$variable_explorar
    df <- datos %>%
      filter(!is.na(.data[[var]]), !is.na(.data[[input$vacuna]])) %>%
      group_by(grupo = as.factor(.data[[var]])) %>%
      summarise(
        n = n(),
        p_hat = mean(.data[[input$vacuna]]),
        ee = sqrt(p_hat * (1 - p_hat) / n),
        .groups = "drop"
      ) %>%
      mutate(
        ic_inf = pmax(0, p_hat - qnorm(0.975) * ee),
        ic_sup = pmin(1, p_hat + qnorm(0.975) * ee)
      )
    
    ggplot(df, aes(x = reorder(grupo, p_hat), y = p_hat)) +
      geom_col(fill = "#8FC1B5", width = 0.65) +
      geom_errorbar(aes(ymin = ic_inf, ymax = ic_sup), width = 0.2, color = "#2C6E63") +
      coord_flip() +
      labs(x = NULL, y = "Proporción vacunados (p̂)") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank())
  })
  
  # ---- Descarga del resumen ----
  output$descargar_resumen <- downloadHandler(
    filename = function() paste0("resumen_mle_", input$vacuna, ".csv"),
    content = function(file) {
      mle <- mle_reactivo()
      var <- input$variable_explorar
      resumen_grupo <- datos %>%
        filter(!is.na(.data[[var]]), !is.na(.data[[input$vacuna]])) %>%
        group_by(grupo = .data[[var]]) %>%
        summarise(n = n(), p_hat = mean(.data[[input$vacuna]]), .groups = "drop")
      
      resumen_general <- data.frame(
        vacuna = input$vacuna,
        n = mle$n,
        p_hat_MLE = mle$p_hat,
        ee = mle$ee,
        ic_95_inf = mle$ic_95[1],
        ic_95_sup = mle$ic_95[2]
      )
      
      writeLines("# Resumen general", file)
      write.table(resumen_general, file, sep = ",", row.names = FALSE, append = TRUE)
      cat("\n# Resumen por grupo:", var, "\n", file = file, append = TRUE)
      write.table(resumen_grupo, file, sep = ",", row.names = FALSE, append = TRUE)
    }
  )
  
}

# --------------------------------------------------------------
# 5. Run
# --------------------------------------------------------------
shinyApp(ui = ui, server = server)