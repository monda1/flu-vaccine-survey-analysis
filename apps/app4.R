#
# Shiny App: (1) Distribución conjunta, Bernoulli, IC y estimación
# bayesiana de p  ·  (2) Regresión logística múltiple, bondad de
# ajuste y multicolinealidad (VIF)
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
library(broom)

# Si el paquete 'conflicted' está cargado en la sesión (p. ej. por otro
# script previamente ejecutado), algunas funciones base como update()
# quedan ambiguas entre paquetes. Fijamos aquí explícitamente cuál usar.
if ("conflicted" %in% loadedNamespaces()) {
  conflicted::conflicts_prefer(stats::update, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::filter, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::select, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::mutate, .quiet = TRUE)
}

pkg_car   <- requireNamespace("car", quietly = TRUE)
pkg_proc  <- requireNamespace("pROC", quietly = TRUE)
pkg_hl    <- requireNamespace("ResourceSelection", quietly = TRUE)
pkg_glmnet <- requireNamespace("glmnet", quietly = TRUE)
if (pkg_car) library(car)

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

datos$h1n1_vaccine     <- as.numeric(as.character(datos$h1n1_vaccine))
datos$seasonal_vaccine <- as.numeric(as.character(datos$seasonal_vaccine))
datos$total_vacunas    <- datos$h1n1_vaccine + datos$seasonal_vaccine

datos_biv <- datos %>% filter(!is.na(h1n1_vaccine), !is.na(seasonal_vaccine))
tabla_conjunta <- table(H1N1 = datos_biv$h1n1_vaccine, Estacional = datos_biv$seasonal_vaccine)
prop_conjunta  <- prop.table(tabla_conjunta)

p_h1n1     <- mean(datos_biv$h1n1_vaccine)
p_seasonal <- mean(datos_biv$seasonal_vaccine)
covarianza  <- cov(datos_biv$h1n1_vaccine, datos_biv$seasonal_vaccine)
correlacion <- cor(datos_biv$h1n1_vaccine, datos_biv$seasonal_vaccine)
test_chi_biv <- chisq.test(tabla_conjunta)

# --------------------------------------------------------------
# 2. Modelo de regresión logística (se ajusta una sola vez al arrancar)
# --------------------------------------------------------------
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
    age_group = factor(age_group), education = factor(education),
    income_poverty = factor(income_poverty), sex = factor(sex),
    h1n1_vaccine = factor(h1n1_vaccine)
  )

modelo_logit <- glm(
  h1n1_vaccine ~ h1n1_concern + h1n1_knowledge +
    opinion_h1n1_vacc_effective + opinion_h1n1_risk + opinion_h1n1_sick_from_vacc +
    opinion_seas_vacc_effective + opinion_seas_risk + opinion_seas_sick_from_vacc +
    doctor_recc_h1n1 + doctor_recc_seasonal +
    chronic_med_condition + health_worker + health_insurance +
    age_group + education + income_poverty + sex,
  data = datos_modelo, family = binomial(link = "logit")
)

odds_ratios <- tidy(modelo_logit, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  arrange(desc(abs(estimate - 1)))

pseudo_r2 <- 1 - (modelo_logit$deviance / modelo_logit$null.deviance)
pred_prob <- predict(modelo_logit, type = "response")

calcular_ame_numerica <- function(modelo, datos_, variable, delta = 1) {
  datos_mas <- datos_
  datos_mas[[variable]] <- datos_mas[[variable]] + delta
  p_base <- predict(modelo, newdata = datos_, type = "response")
  p_mas  <- predict(modelo, newdata = datos_mas, type = "response")
  mean(p_mas - p_base, na.rm = TRUE)
}
vars_ame <- c("h1n1_concern", "h1n1_knowledge", "opinion_h1n1_risk",
              "opinion_h1n1_vacc_effective", "opinion_seas_risk")
tabla_ame <- data.frame(
  variable = vars_ame,
  ame_pp = sapply(vars_ame, function(v) calcular_ame_numerica(modelo_logit, datos_modelo, v) * 100)
) %>% arrange(desc(abs(ame_pp)))

if (pkg_car) {
  vif_resultado <- vif(modelo_logit)
  if (is.matrix(vif_resultado)) {
    vif_df <- data.frame(variable = rownames(vif_resultado),
                         VIF_ajustado = vif_resultado[, "GVIF^(1/(2*Df))"]^2)
  } else {
    vif_df <- data.frame(variable = names(vif_resultado), VIF_ajustado = vif_resultado)
  }
  vif_df <- vif_df %>% arrange(desc(VIF_ajustado))
}

vars_numericas_modelo <- c("h1n1_concern", "h1n1_knowledge",
                           "opinion_h1n1_vacc_effective", "opinion_h1n1_risk", "opinion_h1n1_sick_from_vacc",
                           "opinion_seas_vacc_effective", "opinion_seas_risk", "opinion_seas_sick_from_vacc")
cor_matrix_num <- cor(datos_modelo[vars_numericas_modelo])
autovalores <- eigen(cor_matrix_num)$values
indice_condicion <- sqrt(max(autovalores) / autovalores)

modelo_interaccion <- glm(
  h1n1_vaccine ~ doctor_recc_h1n1 * health_worker +
    h1n1_concern + chronic_med_condition + health_insurance + age_group,
  data = datos_modelo, family = binomial(link = "logit")
)
modelo_sin_interaccion <- stats::update(modelo_interaccion, . ~ . - doctor_recc_h1n1:health_worker)
lr_test_interaccion <- anova(modelo_sin_interaccion, modelo_interaccion, test = "Chisq")

modelo_sin_opiniones <- stats::update(
  modelo_logit,
  . ~ . - opinion_h1n1_vacc_effective - opinion_h1n1_risk - opinion_h1n1_sick_from_vacc -
    opinion_seas_vacc_effective - opinion_seas_risk - opinion_seas_sick_from_vacc
)
lr_test_opiniones <- anova(modelo_sin_opiniones, modelo_logit, test = "Chisq")

if (pkg_proc) {
  library(pROC)
  roc_obj <- roc(datos_modelo$h1n1_vaccine, pred_prob, quiet = TRUE)
}
if (pkg_hl) {
  library(ResourceSelection)
  hl_test <- hoslem.test(as.numeric(as.character(datos_modelo$h1n1_vaccine)), pred_prob, g = 10)
}
if (pkg_glmnet) {
  library(glmnet)
  x_ridge <- model.matrix(modelo_logit)[, -1]
  y_ridge <- datos_modelo$h1n1_vaccine
  set.seed(123)
  cv_ridge <- cv.glmnet(x_ridge, y_ridge, family = "binomial", alpha = 0)
  coef_ridge <- as.matrix(coef(cv_ridge, s = "lambda.min"))
  coef_glm <- coef(modelo_logit)
  comparacion_coef <- data.frame(
    variable = rownames(coef_ridge),
    coef_glm_original = coef_glm[rownames(coef_ridge)],
    coef_ridge = coef_ridge[, 1]
  ) %>% filter(!is.na(coef_glm_original))
}

# --------------------------------------------------------------
# 3. UI
# --------------------------------------------------------------
ui <- page_sidebar(
  title = "Distribución conjunta & Regresión logística · Vacunación H1N1",
  theme = bs_theme(
    version = 5, bootswatch = "minty",
    base_font = font_google("Inter"), heading_font = font_google("Inter"),
    primary = "#2C6E63"
  ),
  
  sidebar = sidebar(
    width = 300,
    h6("Distribución conjunta y estimación"),
    sliderInput("nivel_conf", "Nivel de confianza (IC)", min = 0.80, max = 0.99, value = 0.95, step = 0.01),
    sliderInput("alpha_prior", "Prior Beta: α", min = 0.5, max = 10, value = 1, step = 0.5),
    sliderInput("beta_prior", "Prior Beta: β", min = 0.5, max = 10, value = 1, step = 0.5),
    hr(),
    h6("Regresión logística"),
    selectInput("var_ame", "Variable para efecto marginal (AME):",
                choices = vars_ame, selected = vars_ame[1]),
    hr(),
    p(class = "text-muted small",
      "El modelo logístico y sus diagnósticos se ajustan una sola vez al iniciar
       la app (no dependen de los controles de arriba, salvo la selección de AME).")
  ),
  
  navset_card_underline(
    nav_panel(
      "Distribución conjunta",
      card_body(
        plotOutput("plot_conjunta", height = "320px"),
        tableOutput("tabla_marginales"),
        radioButtons("direccion_cond", "Distribución condicionada:",
                     choices = c("P(Estacional | H1N1)" = "h1n1", "P(H1N1 | Estacional)" = "seasonal"),
                     inline = TRUE),
        tableOutput("tabla_condicionada"),
        p(class = "text-muted small mt-2",
          "Mapa de calor de la distribución conjunta P(H1N1, Estacional), con las
           marginales de cada vacuna y la distribución condicionada en la dirección
           elegida.")
      )
    ),
    nav_panel(
      "Covarianza e independencia",
      card_body(
        tableOutput("tabla_covar"),
        p(class = "text-muted small mt-2",
          "Covarianza, correlación (coeficiente φ para dos binarias), test χ² de
           independencia y comprobación directa de la definición P(A∩B) = P(A)·P(B).")
      )
    ),
    nav_panel(
      "Bernoulli: momentos",
      card_body(
        tableOutput("tabla_bernoulli"),
        p(class = "text-muted small mt-2",
          "Para una Bernoulli, E[X]=p, Var(X)=p(1-p), y E[Xᵏ]=p para cualquier k≥1
           (porque 0ᵏ=0 y 1ᵏ=1). La tabla comprueba esto último empíricamente.")
      )
    ),
    nav_panel(
      "Intervalos de confianza",
      card_body(
        tableOutput("tabla_ic"),
        p(class = "text-muted small mt-2",
          "Comparación del IC para p con la aproximación normal (TCL) frente al
           intervalo binomial exacto de Clopper-Pearson, al nivel de confianza
           elegido en el sidebar.")
      )
    ),
    nav_panel(
      "Estimación bayesiana",
      card_body(
        plotOutput("plot_bayes", height = "320px"),
        tableOutput("tabla_bayes"),
        p(class = "text-muted small mt-2",
          "Posterior Beta(α + éxitos, β + fracasos) tras observar los datos, partiendo
           del prior Beta(α, β) elegido en el sidebar. Beta(1,1) equivale a un prior
           uniforme (no informativo).")
      )
    ),
    nav_panel(
      "Z = X + Y",
      card_body(
        plotOutput("plot_Z", height = "300px"),
        tableOutput("tabla_Z"),
        p(class = "text-muted small mt-2",
          "Z = nº total de vacunas recibidas (0, 1 o 2). E[Z]=E[X]+E[Y] siempre se
           cumple; Var[Z]=Var(X)+Var(Y)+2·Cov(X,Y) depende de la covarianza entre
           ambas vacunas.")
      )
    ),
    nav_panel(
      "Odds ratios y AME",
      card_body(
        plotOutput("plot_odds", height = "380px"),
        tableOutput("tabla_ame"),
        p(class = "text-muted small mt-2",
          "Odds ratios del modelo logístico completo (excluye el intercepto), con IC
           95%; la línea vertical marca OR=1 (sin efecto). Abajo, el efecto marginal
           promedio (AME) de la variable elegida, en puntos de probabilidad reales de
           vacunarse — más comparable entre variables que el odds ratio.")
      )
    ),
    nav_panel(
      "Bondad de ajuste",
      card_body(
        uiOutput("bloque_roc"),
        tableOutput("tabla_ajuste"),
        p(class = "text-muted small mt-2",
          "Pseudo R² de McFadden (no comparable a un R² lineal, valores 0.1-0.3 son
           normales en logística), AUC de la curva ROC y test de Hosmer-Lemeshow
           (p > 0.05 = no hay evidencia de falta de ajuste).")
      )
    ),
    nav_panel(
      "Multicolinealidad (VIF)",
      card_body(
        plotOutput("plot_vif", height = "380px"),
        tableOutput("tabla_indice_condicion"),
        p(class = "text-muted small mt-2",
          "VIF ajustado por predictor (líneas de referencia en 5 y 10). El índice de
           condición, basado en los autovalores de la matriz de correlación de las
           variables de opinión, detecta colinealidad estructural que el VIF por sí
           solo puede pasar por alto (regla: >30 indica colinealidad grave).")
      )
    ),
    nav_panel(
      "Interacción y modelos anidados",
      card_body(
        h6("¿Depende el efecto de la recomendación médica de ser sanitario?"),
        tableOutput("tabla_interaccion"),
        hr(),
        h6("¿Aportan las 6 variables de opinión, como grupo, al modelo?"),
        tableOutput("tabla_opiniones"),
        p(class = "text-muted small mt-2",
          "Ambos son tests de razón de verosimilitud entre modelos anidados de
           verdad (mismo set de datos, uno es subconjunto exacto de variables del
           otro). p-valor bajo = el término sí aporta información significativa.")
      )
    ),
    nav_panel(
      "Ridge",
      card_body(
        uiOutput("bloque_ridge"),
        p(class = "text-muted small mt-2",
          "La regresión Ridge penaliza (encoge) los coeficientes hacia 0 en lugar de
           eliminar variables, estabilizando las estimaciones cuando hay
           multicolinealidad — el encogimiento es mayor en las variables más
           correlacionadas entre sí (las de opinión).")
      )
    )
  )
)

# --------------------------------------------------------------
# 4. Server
# --------------------------------------------------------------
server <- function(input, output, session) {
  
  # ---- Distribución conjunta ----
  output$plot_conjunta <- renderPlot({
    df <- as.data.frame(prop_conjunta)
    ggplot(df, aes(x = H1N1, y = Estacional, fill = Freq)) +
      geom_tile(color = "white") +
      geom_text(aes(label = scales::percent(Freq, accuracy = 0.1))) +
      scale_fill_gradient(low = "white", high = "#2C6E63") +
      labs(x = "Vacuna H1N1 (0=No, 1=Sí)", y = "Vacuna Estacional (0=No, 1=Sí)", fill = "Proporción") +
      theme_minimal(base_size = 13)
  })
  
  output$tabla_marginales <- renderTable({
    data.frame(
      Vacuna = c("H1N1", "Estacional"),
      `P(No vacunado)` = c(1 - p_h1n1, 1 - p_seasonal),
      `P(Vacunado)` = c(p_h1n1, p_seasonal),
      check.names = FALSE
    )
  }, digits = 4)
  
  output$tabla_condicionada <- renderTable({
    if (input$direccion_cond == "h1n1") {
      m <- prop.table(tabla_conjunta, margin = 1)
      titulo <- "H1N1"
    } else {
      m <- prop.table(tabla_conjunta, margin = 2)
      titulo <- "Estacional"
    }
    df <- as.data.frame.matrix(m)
    df <- cbind(Dado = rownames(df), df)
    df
  }, digits = 4)
  
  # ---- Covarianza e independencia ----
  output$tabla_covar <- renderTable({
    p_11_obs <- prop_conjunta["1", "1"]
    p_11_esp <- p_h1n1 * p_seasonal
    data.frame(
      Métrica = c("Covarianza", "Correlación (φ)", "Estadístico χ²", "p-valor χ²",
                  "P(H1N1=1, Estacional=1) observada", "P(H1N1=1)·P(Estacional=1) esperada",
                  "Diferencia", "¿Independientes?"),
      Valor = c(sprintf("%.4f", covarianza), sprintf("%.4f", correlacion),
                sprintf("%.4f", unname(test_chi_biv$statistic)),
                format.pval(test_chi_biv$p.value, digits = 4),
                sprintf("%.4f", p_11_obs), sprintf("%.4f", p_11_esp),
                sprintf("%.4f", p_11_obs - p_11_esp),
                ifelse(test_chi_biv$p.value > 0.05, "Sí (no se rechaza H0)", "No (se rechaza H0)"))
    )
  })
  
  # ---- Bernoulli ----
  output$tabla_bernoulli <- renderTable({
    var_h1n1 <- p_h1n1 * (1 - p_h1n1)
    var_seasonal <- p_seasonal * (1 - p_seasonal)
    base <- data.frame(
      Métrica = c("E[X] = p", "Var(X) = p(1-p)",
                  paste0("E[X^", 2:4, "]")),
      H1N1 = c(p_h1n1, var_h1n1,
               sapply(2:4, function(k) mean(datos_biv$h1n1_vaccine^k))),
      Estacional = c(p_seasonal, var_seasonal,
                     sapply(2:4, function(k) mean(datos_biv$seasonal_vaccine^k)))
    )
    base %>% mutate(across(c(H1N1, Estacional), ~ round(.x, 4)))
  })
  
  # ---- Intervalos de confianza ----
  output$tabla_ic <- renderTable({
    conf <- input$nivel_conf
    z <- qnorm(1 - (1 - conf) / 2)
    n <- nrow(datos_biv)
    
    calcular_fila <- function(p, x, nombre) {
      se <- sqrt(p * (1 - p) / n)
      ic_normal <- p + c(-1, 1) * z * se
      ic_exacto <- binom.test(x, n, conf.level = conf)$conf.int
      data.frame(Vacuna = nombre, p_hat = round(p, 4),
                 `IC Normal inf` = round(ic_normal[1], 4), `IC Normal sup` = round(ic_normal[2], 4),
                 `IC Exacto inf` = round(ic_exacto[1], 4), `IC Exacto sup` = round(ic_exacto[2], 4),
                 check.names = FALSE)
    }
    bind_rows(
      calcular_fila(p_h1n1, sum(datos_biv$h1n1_vaccine), "H1N1"),
      calcular_fila(p_seasonal, sum(datos_biv$seasonal_vaccine), "Estacional")
    )
  })
  
  # ---- Estimación bayesiana ----
  posterior_reactivo <- reactive({
    a <- input$alpha_prior; b <- input$beta_prior
    n <- nrow(datos_biv)
    post_h1n1 <- c(a + sum(datos_biv$h1n1_vaccine), b + n - sum(datos_biv$h1n1_vaccine))
    post_seasonal <- c(a + sum(datos_biv$seasonal_vaccine), b + n - sum(datos_biv$seasonal_vaccine))
    list(post_h1n1 = post_h1n1, post_seasonal = post_seasonal)
  })
  
  output$plot_bayes <- renderPlot({
    r <- posterior_reactivo()
    curva <- data.frame(p = seq(0, 1, length.out = 500)) %>%
      mutate(H1N1 = dbeta(p, r$post_h1n1[1], r$post_h1n1[2]),
             Estacional = dbeta(p, r$post_seasonal[1], r$post_seasonal[2])) %>%
      tidyr::pivot_longer(-p, names_to = "Vacuna", values_to = "Densidad")
    
    ggplot(curva, aes(p, Densidad, color = Vacuna)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c(H1N1 = "#2C6E63", Estacional = "#D9534F")) +
      labs(x = "p", y = "Densidad posterior") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
  })
  
  output$tabla_bayes <- renderTable({
    r <- posterior_reactivo()
    conf <- input$nivel_conf
    media_h1n1 <- r$post_h1n1[1] / sum(r$post_h1n1)
    media_seasonal <- r$post_seasonal[1] / sum(r$post_seasonal)
    ic_h1n1 <- qbeta(c((1 - conf) / 2, 1 - (1 - conf) / 2), r$post_h1n1[1], r$post_h1n1[2])
    ic_seasonal <- qbeta(c((1 - conf) / 2, 1 - (1 - conf) / 2), r$post_seasonal[1], r$post_seasonal[2])
    data.frame(
      Vacuna = c("H1N1", "Estacional"),
      `Media posterior` = round(c(media_h1n1, media_seasonal), 4),
      `IC creíble inf` = round(c(ic_h1n1[1], ic_seasonal[1]), 4),
      `IC creíble sup` = round(c(ic_h1n1[2], ic_seasonal[2]), 4),
      check.names = FALSE
    )
  })
  
  # ---- Z = X + Y ----
  output$plot_Z <- renderPlot({
    dist_Z <- as.data.frame(prop.table(table(Z = datos_biv$total_vacunas)))
    ggplot(dist_Z, aes(Z, Freq)) +
      geom_col(fill = "#2C6E63", width = 0.5) +
      geom_text(aes(label = scales::percent(Freq, accuracy = 0.1)), vjust = -0.5) +
      labs(x = "Z (nº de vacunas recibidas)", y = "Proporción") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
  })
  
  output$tabla_Z <- renderTable({
    E_Z_teorico <- p_h1n1 + p_seasonal
    E_Z_empirico <- mean(datos_biv$total_vacunas)
    Var_Z_teorico <- p_h1n1 * (1 - p_h1n1) + p_seasonal * (1 - p_seasonal) + 2 * covarianza
    Var_Z_empirico <- var(datos_biv$total_vacunas)
    data.frame(
      Métrica = c("E[Z] teórico", "E[Z] empírico", "Var[Z] teórico", "Var[Z] empírico"),
      Valor = sprintf("%.4f", c(E_Z_teorico, E_Z_empirico, Var_Z_teorico, Var_Z_empirico))
    )
  })
  
  # ---- Odds ratios y AME ----
  output$plot_odds <- renderPlot({
    ggplot(odds_ratios, aes(x = reorder(term, estimate), y = estimate)) +
      geom_point(color = "#2C6E63", size = 2.5) +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, color = "#2C6E63") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "#D9534F") +
      coord_flip() +
      labs(x = NULL, y = "Odds ratio (IC 95%)") +
      theme_minimal(base_size = 11) +
      theme(panel.grid.minor = element_blank())
  })
  
  output$tabla_ame <- renderTable({
    v <- input$var_ame
    ame <- calcular_ame_numerica(modelo_logit, datos_modelo, v)
    data.frame(
      Métrica = c(paste0("AME de ", v, " (p.p. de probabilidad)")),
      Valor = sprintf("%.2f", ame * 100)
    )
  })
  
  # ---- Bondad de ajuste ----
  output$bloque_roc <- renderUI({
    if (pkg_proc) plotOutput("plot_roc", height = "340px")
    else p(class = "text-muted", "Paquete 'pROC' no disponible: instálalo para ver la curva ROC.")
  })
  
  if (pkg_proc) {
    output$plot_roc <- renderPlot({
      ggroc(roc_obj, color = "#D9534F", linewidth = 1) +
        geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "gray50") +
        labs(subtitle = paste0("AUC = ", round(auc(roc_obj), 3)),
             x = "Especificidad", y = "Sensibilidad") +
        theme_minimal(base_size = 13)
    })
  }
  
  output$tabla_ajuste <- renderTable({
    filas <- data.frame(
      Métrica = "Pseudo R² (McFadden)",
      Valor = sprintf("%.4f", pseudo_r2)
    )
    if (pkg_proc) {
      filas <- bind_rows(filas, data.frame(Métrica = "AUC", Valor = sprintf("%.4f", as.numeric(auc(roc_obj)))))
    }
    if (pkg_hl) {
      filas <- bind_rows(filas, data.frame(
        Métrica = c("Hosmer-Lemeshow: estadístico", "Hosmer-Lemeshow: p-valor"),
        Valor = c(sprintf("%.4f", hl_test$statistic), format.pval(hl_test$p.value, digits = 4))
      ))
    }
    filas
  })
  
  # ---- VIF e índice de condición ----
  output$plot_vif <- renderPlot({
    validate(need(pkg_car, "Paquete 'car' no disponible: instálalo para calcular el VIF."))
    ggplot(vif_df, aes(x = reorder(variable, VIF_ajustado), y = VIF_ajustado)) +
      geom_col(fill = "#2C6E63", width = 0.6) +
      geom_hline(yintercept = 5, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = 10, linetype = "dashed", color = "#D9534F") +
      coord_flip() +
      labs(x = NULL, y = "VIF ajustado") +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank())
  })
  
  output$tabla_indice_condicion <- renderTable({
    data.frame(
      Componente = seq_along(autovalores),
      Autovalor = round(autovalores, 4),
      `Índice de condición` = round(indice_condicion, 2),
      check.names = FALSE
    )
  })
  
  # ---- Interacción y modelos anidados ----
  output$tabla_interaccion <- renderTable({
    as.data.frame(lr_test_interaccion) %>%
      tibble::rownames_to_column("Modelo") %>%
      mutate(across(where(is.numeric), ~ round(.x, 4)))
  })
  
  output$tabla_opiniones <- renderTable({
    as.data.frame(lr_test_opiniones) %>%
      tibble::rownames_to_column("Modelo") %>%
      mutate(across(where(is.numeric), ~ round(.x, 4)))
  })
  
  # ---- Ridge ----
  output$bloque_ridge <- renderUI({
    if (pkg_glmnet) tagList(plotOutput("plot_ridge", height = "420px"))
    else p(class = "text-muted", "Paquete 'glmnet' no disponible: instálalo para ver la comparación Ridge.")
  })
  
  if (pkg_glmnet) {
    output$plot_ridge <- renderPlot({
      df <- comparacion_coef %>%
        tidyr::pivot_longer(cols = c(coef_glm_original, coef_ridge),
                            names_to = "tipo", values_to = "coef")
      ggplot(df, aes(x = reorder(variable, coef), y = coef, fill = tipo)) +
        geom_col(position = "dodge", width = 0.6) +
        coord_flip() +
        scale_fill_manual(values = c(coef_glm_original = "#8FC1B5", coef_ridge = "#2C6E63"),
                          labels = c(coef_glm_original = "GLM original", coef_ridge = "Ridge")) +
        labs(x = NULL, y = "Coeficiente", fill = NULL) +
        theme_minimal(base_size = 11) +
        theme(panel.grid.minor = element_blank())
    })
  }
}

# --------------------------------------------------------------
# 5. Run
# --------------------------------------------------------------
shinyApp(ui = ui, server = server)