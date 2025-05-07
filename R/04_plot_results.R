#! Load packages
library(tidyverse)
library(patchwork)

#! Main directories and paths
#! --------------------------

#! Input data
dir_input <- file.path("C:/Gustavo/Etudes/20241104_Post-traitement_Previsions_GR5H_RI/02_DATA/input")

#! Output data
dir_output <- file.path("C:/Gustavo/Etudes/20241104_Post-traitement_Previsions_GR5H_RI/02_DATA/output")
cal_method <- "WsRf"
path_results <- file.path(dir_output, "02_post_processing", cal_method)

#! Plots
dir_plot <- file.path("C:/Gustavo/Etudes/20241104_Post-traitement_Previsions_GR5H_RI/04_COMM/02_RAPPORT/Figures", cal_method)
if (!dir.exists(dir_plot )) dir.create(dir_plot, recursive = TRUE)

#! Read catchments data
#! --------------------

#! Catchments list
catchments_list_path <- file.path(dir_input, "liste_bv.rds")
list_bv <- readr::read_rds(catchments_list_path)

#! Catchments properties data (code, surface, altitude)
synth_fin <- get(load(file.path(dir_input, "catchments_info.RData")))

data_list_bv <- synth_fin |> 
  dplyr::filter(CODE_H3 %in% list_bv) |> 
  dplyr::select(CODE_H3, LABEL_H3, SURF_IR, XL93_IR, YL93_IR, ALTI_IR)

#! Response time
response_time <- readr::read_rds(file.path(dir_input, "response_times.rds")) |> 
  dplyr::as_tibble()

#! Read results
#!-------------

models <- c("xgboost", "random_forest", "mlp")

results_bv <- tibble::tibble()
for (model_name in models) {
  files <- list.files(file.path(path_results, model_name), full.names = TRUE)
  for (file in files) {
    bv <- stringr::str_sub(file, -14, -5)
    results_horizons <- readr::read_rds(file) |> 
      dplyr::mutate(Code = bv, model_name = model_name) |> 
      dplyr::left_join(response_time, dplyr::join_by(Code)) |> 
      dplyr::select(Code, Hprev, model_name, TimeC, eval) 
  
    results_bv <- results_bv |> 
      dplyr::bind_rows(results_horizons)
  
  }
}

#! NSE and KGE
results_NSE_KGE <- results_bv |> 
    tidyr::unnest(eval) |> 
    tidyr::unnest(NSE_KGE_bounded) |> 
    dplyr::select(-cont_table) |> 
    dplyr::filter(Mask %in% c("Overall Performance", "Flood events")) |> 
    dplyr::mutate(
      Variable = dplyr::case_when(
        Variable == "Qprev" ~ "GR5H_RI",
        Variable == "Qcorr" ~ model_name,
        Variable == "Qtan" ~ "Tangara"
      ),
      Variable = factor(
        Variable, 
        levels = c("GR5H_RI", "Tangara", "xgboost", "random_forest", "mlp")
      ),
      Variable = forcats::fct_recode(
        Variable,
        "XGBoost" = "xgboost", 
        "Random Forest" = "random_forest", 
        "MLP" = "mlp"
      ),
      Hprev = factor(
        Hprev, 
        levels = c("H3", "H6", "H12", "H24")
      ),
      Hprev = forcats::fct_recode(
        Hprev,
        "3 h" = "H3", 
        "6 h" = "H6", 
        "12 h" = "H12", 
        "24 h" = "H24", 
      ),
      Mode = dplyr::case_when(
        Mode == "Calage" ~ "Calage",
        Mode == "Test" ~ "Évaluation"
      ),
      Mask = dplyr::case_when(
        Mask == "Overall Performance" ~ "Performance globale",
        Mask == "Flood events" ~ "Évènements de crue"
      ),
      TimeC = dplyr::if_else(
        dplyr::between(TimeC, 0, 9), 
        "Tr [1, 9] h",
        dplyr::if_else(
          dplyr::between(TimeC, 10, 18),
            "Tr [10, 18] h",
            "Tr [19, 128] h",
        )
      )
    ) |> 
    dplyr::select(-model_name) |> 
    dplyr::distinct(Code, Hprev, TimeC, Mask, Variable, Mode, .keep_all = TRUE) |> 
    dplyr::arrange(Code, Hprev) 

#! Contingency table
results_cont_table <- results_bv |> 
  tidyr::unnest(eval) |> 
  tidyr::unnest(cont_table) |> 
  dplyr::select(-NSE_KGE_bounded) |> 
  dplyr::mutate(
    Variable = dplyr::case_when(
      Variable == "Qprev" ~ "GR5H_RI",
      Variable == "Qcorr" ~ model_name,
      Variable == "Qtan" ~ "Tangara"
    ),
    Variable = factor(
      Variable, 
      levels = c("GR5H_RI", "Tangara", "xgboost", "random_forest", "mlp")
    ),
    Variable = forcats::fct_recode(
      Variable,
      "XGBoost" = "xgboost", 
      "Random Forest" = "random_forest", 
      "MLP" = "mlp"
    ),
    Hprev = factor(
      Hprev, 
      levels = c("H3", "H6", "H12", "H24")
    ),
    Hprev = forcats::fct_recode(
      Hprev,
      "3 h" = "H3", 
      "6 h" = "H6", 
      "12 h" = "H12", 
      "24 h" = "H24", 
    ),
    Mode = dplyr::case_when(
      Mode == "Calage" ~ "Calage",
      Mode == "Test" ~ "Évaluation"
    )
  ) |> 
  dplyr::select(-model_name) |> 
  dplyr::distinct(Code, Hprev, TimeC, Metric, Variable, Mode, .keep_all = TRUE) |> 
  dplyr::arrange(Code, Hprev, Metric) 

#! Compare performances for NSE or KGE (Qprev x Qcorr x Qtan)
#! ----------------------------------------------------------
scatter_nse_kge <- function(mode, metric, model_x, model_y) {

  data <- results_NSE_KGE |>
    dplyr::filter(Mode == mode) |> 
    dplyr::select(Code, Hprev, Mask, {{ metric }}, Variable) |> 
    dplyr::filter(Variable  %in% c(model_x, model_y)) |> 
    tidyr::pivot_wider(
        names_from = Variable,
        values_from = {{ metric }}
    )

  count <- data |> 
    dplyr::mutate(
      compare_model = dplyr::case_when(
        {{ model_y }} > {{ model_x }}  ~ -1, 
        {{ model_y }} < {{ model_x }}  ~  1, 
        {{ model_y }} == {{ model_x }} ~  0  
      )
    ) |> 
    dplyr::group_by(Hprev, Mask, compare_model) |> 
    dplyr::summarise(n = n(), .groups = "drop")

  data |> 
    ggplot2::ggplot(ggplot2::aes(x = !!sym(model_x), y = !!sym(model_y))) +
    ggplot2::geom_point(alpha = 0.15) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red", linetype = 2) +
    ggplot2::facet_grid(Mask~Hprev) +
    ggplot2::coord_equal(xlim = c(-0.5, 1), ylim = c(-0.5, 1)) +
    ggplot2::labs(
      x = model_x,
      y = model_y
    ) +
    ggplot2::geom_text(
      data = count, 
      ggplot2::aes(
        x = dplyr::case_when(
          compare_model == -1 ~ -0.3, 
          compare_model == 1  ~ 0.8,
          compare_model == 0  ~ 0 
        ),
        y = case_when(
          compare_model == -1 ~ 0.8, 
          compare_model == 1  ~ -0.3,
          compare_model == 0  ~ 0 
        ),
        label = paste0("n = ", n),
        color = factor(compare_model)
      ),
      size = 3
    ) +
    ggplot2::scale_color_manual(
      values = c("-1" = "black", "1" = "black", "0" = "red")
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      panel.spacing = ggplot2::unit(0.3, "cm")
    )
}

scatter_nse_kge("Évaluation", "NSE", "Random Forest", "GR5H_RI")
ggplot2::ggsave(file.path(dir_plot, "/NSE_RF_GR5H_RI_test.png"), width = 9, height = 7, dpi = 600)

scatter_nse_kge("Évaluation", "KGE", "Random Forest", "GR5H_RI")
ggplot2::ggsave(file.path(dir_plot, "/KGE_RF_GR5H_RI_test.png"), width = 9, height = 7, dpi = 600)


scatter_nse_kge("Évaluation", "NSE", "Random Forest", "Tangara")
ggplot2::ggsave(file.path(dir_plot, "/NSE_RF_Tan_test.png"), width = 9, height = 7, dpi = 600)

scatter_nse_kge("Évaluation", "KGE", "Random Forest", "Tangara")
ggplot2::ggsave(file.path(dir_plot, "/KGE_RF_Tan_test.png"), width = 9, height = 7, dpi = 600)


#! Compare performances for contingency table metrics (Qprev x Qcorr x Qtan)
#! -------------------------------------------------------------------------

scatter_cont_table <- function(mode, model_x, model_y) {

  data <- results_cont_table |>
    dplyr::filter(Variable  %in% c(model_x, model_y)) |> 
    tidyr::pivot_wider(
        names_from = Variable,
        values_from = Value
    )

  count <- data |> 
    dplyr::mutate(
      compare_model = dplyr::case_when(
        {{ model_y }} > {{ model_x }}  ~ -1, 
        {{ model_y }} < {{ model_x }}  ~  1, 
        {{ model_y }} == {{ model_x }} ~  0  
      )
    ) |> 
    dplyr::group_by(Hprev, Metric, compare_model) |> 
    dplyr::summarise(n = n(), .groups = "drop")

  data|> 
    ggplot2::ggplot(ggplot2::aes(x = !!sym(model_x) , y = !!sym(model_y))) +
    ggplot2::geom_point(size = 1, alpha = 0.15) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red", linetype = 2) +
    ggplot2::facet_grid(Metric~Hprev) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(scale = 1) 
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_percent(scale = 1) 
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = model_x,
      y = model_y
    ) +
    ggplot2::geom_text(
      data = count, 
      ggplot2::aes(
        x = dplyr::case_when(
          compare_model == -1 ~ 10, 
          compare_model == 1  ~ 90,  
          compare_model == 0  ~ 10 
        ),
        y = case_when(
          compare_model == -1 ~ 90, 
          compare_model == 1  ~ 10,
          compare_model == 0  ~ 25 
        ),
        label = paste0("n = ", n),
        color = factor(compare_model)
      ),
      size = 3
    ) +
    ggplot2::scale_color_manual(
      values = c("-1" = "black", "1" = "black", "0" = "red")
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      panel.spacing = ggplot2::unit(0.3, "cm"),
      axis.text.x = element_text(angle = 90, hjust = 1)
    )
}

scatter_cont_table("Évaluation", "Random Forest", "GR5H_RI")
ggplot2::ggsave(file.path(dir_plot, "/CT_RF_GR5H_RI_test.png"), width = 9, height = 7, dpi = 600)

scatter_cont_table("Évaluation", "Random Forest", "Tangara")
ggplot2::ggsave(file.path(dir_plot, "/CT_RF_Tan_test.png"), width = 9, height = 7, dpi = 600)

#! Boxplot NSE/KGE
#! ---------------

boxplot_plot <- function(metric) {
  results_NSE_KGE |> 
    dplyr::select(Code, Hprev, Mask, {{ metric }}, Variable, Mode) |> 
    dplyr::group_by(Variable, Hprev, Mask, Mode) |> 
    dplyr::summarise(
      Q05 = quantile(!!sym(metric), 0.05, na.rm = TRUE),
      Q25 = quantile(!!sym(metric), 0.25, na.rm = TRUE),
      Q50 = quantile(!!sym(metric), 0.50, na.rm = TRUE),
      Q75 = quantile(!!sym(metric), 0.75, na.rm = TRUE),
      Q95 = quantile(!!sym(metric), 0.95, na.rm = TRUE)
    ) |> 
    ggplot2::ggplot(ggplot2::aes(x = Variable, fill = Mode)) +
    ggplot2::geom_boxplot(
      ggplot2::aes(
        ymin = Q05, lower = Q25, middle = Q50, upper = Q75, ymax = Q95
      ), 
      alpha = 0.5, 
      stat = 'identity'
    ) +
    ggplot2::facet_grid(Hprev ~ Mask, scales = "free_y") +
    ggplot2::labs(
      x = NULL,
      y = paste0(metric, " [-]"),
      fill = NULL,
      color = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Calage" = "#1B9E77",
        "Évaluation"= "#7570B3"
      )
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

boxplot_plot("NSE")
ggplot2::ggsave(file.path(dir_plot, "boxplot_NSE.png"), width = 9, height = 7, dpi = 300)

boxplot_plot("KGE")
ggplot2::ggsave(file.path(dir_plot, "boxplot_KGE.png"), width = 9, height = 7, dpi = 300)

#! Boxplot NSE/KGE  x response time
#! ---------------------------------

boxplot_tr <- function(metric, temps_reponse, show_strip) {
  p <- results_NSE_KGE |> 
    dplyr::filter(Mask == "Évènements de crue") |> 
    dplyr::filter(TimeC == temps_reponse) |> 
    dplyr::group_by(Variable, Hprev, TimeC) |> 
    dplyr::summarise(
      Q05 = quantile(!!sym(metric), 0.05, na.rm = TRUE),
      Q25 = quantile(!!sym(metric), 0.25, na.rm = TRUE),
      Q50 = quantile(!!sym(metric), 0.50, na.rm = TRUE),
      Q75 = quantile(!!sym(metric), 0.75, na.rm = TRUE),
      Q95 = quantile(!!sym(metric), 0.95, na.rm = TRUE)
    ) |> 
    ggplot2::ggplot(ggplot2::aes(y = factor(Variable, levels = rev(levels(Variable))), fill = Variable)) +
    ggplot2::geom_boxplot(
      ggplot2::aes(
        xmin = Q05, xlower = Q25, xmiddle = Q50, xupper = Q75, xmax = Q95,
        group = interaction(Variable, Hprev)
      ), 
      alpha = 0.5, 
      stat = 'identity'
    ) +
    ggplot2::facet_grid(TimeC ~ Hprev, scales = "free") +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      fill = NULL,
      color = NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      strip.text.x = if (!show_strip) ggplot2::element_blank() else ggplot2::element_text(),
      plot.margin = ggplot2::unit(c(0.3, 0, 0, 0), "cm"),
      panel.spacing.x = ggplot2::unit(0.5, "cm")
    )
}

p1 <- boxplot_tr("NSE", "Tr [1, 9] h", show_strip = TRUE)
p2 <- boxplot_tr("NSE", "Tr [10, 18] h", show_strip = FALSE)
p3 <- boxplot_tr("NSE", "Tr [19, 128] h", show_strip = FALSE)

(p1 / p2 / p3)
ggplot2::ggsave(file.path(dir_plot, "boxplot_NSE_Tr.png"), width = 9, height = 7, dpi = 300)

p1 <- boxplot_tr("KGE", "Tr [1, 9] h", show_strip = TRUE)
p2 <- boxplot_tr("KGE", "Tr [10, 18] h", show_strip = FALSE)
p3 <- boxplot_tr("KGE", "Tr [19, 128] h", show_strip = FALSE)

(p1 / p2 / p3)
ggplot2::ggsave(file.path(dir_plot, "boxplot_KGE_Tr.png"), width = 9, height = 7, dpi = 300)


#! Boxplot CT results
#! ------------------
results_cont_table |> 
  dplyr::filter(Mode == "Évaluation") |> 
  dplyr::group_by(Variable, Hprev, Metric) |> 
  dplyr::summarise(
    Q05 = quantile(Value, 0.05, na.rm = TRUE),
    Q25 = quantile(Value, 0.25, na.rm = TRUE),
    Q50 = quantile(Value, 0.50, na.rm = TRUE),
    Q75 = quantile(Value, 0.75, na.rm = TRUE),
    Q95 = quantile(Value, 0.95, na.rm = TRUE)
  ) |> 
  ggplot2::ggplot(ggplot2::aes(x = Variable, fill = Variable)) +
  ggplot2::geom_boxplot(
    ggplot2::aes(
      ymin = Q05, lower = Q25, middle = Q50, upper = Q75, ymax = Q95
    ), 
    alpha = 0.5, 
    stat = 'identity'
  ) +
  ggplot2::facet_grid(Metric ~ Hprev, scales = "free_y") +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    fill = NULL,
    color = NULL
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_percent(scale = 1) 
  ) + 
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot2::ggsave(file.path(dir_plot, "boxplot_ct.png"), width = 9, height = 7, dpi = 300)

#! Hydrographs
#! -----------

bv <- "M601401010"
model_name <- "random_forest"
event_id <- 3

plot_hydrograph_event <- function(bv, event_id, model_name) {

  #! Catchment flow events
  path_bv_events <- file.path(dir_input, "events_selection", paste0(bv, ".rds"))

  events <- readr::read_rds(path_bv_events)

  events_flows_info <- purrr::map_df(
    events$Date, 
    \(x) dplyr::tibble(start = min(x), end = max(x))
    ) |>
    dplyr::mutate(
      ID = dplyr::row_number()
    )

  #! COMEPHORE data
  path_Prec <- file.path(dir_input, paste0("hourly_hydro_database/P_TS_hourly/", "chroPrecip_", bv, ".Rdata"))

  chroPrecip <- get(load(path_Prec))
  Pmm <- tidyr::as_tibble(chroPrecip) |> 
    dplyr::select(dateR, pmm) |> 
    dplyr::rename(
      Date = dateR,
      Pmm = pmm
    )

  #! Forecasts results 
  path_forecasts <- file.path(path_results, model_name, paste0(bv, ".rds"))

  prev <- readr::read_rds(path_forecasts) |> 
    tidyr::unnest(test) |> 
    dplyr::select(Date, Hprev, Qobs, Qprev, Qcorr, Qtan) |> 
    dplyr::left_join(Pmm, join_by(Date)) |> 
    dplyr::mutate(
      Hprev = factor(Hprev, levels = c("H3", "H6", "H12", "H24")),
      Hprev = fct_recode(Hprev, "3 h" = "H3", "6 h" = "H6", "12 h" = "H12", "24 h" = "H24")
    )

  Q95 <- quantile(prev$Qobs, probs = 0.95, na.rm = TRUE)
  lim_y <- max(prev$Qobs, prev$Qprev, prev$Qcorr, na.rm = TRUE) * 1.1
  label_bv <- data_list_bv |> filter(CODE_H3 == bv) |> dplyr::pull(LABEL_H3)
  surf_bv <- data_list_bv |> filter(CODE_H3 == bv) |> dplyr::pull(SURF_IR)
  surf_bv <- round(surf_bv, 0)

  results_H3_H6 <- prev |> 
    dplyr::filter(Hprev %in% c("3 h", "6 h"))

  results_H12_H24 <- prev |> 
    dplyr::filter(Hprev %in% c("12 h", "24 h"))

  plot_hydrograph <- function(data, legend_position, Q95) {

    start <- events_flows_info |> filter(ID == event_id) |> dplyr::pull(start)
    end <- events_flows_info |> filter(ID == event_id) |> dplyr::pull(end)
    
    data <- data |> 
      dplyr::filter(Date >= start & Date <= end)

    hydrograph <- data |> 
      dplyr::select(-Qtan) |> 
      tidyr::pivot_longer(
        cols = c(Qobs, Qprev, Qcorr), 
        names_to = "Variable", 
        values_to = "Value"
      ) |> 
      dplyr::mutate(
        Variable = forcats::fct_recode(
          Variable,
          "Observé" = "Qobs", 
          "GR5H_RI" = "Qprev", 
          "Random Forest" = "Qcorr",
        ),
        Variable = factor(Variable, levels = c("Observé", "GR5H_RI", "Random Forest")),
      ) |> 
      ggplot2::ggplot() +
      ggplot2::geom_line(aes(x = Date, y = Value, color = Variable, linewidth = Variable)) +
      ggplot2::geom_hline(aes(yintercept = Q95, color = "Q95"), linetype = 2) + 
      ggplot2::annotate("segment", x = start, xend = end, y = 0, yend = 0, color = "black") +
      ggplot2::annotate("segment", x = start, xend = start, y = -Inf, yend = Inf, color = "black") +
      ggplot2::annotate("segment", x = end, xend = end, y = -Inf, yend = Inf, color = "black") +
      ggplot2::facet_wrap(~ Hprev, nrow = 1) +
      ggplot2::scale_x_datetime(date_labels = "%d/%m/%y") + 
      ggplot2::scale_y_continuous(limits = c(0, lim_y)) +
      ggplot2::coord_cartesian(xlim = c(start, end), expand = FALSE) +
      ggplot2::scale_color_manual(values = c(
        "Observé" = scales::alpha("gray30", 0.4), 
        "GR5H_RI" = "#d410aa",  
        "Random Forest" = "#0677b8",
        "Q95" = "black"
      )) +
      ggplot2::scale_linewidth_manual(values = c(
        "Observé" = 1.5,  
        "GR5H_RI" = 0.7, 
        "Random Forest" = 0.7
      ), guide = "none") +
      ggplot2::labs(
        x = NULL,
        y = "Q (mm/h)",
        color = NULL
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = legend_position, 
        legend.text = element_text(size = 11),
        panel.border = element_blank(),
        plot.margin = unit(c(0, 0, 0.5, 0),"cm"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x = unit(3, "lines"),
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9)
      )

    hyetograph <- data |> 
      filter(Date >= start & Date <= end) |> 
      ggplot2::ggplot() +
      ggplot2::geom_col(aes(x = Date, y = Pmm), fill = "dodgerblue4", color = "dodgerblue4") +
      ggplot2::annotate("segment", x = start, xend = start, y = -Inf, yend = Inf, color = "black") +
      ggplot2::annotate("segment", x = end, xend = end, y = -Inf, yend = Inf, color = "black") +
      ggplot2::labs(
        x = NULL,
        y = "P (mm/h)"
      ) +
      ggplot2::scale_y_reverse() +
      ggplot2::facet_wrap(~ Hprev, nrow = 1) +
      ggplot2::scale_x_datetime(date_labels = "%d/%m/%y") +
      ggplot2::coord_cartesian(xlim = c(start, end), expand = FALSE) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.margin = unit(c(0, 0, 0, 0),"cm"),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x = unit(3, "lines"),
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 11)
      )

    final_plot <- hyetograph + hydrograph +
      patchwork::plot_layout(nrow = 2, heights = c(0.4, 1))
  }

  (plot_hydrograph(results_H3_H6, "none", Q95) / plot_hydrograph(results_H12_H24, "bottom", Q95)) + 
    patchwork::plot_annotation(
      title = paste0(label_bv, " (", round(surf_bv, 2), " km²)"),
      theme = theme(plot.title = element_text(hjust = 0.5))
    )
}

plot_hydrograph_event(
  bv = bv,
  model_name =  model_name,
  event_id = event_id
)

ggplot2::ggsave(file.path(dir_plot, paste0("events/", bv, ".png")), width = 9, height = 7, dpi = 300)
