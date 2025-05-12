#' Evaluate model performance (NSE, KGE, C2MP, POD, FAR, CSI)
#'
#' @param results_cal A tibble containing the model results in calibration mode
#' @param results_test A tibble containing the model results in evaluation mode
#' @param events A tibble containing information about the flow events: start, end dates, and event IDs
#' @param Qthr Streamflow threshold to evaluate contingency table
#' @param Hprev Forecast horizon ("H3", "H6", "H12" or "H24")
#' @return A tibble containing the calculated performance metrics for the specified period, including:
#'   - `NSE_KGE bounded`: NSE and KGE metrics for overall performance and flood events
#'   - `C2MP`: C2MP metrics for overall performance and flood events
#'   - `cont_table`: Contingency table metrics (POD, FAR, CSI) for each model variable (Qprev, Qcorr, Qtan)
#' @export

evaluate_results <- function(
  results_cal, results_test, events, Qthr, Hprev) {
  
  #! Forecast horizon in hours
  Hprev_hours <- as.integer(substr(Hprev, 2, nchar(Hprev)))

  #! Forecast horizon in date-time format
  Hprev_datetime <- lubridate::hours(Hprev_hours) 

  #! Persistance (Qobs_tmH)
  results_cal <- results_cal |> 
    dplyr::mutate(
      Qobs_tmH = Qobs[match(Date - Hprev_datetime, Date)]
    ) |> 
    tidyr::drop_na()

  results_test <- results_test |> 
    dplyr::mutate(
      Qobs_tmH = Qobs[match(Date - Hprev_datetime, Date)]
    ) |> 
    tidyr::drop_na()
  
  Qobs <- results_cal |> dplyr::pull(Qobs)
  Dates <- results_cal |> dplyr::pull(Date)

  Qprev_cal <- results_cal |> dplyr::pull(Qprev)
  Qprev_test <- results_test |> dplyr::pull(Qprev)
  Qcorr_cal <- results_cal |> dplyr::pull(Qcorr)
  Qcorr_test <- results_test |> dplyr::pull(Qcorr)
  Qtan_cal <- results_cal |> dplyr::pull(Qtan)
  Qtan_test <- results_test |> dplyr::pull(Qtan)

  Qpers <- results_cal |> dplyr::pull(Qobs_tmH)
  
  #! Temporal masks
  mask_events <- purrr::map_lgl(
    Dates, 
    \(x) any(x >= events$start & x <= events$end)
  )

  mask_overall <- rep(TRUE, length(mask_events))

  #! Generate masks for each event
  mask_individual_events <- purrr::map(
    events$ID, 
    \(x) {
      start <- events |> 
        dplyr::filter(ID == x) |> 
        dplyr::pull(start)
      end <- events |> 
        dplyr::filter(ID == x) |> 
        dplyr::pull(end)
      Dates |> 
        dplyr::between(start, end)
    }) |>
    #! Combine the results into a matrix
    purrr::reduce(rbind)
  
  #! Bind all masks
  mask <- rbind(
    mask_overall, mask_events, mask_individual_events
  )

  #! Minimum time-steps required to validate the event performance (in hours)
  min_ts_event_threshold <- 12
  
  exclude_short_events <- function(perf_Qvar_mode) {
    NSE_Qvar_mode_values <- as.vector(perf_Qvar_mode[[1]])
    KGE_Qvar_mode_values <- as.vector(perf_Qvar_mode[[2]])
    completeness <- as.vector(perf_Qvar_mode[[3]])
    #! Identify events with time-steps below the threshold
    values_to_exclude <- completeness <= min_ts_event_threshold
    
    #! Assign NA to values corresponding to events that don't meet the minimum time-step threshold
    NSE_Qvar_mode_values[values_to_exclude] <- NA
    KGE_Qvar_mode_values[values_to_exclude] <- NA

    return(
      tibble::tibble(
        Mask = c("Overall performance", "Flood events", paste0("Event ", 1:length(events$ID))),
        NSE = NSE_Qvar_mode_values,
        KGE = KGE_Qvar_mode_values
      )
    ) 
  }

  #! Calculate NSE and KGE model efficiency coefficients 
  perf_Qprev_cal <- suppressWarnings(evalhyd::evald(
    q_obs = array(Qobs), 
    q_prd = array(Qprev_cal), 
    metrics = c("NSE", "KGE"),
    t_msk = array(mask, dim = c(1, dim(mask)[1], length(Qobs))),
    diagnostics = c("completeness")
  )) |>
    exclude_short_events() |> 
    dplyr::mutate(
      Variable = "Qprev",
      Mode = "Calage"
    )

  perf_Qprev_test <- suppressWarnings(evalhyd::evald(
    q_obs = array(Qobs), 
    q_prd = array(Qprev_test), 
    metrics = c("NSE", "KGE"),
    t_msk = array(mask, dim = c(1, dim(mask)[1], length(Qobs))),
    diagnostics = c("completeness")
  )) |>
    exclude_short_events() |> 
    dplyr::mutate(
      Variable = "Qprev",
      Mode = "Test"
    )

  perf_Qcorr_cal <- suppressWarnings(evalhyd::evald(
    q_obs = array(Qobs), 
    q_prd = array(Qcorr_cal), 
    metrics = c("NSE", "KGE"),
    t_msk = array(mask, dim = c(1, dim(mask)[1], length(Qobs))),
    diagnostics = c("completeness")
  )) |>
    exclude_short_events() |> 
    dplyr::mutate(
      Variable = "Qcorr",
      Mode = "Calage"
    )

  perf_Qcorr_test <- suppressWarnings(evalhyd::evald(
    q_obs = array(Qobs), 
    q_prd = array(Qcorr_test), 
    metrics = c("NSE", "KGE"),
    t_msk = array(mask, dim = c(1, dim(mask)[1], length(Qobs))),
    diagnostics = c("completeness")
  )) |>
    exclude_short_events() |> 
    dplyr::mutate(
      Variable = "Qcorr",
      Mode = "Test"
    )

  perf_Qtan_cal <- suppressWarnings(evalhyd::evald(
    q_obs = array(Qobs), 
    q_prd = array(Qtan_cal), 
    metrics = c("NSE", "KGE"),
    t_msk = array(mask, dim = c(1, dim(mask)[1], length(Qobs))),
    diagnostics = c("completeness")
  )) |>
    exclude_short_events() |> 
    dplyr::mutate(
      Variable = "Qtan",
      Mode = "Calage"
    )

  perf_Qtan_test <- suppressWarnings(evalhyd::evald(
    q_obs = array(Qobs), 
    q_prd = array(Qtan_test), 
    metrics = c("NSE", "KGE"),
    t_msk = array(mask, dim = c(1, dim(mask)[1], length(Qobs))),
    diagnostics = c("completeness")
  )) |>
    exclude_short_events() |> 
    dplyr::mutate(
      Variable = "Qtan",
      Mode = "Test"
    )
  
  #! Persistance criteria
  #! --------------------

  persistance <- function(Qforecast, mask) {

    data <- tibble::tibble(
      Qobs = Qobs,
      Qforecast = Qforecast,
      Qpers = Qpers,
      mask = mask
    ) |> 
    dplyr::filter(mask)

    Eff <- 1 - (sum((data$Qobs - data$Qforecast)^2) / sum((data$Qobs - data$Qpers)^2))
    C2MP <- Eff / (2 - Eff)
  }

  C2MP_Qprev_cal_events <- persistance(Qprev_cal, mask_events)
  C2MP_Qprev_cal_overall <- persistance(Qprev_cal, mask_overall)
  C2MP_Qprev_test_events <- persistance(Qprev_test, mask_events)
  C2MP_Qprev_test_overall <- persistance(Qprev_test, mask_overall)

  C2MP_Qcorr_cal_events <- persistance(Qcorr_cal, mask_events)
  C2MP_Qcorr_cal_overall <- persistance(Qcorr_cal, mask_overall)
  C2MP_Qcorr_test_events <- persistance(Qcorr_test, mask_events)
  C2MP_Qcorr_test_overall <- persistance(Qcorr_test, mask_overall)

  C2MP_Qtan_cal_events <- persistance(Qtan_cal, mask_events)
  C2MP_Qtan_cal_overall <- persistance(Qtan_cal, mask_overall)
  C2MP_Qtan_test_events <- persistance(Qtan_test, mask_events)
  C2MP_Qtan_test_overall <- persistance(Qtan_test, mask_overall)

  #! Contingency table
  #! -----------------
  
  #! Metric functions
  POD_function <- function(cont_table) {
    a <- cont_table[[1]][1]
    c <- cont_table[[1]][3]
    POD <- a / (a + c) * 100
  }
  
  FAR_function <- function(cont_table) {
    a <- cont_table[[1]][1]
    b <- cont_table[[1]][2]
    FAR <- b / (a + b) * 100
  }
  
  CSI_function <- function(cont_table) {
    a <- cont_table[[1]][1]
    b <- cont_table[[1]][2]
    c <- cont_table[[1]][3]
    CSI <- a / (a + b + c) * 100
  }
  
  #! Qprev
  ct_Qprev_cal <- evalhyd::evald(
    q_obs = array(Qobs), 
    q_prd = array(Qprev_cal), 
    q_thr = Qthr,
    events = "high",
    metrics = c("CONT_TBL")
  )
  POD_Qprev_cal <- POD_function(ct_Qprev_cal)
  FAR_Qprev_cal  <- FAR_function(ct_Qprev_cal)
  CSI_Qprev_cal  <- CSI_function(ct_Qprev_cal)
  
  ct_Qprev_test <- evalhyd::evald(
    q_obs = array(Qobs), 
    q_prd = array(Qprev_test), 
    q_thr = Qthr,
    events = "high",
    metrics = c("CONT_TBL")
  )
  POD_Qprev_test <- POD_function(ct_Qprev_test)
  FAR_Qprev_test  <- FAR_function(ct_Qprev_test)
  CSI_Qprev_test <- CSI_function(ct_Qprev_test)

  #! Qcorr
  ct_Qcorr_cal <- evalhyd::evald(
    q_obs = array(Qobs), 
    q_prd = array(Qcorr_cal), 
    q_thr = Qthr,
    events = "high",
    metrics = c("CONT_TBL")
  )
  POD_Qcorr_cal <- POD_function(ct_Qcorr_cal)
  FAR_Qcorr_cal  <- FAR_function(ct_Qcorr_cal)
  CSI_Qcorr_cal  <- CSI_function(ct_Qcorr_cal)
  
  ct_Qcorr_test <- evalhyd::evald(
    q_obs = array(Qobs), 
    q_prd = array(Qcorr_test), 
    q_thr = Qthr,
    events = "high",
    metrics = c("CONT_TBL")
  )
  POD_Qcorr_test <- POD_function(ct_Qcorr_test)
  FAR_Qcorr_test  <- FAR_function(ct_Qcorr_test)
  CSI_Qcorr_test <- CSI_function(ct_Qcorr_test)

  #! Qtan
  ct_Qtan_cal <- evalhyd::evald(
    q_obs = array(Qobs), 
    q_prd = array(Qtan_cal), 
    q_thr = Qthr,
    events = "high",
    metrics = c("CONT_TBL")
  )
  POD_Qtan_cal <- POD_function(ct_Qtan_cal)
  FAR_Qtan_cal  <- FAR_function(ct_Qtan_cal)
  CSI_Qtan_cal  <- CSI_function(ct_Qtan_cal)
  
  ct_Qtan_test <- evalhyd::evald(
    q_obs = array(Qobs), 
    q_prd = array(Qtan_test), 
    q_thr = Qthr,
    events = "high",
    metrics = c("CONT_TBL")
  )
  POD_Qtan_test <- POD_function(ct_Qtan_test)
  FAR_Qtan_test  <- FAR_function(ct_Qtan_test)
  CSI_Qtan_test <- CSI_function(ct_Qtan_test)

  #! Bounded version
  NSE_bounded <- function(NSE) {
    NSE_b <- NSE / (2 - NSE)
  }

  KGE_bounded <- function(KGE) {
    KGE_b <- KGE / (2 - KGE)
  }
  
  df_NSE_KGE_bounded <- dplyr::bind_rows(
    perf_Qprev_cal, 
    perf_Qprev_test,
    perf_Qcorr_cal,
    perf_Qcorr_test,
    perf_Qtan_cal,
    perf_Qtan_test
  ) |> 
  dplyr::mutate(
    NSE = NSE_bounded(NSE),
    KGE = KGE_bounded(KGE)
  )

  C2MP_Qprev_cal_events <- persistance(Qprev_cal, mask_events)
  C2MP_Qprev_cal_overall <- persistance(Qprev_cal, mask_overall)
  C2MP_Qprev_test_events <- persistance(Qprev_test, mask_events)
  C2MP_Qprev_test_overall <- persistance(Qprev_test, mask_overall)

  C2MP_Qcorr_cal_events <- persistance(Qcorr_cal, mask_events)
  C2MP_Qcorr_cal_overall <- persistance(Qcorr_cal, mask_overall)
  C2MP_Qcorr_test_events <- persistance(Qcorr_test, mask_events)
  C2MP_Qcorr_test_overall <- persistance(Qcorr_test, mask_overall)

  C2MP_Qtan_cal_events <- persistance(Qtan_cal, mask_events)
  C2MP_Qtan_cal_overall <- persistance(Qtan_cal, mask_overall)
  C2MP_Qtan_test_events <- persistance(Qtan_test, mask_events)
  C2MP_Qtan_test_overall <- persistance(Qtan_test, mask_overall)

  df_persistance <- tibble::tibble(
    Mask = rep(c("Overall performance", "Flood Events"), times = 6),
    Variable = rep(c("Qprev", "Qcorr", "Qtan"), each = 4),
    Mode = rep(c("Calage", "Test"), each = 2, times = 3),
    Value = c(
      C2MP_Qprev_cal_overall, C2MP_Qprev_cal_events,
      C2MP_Qprev_test_overall, C2MP_Qprev_test_events,

      C2MP_Qcorr_cal_overall, C2MP_Qcorr_cal_events,
      C2MP_Qcorr_test_overall, C2MP_Qtan_test_events,

      C2MP_Qtan_cal_overall, C2MP_Qtan_cal_events,
      C2MP_Qtan_test_overall, C2MP_Qtan_test_events
    )
  )

  df_cont_table <- tibble::tibble(
    Metric = rep(c("POD", "FAR", "CSI"), each = 6),
    Variable = rep(c("Qprev", "Qprev", "Qcorr", "Qcorr", "Qtan", "Qtan"), times = 3),
    Mode = rep(c("Calage", "Test"), times = 3 * 3), 
    Value = c(
      POD_Qprev_cal, POD_Qprev_test,
      POD_Qcorr_cal, POD_Qcorr_test,
      POD_Qtan_cal, POD_Qtan_test,
  
      FAR_Qprev_cal, FAR_Qprev_test,
      FAR_Qcorr_cal, FAR_Qcorr_test,
      FAR_Qtan_cal, FAR_Qtan_test,
  
      CSI_Qprev_cal, CSI_Qprev_test,
      CSI_Qcorr_cal, CSI_Qcorr_test,
      CSI_Qtan_cal, CSI_Qtan_test
    )
  )
  
  #! Save results
  df_eval_results <- tidyr::tibble(
    NSE_KGE_bounded = list(df_NSE_KGE_bounded),
    C2MP = list(df_persistance),
    cont_table = list(df_cont_table)
  )
  return(df_eval_results)
}
