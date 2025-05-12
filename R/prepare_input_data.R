#' Prepare Input Data for Post-Processing Model
#' 
#' @param path_prev Path to the GR5H_RI forecast data (RDS file).
#' @param time_steps Time window, in hours, used to compute past streamflow prediction errors for the post-processing model.
#' @param Hprev Forecast horizon ("H3", "H6", "H12" or "H24").
#' @param Qobs_ts Observed streamflow time series tibble with columns "Date" (dttm) and "Qobs" (dbl).
#' @param cal_method Calibration method ("WsRf, "Ref" or "OL")
#'
#' @return A tibble containing pre-processed data for use in the post-processing model, including time series 
#' of observed and forecasted streamflow.
#' @export

prepare_input_data <- function(
  path_prev, time_steps, Hprev, Qobs_ts, cal_method
) {
  
  #! Number of lags calculated from the error at time t predicted at time t-1 (Err_t_t-1)
  n_lags <- time_steps - 1

  #! Forecast horizon in hours
  Hprev_hours <- as.integer(substr(Hprev, 2, nchar(Hprev)))

  #! Forecast horizon in date-time format
  Hprev_datetime <- lubridate::hours(Hprev_hours) 

  #! Read forecasts into a list and then convert to a tibble
  prev <- readr::read_rds(path_prev) 
  prev <- tidyr::as_tibble(prev$out_Run[1:26]) |> 
    dplyr::filter(H1 >= 0) |> #! Remove GR5H_RI missing predictions
    dplyr::rename(Date = DatesR)
  
  #! Join observation data (Qobs and Pmm) with GR5H_RI model predictions
  prev <- prev |> 
    dplyr::select(Date, H1, OL, tidyr::all_of(Hprev)) |> 
    dplyr::left_join(Qobs_ts, dplyr::join_by(Date))
  
  #! Exclude data from October, November, and December 2021 (issues with Comephore data)
  prev <- prev |> 
    dplyr::mutate(
      year = lubridate::year(Date),
      month = lubridate::month(Date)
    ) |> 
    dplyr::filter(!(year == 2021 & month %in% c(10, 11, 12))) |> 
    dplyr::select(-year, -month) |> 
    dplyr::rename(Qobs_t = Qobs)
  
  #! Calculate streamflow threshold (to evaluate contingency table) 
  Qthr <- stats::quantile(
    prev |> 
      dplyr::pull(Qobs_t), 
      probs = prob_thr,
      na.rm = TRUE
  )
  
  #! Prepare model input data 
  if (cal_method != "OL") {
    prev <- prev |>
      dplyr::mutate(
        #! Observed flow at time t+H: Qobs(t+H)
        Qobs_tpH = Qobs_t[match(Date + Hprev_datetime, Date)],  
        #! Observed flow at time t-1: Qobs(t-1)
        Qobs_tm1 = Qobs_t[match(Date - lubridate::hours(1), Date)],  
        #! Observed gradient flow between t and t-1
        dQobs_t_tm1 = Qobs_t - Qobs_tm1,
        #! Predicted flow at time t forecasted at time t-H: Qprev(t, t-H)
        Qprev_t_tmH = prev[[Hprev]][match(Date - Hprev_datetime, Date)],
        #! Predicted flow at time t+H forecasted at time t: Qprev(t+H, t)
        Qprev_tpH_t = prev[[Hprev]],   
        #! Error at time t+H predicted at time t: Err(t+H, t)
        Err_tpH_t = Qprev_tpH_t - Qobs_tpH,     
        #! Predicted flow at time t forecasted at time t-1: Qprev(t, t-1)
        Qprev_t_tm1 = H1[match(Date - lubridate::hours(1), Date)],
        #! Error at time t predicted at time t-1: Err(t, t-1)
        Err_t_tm1 = Qprev_t_tm1 - Qobs_t   
      ) 
  } else {
    prev <- prev |> 
      dplyr::mutate(
        #! Observed flow at time t+H: Qobs(t+H)
        Qobs_tpH = Qobs_t[match(Date + Hprev_datetime, Date)],  
        #! Observed flow at time t-1: Qobs(t-1)
        Qobs_tm1 = Qobs_t[match(Date - lubridate::hours(1), Date)],  
        #! Observed gradient flow between t and t-1
        dQobs_t_tm1 = Qobs_t - Qobs_tm1,
        #! Predicted flow at time t forecasted at time t-H: Qprev(t, t-H)
        Qprev_t_tmH = OL,
        #! Predicted flow at time t+H forecasted at time t: Qprev(t+H, t)
        Qprev_tpH_t = OL[match(Date + Hprev_datetime, Date)],   
        #! Error at time t+H predicted at time t: Err(t+H, t)
        Err_tpH_t = Qprev_tpH_t - Qobs_tpH,   
        #! Predicted flow at time t forecasted at time t-1: Qprev(t, t-1)
        Qprev_t_tm1 = OL,
        #! Error at time t predicted at time t-1: Err(t, t-1)
        Err_t_tm1 = Qprev_t_tm1 - Qobs_t  
      )
  }     

  #! Apply Tangara correction method
  #! The correction factor should be within the range of 0.2 to 5 to prevent excessive corrections
  prev <- prev |> 
    dplyr::mutate(
      #! Tangara correction for predicted flow at time t+H forecasted at time t: Qtan(t+H, t)
      Qtan_tpH_t = dplyr::if_else(
        Qprev_t_tm1 <= 0, 
        NA, 
        dplyr::if_else(
          (Qobs_t / Qprev_t_tm1) ^ 0.45 >= 5, 
          Qprev_tpH_t * 5, 
          dplyr::if_else(
            (Qobs_t / Qprev_t_tm1) ^ 0.45 <= 0.2, 
            Qprev_tpH_t * 0.2, 
            Qprev_tpH_t * (Qobs_t / Qprev_t_tm1) ^ 0.45
          )
        ),
      ), 
      #! Tangara correction for predicted flow at time t forecasted at time t-H: Qtan(t, t-H)
      Qtan_t_tmH = Qtan_tpH_t[match(Date - Hprev_datetime, Date)]
    ) 

  #! Create input variables errors
  for (i in 1:n_lags) {
    prev <- prev |>
      dplyr::mutate(
        !!paste0("Err_tm", i, "_tm", i+1) := Err_t_tm1[match(Date - lubridate::hours(i), Date)]
      )
  }
  
  #! Drop NA's after calculate input errors
  prev <- prev |>
    tidyr::drop_na()

  return(prev)
}
