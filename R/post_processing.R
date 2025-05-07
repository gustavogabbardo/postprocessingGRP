#' Post-processing Model Application (MLP, Random Forest, or XGBoost)
#'
#' @param input_data_cal A tibble containing the pre-processed GR5H_RI forecast data for the calibration period.
#' @param input_data_eval A tibble containing the pre-processed GR5H_RI forecast data for the evaluation period.
#' @param model_name A string indicating the type of post-processing model to apply ("mlp", "random_forest" or "xgboost").
#' @param Hprev Forecast horizon ("H3", "H6", "H12", "H24").
#' 
#' @return A list containing:
#'   - `cal`: The model's predictions for the calibration period.
#'   - `test`: The model's predictions for the evaluation period.
#'   - `model`: Fitted model and calibrated hyperparameters.
#'
#' @export

post_processing <- function(
  input_data_cal, input_data_eval, model_name, Hprev
) {
  
  #! Forecast horizon in hours
  Hprev_hours <- as.integer(substr(Hprev, 2, nchar(Hprev)))

  #! Split data (training, validation and testing)
  #! Initial split
  data_split <- rsample::initial_time_split(input_data_cal, prop = 0.5)
  
  #! Datasets
  data_train <- data_split |> rsample::training()
  data_test <- data_split |> rsample::testing()
  
  #! Resampling
  resample <- rsample::manual_rset(
    splits = list(rsample::make_splits(data_train, assessment = data_test)),
    ids = "Validation"
  )

  #! Model recipe
  recipe <- recipes::recipe(
    #! Predict Err_tpH_t based on all selected variables from "data"
    Err_tpH_t ~ .,         
    data = data_train |> 
      #! Select variables (output = error at t+H ; inputs = errors, observed flows)
      dplyr::select(Qobs_t, Qobs_tm1, dQobs_t_tm1, tidyselect::starts_with("Err"))
  )

  #! Normalize all numerical predictor variables (only for MLP model)
  if (model_name == "mlp") {
    recipe <- recipe |>
      recipes::step_normalize(recipes::all_predictors())  
  }

  #! Machine learning model
  model_mlp <- parsnip::mlp(
    hidden_units = c(64),
    learn_rate = c(0.1),
    epochs = c(10)
  ) |> 
    parsnip::set_engine("brulee") |> 
    parsnip::set_mode("regression")

  model_random_forest <- parsnip::rand_forest(
    trees = 100,
    min_n = 5
  ) |> 
    parsnip::set_engine("ranger") |> 
    parsnip::set_mode("regression")

  model_xgboost <- parsnip::boost_tree(
    trees = 100,
    tree_depth = 6,
    learn_rate = 0.3,
    loss_reduction = 0.2,
    min_n = 5
  ) |> 
    parsnip::set_engine("xgboost") |> 
    parsnip::set_mode("regression")

  model <- switch(
    model_name,
    "mlp" = model_mlp,
    "random_forest" = model_random_forest,
    "xgboost" = model_xgboost,
    stop("Post-processing model must be 'mlp', 'random_forest' or 'xgboost'.")
  )

  #! Workflow
  workflow <- workflows::workflow() |> 
    workflows::add_model(model) |> 
    workflows::add_recipe(recipe)

  #! Fit model
  #! Model will be fit using GR5H_RI results in calibration mode (ensemble of the training and validation data)
  model_fit <- workflow |> 
    parsnip::fit(input_data_cal)

  #! Model predictions function
  make_predictions <- function(input_data) {
    
    #! Select variables (output = error at t+H ; inputs = errors at t, t-1, ..., t - time_steps)
    input_data_selected <- input_data |> 
      dplyr::select(Qobs_t, Qobs_tm1, dQobs_t_tm1, tidyselect::starts_with("Err")) 
    
    #! Predict errors at t+H using the fitted model
    predictions <- stats::predict(model_fit, new_data = input_data_selected)
    
    #! Create a tibble with corresponding dates and predicted errors at t+H
    error_predictions <- tidyr::tibble(
      Date = input_data$Date,
      Err_prev_tpH_t = predictions$.pred
    )
    
    #! Use predict errors to correct streamflow predictions
    results <- input_data |> 
      dplyr::inner_join(error_predictions, dplyr::join_by(Date)) |> 
      dplyr::mutate(
        #! Compute the preliminary corrected flow
        Qcorr_raw = Qprev_tpH_t - Err_prev_tpH_t,
        #! Apply constraints to prevent negative values and excessive corrections
        Qcorr_tpH_t = dplyr::if_else(
          Qprev_tpH_t == 0, 
          0,
          dplyr::if_else(
            Qcorr_raw / Qprev_tpH_t >= 1.5,
            Qprev_tpH_t * 1.5,
            dplyr::if_else(
              Qcorr_raw / Qprev_tpH_t <= 0.5,
              Qprev_tpH_t * 0.5,
              Qcorr_raw
            )
          ) 
        ),
        #! Add Hprev hours to shift predictions accordingly
        Date = Date + lubridate::hours(Hprev_hours)
      ) |> 
      dplyr::select(Date, Err_tpH_t, Qobs_tpH, Qprev_tpH_t, Qcorr_tpH_t, Qtan_tpH_t) |> 
      dplyr::rename(
        Qobs = Qobs_tpH,
        Qprev = Qprev_tpH_t,
        Qtan = Qtan_tpH_t,
        Qcorr = Qcorr_tpH_t
      )
  }
  
  #! Predict errors in the calibration period
  res_cal <- make_predictions(input_data_cal)

  #! Predict errors in the test period
  res_test <-  make_predictions(input_data_eval)

  return(list(
    cal = res_cal,
    test = res_test,
    model = model_fit
  ))
}