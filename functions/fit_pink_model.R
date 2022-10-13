fit_pink_model <- function(test_year,
                           model_type = "rand_forest",
                           run_wide = TRUE,
                           data,
                           scalar = 1,
                           trees = 1000,
                           log_returns = FALSE,
                           weight = FALSE,
                           use_years = TRUE,
                           delta_returns = FALSE,
                           n_mtry = 3,
                           initial_prop = 0.8,
                           assess = 1,
                           factor_years = FALSE,
                           produce = "summary",
                           forecast = FALSE,
                           p_done = "who knows how close to") {
  test_year <- lubridate::ymd(paste(test_year, "07", "03"))
  
  data$year <-
    lubridate::ymd(paste(data$year, "07", "03"))
  
  data$stock <- paste(data$subregion, data$cycle, sep = "_")
  
  data$subregion <- NULL
  
  # data$cycle <- NULL
  
  options(dplyr.summarise.inform = FALSE)
  
  if (forecast) {
    test_year <- test_year + 1
  }
  
  if (log_returns) {
    scalar <- 1
  }
  
  data$returns = data$returns / scalar
  
  data <- data %>%
    mutate(across(where(is.numeric), ~ replace_na(.x, -999)))  # prepare for machine learning
  
  if (forecast) {
    forecast_env <- env_cohorts %>%
      filter(brood_yr == max(data$brood_yr[data$age_group == dep_age]),
             year ==  max(data$year[data$age_group == dep_age])) %>%
      mutate(brood_yr = brood_yr + 1,
             year = year + 1)
    
    env_cohorts <- env_cohorts %>%
      bind_rows(forecast_env)
    
  }
  
  
  forecast_data <- data %>%
    filter(year == max(year)) %>%
    mutate(year = year + 1,
           returns = -999)
  
  # deal with deltas
  if (delta_returns) {
    ogret <- data$returns
    
    data <-  data %>%
      group_by(stock) %>%
      mutate(returns = returns - lag(returns, 1)) %>%
      ungroup()
    
    trim <- !is.na(data$returns)
    
    data <- data[trim,]
    
    ogret <- ogret[trim]
    
  }
  
  
  training <- data %>%
    filter(year < test_year)
  
  
  if (weight) {
    train_weights <-
      sqrt(training$returns) / sum(sqrt(training$returns))
    
  } else {
    train_weights <- rep(1, nrow(training)) / nrow(training)
    
  }
  
  # create cross validation splits
  
  
  annual_returns <-
    training %>% group_by(year) %>% nest() %>% ungroup() %>%
    arrange(year)
  
  salmon_rolling_origin <- rsample::rolling_origin(
    annual_returns,
    initial =  round(nrow(annual_returns) * initial_prop),
    assess = assess,
    cumulative = TRUE
  )
  
  
  # tune model
  
  if (use_years == FALSE) {
    training <- training %>%
      select(-year)
  }
  
  
  salmon_recipe <-
    recipe(returns ~ ., data = data) %>%
    {
      if (log_returns == TRUE) {
        step_log(., all_outcomes())
      } else {
        .
      }
    } %>%
    step_lag(returns, lag = 1:3, default = -999) %>% # generate lags going back three years
    {
      if (factor_years) {
        step_mutate(., year = factor(year,
                                     ordered = TRUE,
                                     levels = unique(data$year)))
      } else {
        .
      }
    } %>%
    step_dummy(all_nominal_predictors())
  
  
  stats <- training %>%
    group_by(stock) %>%
    summarise(mean = mean(returns, na.rm = TRUE),
              sd = sd(returns, na.rm = TRUE))
  
  training <- training %>%
    left_join(stats, by = "stock") %>%
    mutate(returns = (returns - mean) / sd) %>%
    select(-mean, -sd)
  
  data <- data %>%
    left_join(stats, by = "stock") %>%
    mutate(returns = (returns - mean) / sd)
  
  if (delta_returns) {
    ogret <- (ogret - data$mean) / data$sd
  }
  
  data <- data %>%
    select(-mean, -sd)
  
  prepped_salmon <-
    prep(salmon_recipe, data = training, retain = TRUE)
  
  baked_salmon <-
    bake(prepped_salmon, new_data = training %>% group_by(stock))
  
  # baked_salmon$stock <- training$stock
  if (run_wide) {
    baked_salmon <- make_wide_pinks(baked_salmon, training$stock)
    
    prepped_salmon <-
      prep(recipe(returns ~ ., data = baked_salmon), new_data = baked_salmon)
  }
  
  
  if (model_type == "rand_forest") {
    tune_grid <-
      parameters(min_n(range(2, 10)), mtry(), trees(range(200, 2000)))  %>%
      dials::finalize(mtry(), x = baked_salmon %>% select(-(returns)))
    
    ranger_grid <- grid_latin_hypercube(tune_grid, size = 30) %>%
      mutate(grid_row = 1:nrow(.))
    
    tune_grid <-
      tidyr::expand_grid(
        grid_row = 1:nrow(ranger_grid),
        id = unique(salmon_rolling_origin$id)
      ) %>%
      left_join(ranger_grid, by = "grid_row") %>%
      left_join(salmon_rolling_origin, by = "id")
    
  }
  
  if (model_type == "boost_tree") {
    tune_grid <-
      parameters(
        min_n(range(2, 3)),
        tree_depth(range(2, 15)),
        learn_rate(c(log10(.0001), log10(.04))),
        mtry(),
        loss_reduction(range = c(-10, 1.5)),
        sample_prop(range(.5, 1)),
        trees(range(100, 1000))
      ) %>%
      dials::finalize(mtry(), x = baked_salmon %>% select(-(1:2)))
    
    xgboost_grid <- grid_latin_hypercube(tune_grid, size = 40) %>%
      mutate(grid_row = 1:nrow(.)) #%>%
    # mutate(trees = trees)
    tune_grid <-
      tidyr::expand_grid(
        grid_row = 1:nrow(xgboost_grid),
        id = unique(salmon_rolling_origin$id)
      ) %>%
      left_join(xgboost_grid, by = "grid_row") %>%
      left_join(salmon_rolling_origin, by = "id")
    
  }
  
  
  tune_pars <- tune_grid[, colnames(tune_grid) != "id"] %>%
    as.list()
  
  tuning_fit <- future_pmap(
    tune_pars,
    tune_pink_salmon,
    salmon_recipe = salmon_recipe,
    model_type = model_type,
    log_returns = log_returns,
    run_wide = run_wide,
    .progress = TRUE,
    .options = furrr_options(seed = 42)
  )
  
  tune_grid$tuning_fit <- tuning_fit
  
  best_params <- tune_grid %>%
    select(-splits) %>%
    unnest(cols = tuning_fit)
  
  tune_vars <-
    colnames(best_params)[!colnames(best_params) %in% c(".pred", "observed", "id", "grid_row")]
  
  best_params <- best_params %>%
    group_by(across({
      {
        tune_vars
      }
    })) %>%
    # yardstick::mae(observed, .pred) %>%
    yardstick::rmse(observed, .pred) %>%
    ungroup()
  
  best <- best_params %>%
    filter(.estimate == min(.estimate)) %>%
    dplyr::slice(1)
  
  if (model_type == "rand_forest") {
    trained_model <-
      parsnip::rand_forest(
        mode = "regression",
        mtry = best$mtry,
        min_n = best$min_n,
        trees = best$trees
      ) %>%
      parsnip::set_engine("ranger",
                          importance = "none") %>%
      parsnip::fit(formula(prepped_salmon), data = baked_salmon)
    
    
  } # close model_type == "random_forest"
  if (model_type == "boost_tree") {
    trained_model <-
      parsnip::boost_tree(
        mode = "regression",
        mtry = best$mtry,
        min_n = best$min_n,
        learn_rate = best$learn_rate,
        tree_depth = best$tree_depth,
        loss_reduction = best$loss_reduction,
        trees = best$trees,
        sample_size = best$sample_size
      ) %>%
      parsnip::set_engine("xgboost") %>%
      parsnip::fit(formula(prepped_salmon), data = baked_salmon)
    
  }
  
  tmp_prep <-
    prep(salmon_recipe, data = data, retain = TRUE)
  
  baked_data <-
    bake(tmp_prep, new_data = data %>% group_by(stock))
  if (run_wide) {
    baked_data <- make_wide_pinks(baked_data, data$stock)
    
  }

  pred <-
    predict(trained_model, new_data = baked_data)
  
  if (delta_returns) {
    data$returns <- ogret
    
  }
  
  
  data$pred <-
    dplyr::case_when(log_returns == TRUE ~ exp(pred$.pred), TRUE ~  pred$.pred)

  data <- data %>%
    left_join(stats, by = "stock") %>%
    mutate(returns = (returns * sd) + mean,
           pred = (pred * sd) + mean) %>%
    select(-sd, -mean)
  
  if (delta_returns) {
    # somewhat convoluted reconstruction of raw forecasts from delta forecasts
    data$pred[data$year < test_year] <-
      data$returns[data$year < test_year] # first, replace training period predictions with raw returns
    
    # because these are pink salmon years aren't continuous across stocks, so setting up stock years
    data <-  data %>%
      group_by(stock) %>%
      mutate(stock_year = 1:length(returns)) %>%
      ungroup()
    
    for (s in unique(data$stock)) {
      # pull out the stock years in the 2:end of the test years
      years <-
        unique(data$stock_year[data$stock == s &
                                 data$year >= test_year])
      
      
      for (y in years) {
        data$pred[data$stock_year == y &
                    data$stock == s] <-
          pmax(0, data$pred[data$stock_year == (y - 1) &
                              data$stock == s] + data$pred[data$stock_year == y &
                                                             data$stock == s]) # convert delta predictions into raw predictions, building off the observed values in the last year of the training data
        
      } # close year loops
      
      
    } # close stock loop
    
    data <- data %>% 
      select(-stock_year)
    
  } # close if delta returns
  
  
  data <- data %>%
    mutate(
      split = dplyr::case_when(
        year < test_year ~ "training",
        year <= max(data$year) ~ "testing",
        TRUE ~ "forecast"
      ),
      returns = returns * scalar,
      pred = pred * scalar
    )
  
  data$returns[data$split == "forecast"] <-  NA
  #
    # data %>%
    #   filter(year == test_year) %>% 
    #   ggplot(aes(returns, pred, color = split)) +
    #   geom_abline(slope = 1, intercept = 0) +
    #   geom_point()
    # 
  if (produce == "summary") {
    out <- list(data = data, best_params = best_params)
  } else {
    out <-
      list(data = data, trained_model = trained_model)
    
  }
  
  # message(paste0(p_done," done"))
  
  return(out)
  
  
  
}