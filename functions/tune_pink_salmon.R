tune_pink_salmon <-
  function(model_type = "rand_forest",
           salmon_recipe,
           log_returns = FALSE,
           run_wide = TRUE,
           ...) {
    arguments <- list(...)
    
    analysis_split <-   rsample::analysis(arguments$splits) %>%
      unnest(data) %>%
      group_by(stock) 
    
    assessment_split <-  rsample::assessment(arguments$splits) %>%
      unnest(data) %>%
      group_by(stock) 
    
    stats <- analysis_split %>% 
      group_by(stock) %>% 
      summarise(mean = mean(returns, na.rm = TRUE),
                sd = sd(returns, na.rm = TRUE))
    
    combo <- analysis_split %>%
      bind_rows(assessment_split) # doing this to deal with lags, no big deal since there are no calculations in the recipe
    
    combo <- combo %>% 
      left_join(stats, by = "stock") %>% 
      mutate(returns = (returns - mean) / sd) %>% 
      select(-mean, -sd)
    
    prepped_salmon <- prep(salmon_recipe, new_data = analysis_split)
    
    baked_data <- bake(prepped_salmon, new_data = combo %>% group_by(stock))
    
    # baked_data$stock <- combo$stock
    
    
    if (run_wide){
      
      baked_data <- make_wide_pinks(baked_data, combo$stock)
      
    }
    
    analysis_data <- baked_data %>%
      filter(year <= max(analysis_split$year))
    
    assessment_data <- baked_data %>%
      filter(year > max(analysis_split$year))
    
    prepped_salmon <- prep(recipe(returns ~ ., data = baked_data), new_data = baked_data)
    
    if (model_type == "rand_forest") {
      salmon_fits <-
        parsnip::rand_forest(
          mode = "regression",
          mtry = arguments$mtry,
          min_n = arguments$min_n,
          trees = arguments$trees
        ) %>%
        parsnip::set_engine("ranger",
                            importance = "none",
                            splitrule = arguments$splitrule) %>%
        parsnip::fit(formula(prepped_salmon), data = analysis_data)
      
    }
    if (model_type == "boost_tree") {
      salmon_fits <-
        parsnip::boost_tree(
          mode = "regression",
          mtry = arguments$mtry,
          min_n = arguments$min_n,
          learn_rate = arguments$learn_rate,
          tree_depth = arguments$tree_depth,
          loss_reduction = arguments$loss_reduction,
          trees = arguments$trees,
          sample_size = arguments$sample_size
        ) %>%
        parsnip::set_engine("xgboost", nthread = 1) %>%
        parsnip::fit(formula(prepped_salmon), data = analysis_data)
      
    }
    
    
    assessment_prediction <-
      predict(salmon_fits, new_data = assessment_data)
    
    if (log_returns) {
      assessment_prediction$.pred <- exp(assessment_prediction$.pred)
    }
    
    if (nrow(assessment_prediction) != nrow(assessment_split)) {
      stop("something really bad hapened in tune-salmon.R")
    }
    assessment_prediction$observed = assessment_split$returns
    
    return(assessment_prediction)
    
  } # close tune forest
