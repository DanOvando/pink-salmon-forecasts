fit_parametric_pinks <-
  function(this_stock,
           test_year,
           run_wide,
           delta_returns,
           data,
           chains = 4,
           cores = 4,
           warmup = 2500,
           iter = 5000,
           refresh = 0,
           model = "gam") {
    stock_data <- data %>%
      ungroup() %>%
      filter(stock == this_stock) %>%
      arrange(year)
    
    if (model == "gam") {
      if (delta_returns) {
        ogret <- stock_data$returns
        
        stock_data <- stock_data %>%
          mutate(returns = returns - lag(returns, 1))
        
        trim <- !is.na(stock_data$returns)
        
        stock_data <- stock_data[trim,]
        
        ogret <- ogret[trim]
        
        ogret <-
          ogret[-c(1:2)] # account for lags below, fix this if you add more filters
      }
      
      stock_data <- stock_data %>%
        # mutate(returns = returns - lag(returns,1)) %>%
        mutate(rl1 = lag(returns, 1),
               rl2 = lag(returns, 2)) %>%
        filter(!is.na(rl2)) %>%
        mutate(training = year < test_year)
      
      
      if (run_wide) {
        buddies <- data %>%
          ungroup() %>% 
          filter(region == unique(stock_data$region)) %>%
          select(year, subregion, returns) %>%
          group_by(subregion) %>%
          arrange(year) %>%
          mutate(rl1 = lag(returns, 1),
                 rl2 = lag(returns, 2)) %>%
          select(-returns) %>%
          filter(!is.na(rl2)) %>%
          pivot_longer(contains("rl")) %>%
          unite(subregion, name, remove = TRUE, col = "id") %>%
          pivot_wider(names_from = id, values_from = value)
        
        stock_data <- stock_data %>%
          select(year, returns, training, juvenile_index, mean_sst) %>%
          left_join(buddies, by = "year")
        
        gam_formula <-
          as.formula(paste0(
            "returns ~ s(year) + s(juvenile_index) + s(mean_sst) + ",
            paste(colnames(buddies %>% select(-year)), collapse = "+")
          ))
        
      } else {
        gam_formula <-
          as.formula("returns ~ s(year) +s(juvenile_index) + s(mean_sst) + s(rl1) + s(rl2)")
        
      }

      stock_data$year <- stock_data$year -2000
      
      if (delta_returns) {
        pink_gam <-
          stan_gamm4(
            gam_formula,
            data = stock_data %>% filter(training),
            family = gaussian(),
            chains = chains,
            cores = cores,
            warmup = warmup,
            iter = iter,
            refresh = refresh,
            prior_smooth = exponential(0.75, autoscale = FALSE)
          )
        
      } else{
        pink_gam <-
          stan_gamm4(
            gam_formula,
            data = stock_data %>% filter(training),
            family = Gamma(link = "log"),
            chains = chains,
            cores = cores,
            warmup = warmup,
            iter = iter,
            refresh = refresh,
            prior_smooth = exponential(rate = 0.75,autoscale = FALSE)
          )
        
      }
      
      
      stock_forecasts <-
        tidybayes::predicted_draws(pink_gam, newdata = stock_data, ndraws = 1000) %>%
        group_by(year) %>%
        summarise(
          returns = unique(returns),
          pred = mean(.prediction),
          median_predicted = median(.prediction)
        ) %>%
        mutate(training = year < test_year)
      
      if (delta_returns) {
        stock_forecasts$returns <- ogret
        
        
        stock_forecasts$pred[stock_forecasts$year < test_year] <-
          stock_forecasts$returns[stock_forecasts$year < test_year] # first, replace training period predictions with raw returns
        
        # because these are pink salmon years aren't continuous across stocks, so setting up stock years
        stock_forecasts <-  stock_forecasts %>%
          mutate(stock_year = 1:length(returns)) %>%
          ungroup()
        
        years <-
          unique(stock_forecasts$stock_year[stock_forecasts$year >= test_year])
        
        check <- stock_forecasts$pred
        for (y in years) {
          stock_forecasts$pred[stock_forecasts$stock_year == y] <-
            pmax(0, stock_forecasts$pred[stock_forecasts$stock_year == (y - 1)] + stock_forecasts$pred[stock_forecasts$stock_year == y]) # convert delta predictions into raw predictions, building off the observed values in the last year of the training data
          
        } # close year loops
        
        
        stock_forecasts <- stock_forecasts %>%
          select(-stock_year)
        
      } # close delta returns
      
      #
      # test %>%
      #   ggplot(aes(year, .prediction)) +
      #   stat_lineribbon()
      
      # pp <- posterior_predict(pink_gam, newdata = stock_data) %>%
      #   as_tibble() %>%
      #   map_df(mean)
      #
      # stock_data$pp <- as.numeric(pp[1,])
      # browser()
      #
      # stock_forecasts %>%
      #   ggplot() +
      #   geom_point(aes(year, returns, color = training)) +
      #   geom_line(aes(year, pred))
      
    } else if (model == "juvenile_cpue") {
      stock_data <- stock_data %>%
        mutate(training = year < test_year)
      
      juvenile_model <-
        lm(returns ~ juvenile_index, data = stock_data %>% filter(training))
      
      forecasts <- predict(juvenile_model, newdata = stock_data)
      
      stock_forecasts <- stock_data %>%
        mutate(
          pred = forecasts,
          delta_returns = FALSE,
          run_wide = FALSE,
          model_type = model
        )
      
      
    }
    return(stock_forecasts)
    
  }