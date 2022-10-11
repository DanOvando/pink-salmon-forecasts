# make_pink_forecast
# early exploratory version of a script to manage machine learning forecasts of pink salmon returns
# in Kodiak, Prince William Sound, and Southeast Alaska.
# pulls in data on  pink salmon abundance from supplied spreadsheets
# collects data into consistent formats
# fits range of models for one-step-ahead forecasting


# setup -------------------------------------------------------------------


library(tidyverse) # kitchen sink of wrangling

library(lubridate) # wrangle dates

library(readr) # read in various formats

library(readxl) # read in xlsx

library(here) # construct robust relative paths

library(tidymodels) # kitchen sink of predictive modeling tools

library(furrr) # parallel processing for purrr

library(ranger) # random forests

library(rstanarm) # bayesian linear models

cores <- 8 # cores for parallel processing

functions <- list.files(here::here("functions")) # load functions

purrr::walk(functions, ~ source(here::here("functions", .x)))

future::plan(future::multisession, workers = cores) # set up parallel processing

results_path <- file.path("results", "v0.4") # location for results

if (!dir.exists(results_path)){
  dir.create(results_path, recursive = TRUE)
}

fit_forecasts <- TRUE # fit forecasts or load saved run

regionalize <- TRUE

theme_set(theme_minimal(base_size = 18)) # set theme for plots

# seak --------------------------------------------------------------------
# main file received from Andy Piston via email, file Southeast Alaska Pink Data for UW 2022
# Also covariates received from Sara Miller primarily in var2021_final


seak_raw <-
  readxl::read_xlsx(
    here("data", "Southeast Alaska Pink Data for UW 2022.xlsx"),
    sheet = "Data1",
    skip = 3,
    # three blank header rows for some reason
    n_max = 145 # weird things happen after 145 rows
  ) %>%
  janitor::clean_names()

# create harvest layer
seak_harvest <- seak_raw %>%
  select(year, sse, nsei, nseo) %>%
  pivot_longer(-year, values_to = "harvest", names_to = "subregion") %>%
  mutate(harvest = as.numeric(harvest)) %>%
  filter(!is.na(harvest))

# create escapement layer
seak_escape <- seak_raw %>%
  select(year, sse_index, nsei_index, nseo_index) %>%
  pivot_longer(-year, values_to = "escapement", names_to = "subregion") %>%
  mutate(escapement = as.numeric(escapement),
         subregion = str_remove_all(subregion, "_index")) %>% # clean up subregion name
  filter(!is.na(escapement))

# join harvest and escapement
seak <- seak_escape %>%
  left_join(seak_harvest, by = c("year", "subregion")) %>%
  mutate(region = "seak",
         returns = escapement + harvest) %>%
  mutate(across(c(escapement, harvest, returns),  ~ .x * 1e6)) # convert out of millions of salmon

seak_index <-   readxl::read_xlsx(here("data", "var2021_final.xlsx"),
                                  sheet = "var2021_final") %>%
  janitor::clean_names() %>%
  select(year, cpu_ecal, isti20_mjj) %>%
  rename(juvenile_index = cpu_ecal,
         juvenile_sst =  isti20_mjj) %>%
  arrange(year)
  
  seak_index_model_data <- seak_index %>% 
  select(year, juvenile_index) %>% 
  left_join(seak, by = "year") %>% 
  filter(!is.na(returns))

seak_index_model <-
  ranger(juvenile_index ~ returns + subregion, data = seak_index_model_data)


seak_index_prediction <- predict(seak_index_model, seak)

seak_index_impute <-
  data.frame(year = seak$year,
             predicted_index = seak_index_prediction$predictions) %>%
  group_by(year) %>%
  summarise(predicted_index = mean(predicted_index))


seak_index <- seak_index %>% 
  right_join(seak_index_impute, by = "year") %>% 
  arrange(year)

seak_index$juvenile_index[is.na(seak_index$juvenile_index)] <- seak_index$predicted_index[is.na(seak_index$juvenile_index)]


seak_index <-  seak_index %>% 
  select(-predicted_index) %>% 
  mutate(
    juvenile_index_l1 = lag(juvenile_index),
    juvenile_index_l2 = lag(juvenile_index, 2),
    juvenile_index_l3 = lag(juvenile_index, 3),
    juvenile_sst_l1 = lag(juvenile_sst),
    juvenile_sst_l2 = lag(juvenile_sst, 2),
    juvenile_sst_l3 = lag(juvenile_sst, 3)
  )
  # rename(year = j_year)



# kodiak ------------------------------------------------------------------
# received from Birch Foster via email, file "pink salmon data request C. Boatright"
#

kod_raw <-
  readxl::read_xlsx(
    here("data", "pink salmon data request C. Boatright.xlsx"),
    sheet = "Kodiak Wild pinks",
    skip = 1
  ) %>%
  janitor::clean_names()


# create harvest layer
kod_harvest <- kod_raw %>%
  select(year_1, afognak_h:kmah) %>%
  rename(year = year_1,
         kma = kmah) %>% # looks like the correct subregion name is kma, the h is supposed to show harvest but got lost because CamelCase
  pivot_longer(-year, names_to = "subregion", values_to = "harvest") %>%
  mutate(subregion = str_remove_all(subregion, "_h"))

# create escapement layer
kod_escape <- kod_raw %>%
  select(year_1:kmae) %>%
  rename(year = year_1,
         kma = kmae) %>% # looks like the correct subregion name is kma, the E is supposed to show escapement but got lost because CamelCase
  pivot_longer(-year, names_to = "subregion", values_to = "escapement") %>%
  mutate(subregion = str_remove_all(subregion, "_e"))



# join together kodiak data
kod <- kod_escape %>%
  left_join(kod_harvest, by = c("year", "subregion")) %>%
  mutate(region = "kod",
         returns = escapement + harvest)


# pws ---------------------------------------------------------------------
# received from Jennifer Morella via email, file " pink forecast request"
# CCP stands for  commercial common property fishery
# HCR is hatchery cost recovery
# Stream index is then "escapement", to a point
# Total is then total return

pws_raw <-
  readxl::read_xlsx(here("data", "pink forecast request.xlsx"),
                    sheet = "Sheet1") %>%
  janitor::clean_names()

pws <- pws_raw %>%
  rename_with( ~ "escapement", starts_with("x134")) %>%
  rename(year = run_year) %>%
  mutate(across(starts_with("wild_"), as.numeric)) %>% # some problems with some manual NAs instead of blanks in raw file
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
  mutate(
    harvest =
      wild_ccp +
      wild_hcr  +
      wild_hatchery_rack_brood_and_post_brood_sales,
    # per documentation wild harvests
    returns = escapement + harvest
  ) %>%
  select(year, harvest, escapement, returns) %>%
  mutate(region = "pws",
         subregion = "pws")


# all together ------------------------------------------------------------
# pull all the pink data together

pinks <- bind_rows(seak, kod, pws) %>%
  mutate(cycle = ifelse(year %% 2 == 0, "even", "odd")) %>%
  mutate(returns = returns / 1e6,
         escapement = escapement / 1e6,
         harvest = harvest / 1e6 )  # turning into millions of salmon

if (regionalize){
pinks <- pinks %>% 
  group_by(year, region, cycle) %>% 
  summarise(returns = sum(returns),
            escapement = sum(escapement), 
            harvest = sum(harvest)) %>% 
  mutate(subregion = region) %>% 
  arrange(region, year)
}

# add in additional data

  pinks <- pinks %>%
  left_join(seak_index, by = "year") 

# quick exploratory plot
pinks %>%
  ggplot(aes(year, returns, fill = subregion)) +
  geom_col() +
  facet_wrap(~ region)

pinks %>% 
  ggplot(aes(juvenile_index, returns, color = subregion)) + 
  geom_point() + 
  geom_smooth(method ="lm") +
  facet_wrap(~region)



# explore -----------------------------------------------------------------


pinks_tmp <- pinks %>% 
  select(-contains("juvenile"))

foo <- function(x, y) {
  # generate arbitrary numbers of lags of returns per stock
  y %>%
    group_by(cycle, subregion, region) %>%
    arrange(year) %>%
    mutate("lag_{{x}}" := dplyr::lag(returns, x, order_by = year)) %>%
    ungroup() %>%
    select(contains("lag_"))
  
}

pinks_tmp <- pinks_tmp %>%
  group_by(cycle, subregion, region) %>%
  arrange(year) %>%
  bind_cols(map_dfc(1:10, foo, y = pinks_tmp)) %>% # bind on columns of lags per stock
  rename_with( ~ str_remove_all(.x, "L"), starts_with("lag_")) %>% # remove L from number indicating how many lags back we are
left_join(seak_index, by = "year") %>% 
  mutate(across(where(is.numeric), ~ replace_na(.x, -999)))  # prepare for machine learning
  

# run a quick exploratory random forest
exploratory_model <-
  ranger(
    returns ~ .,
    data = pinks_tmp %>% select(returns, contains("region"), year, contains("lag_"), contains("juvenile")),
    importance = "impurity_corrected"
  )

vip::vip(exploratory_model) # variable importance plot of random forest returns model


# fit models --------------------------------------------------------------


data <- pinks %>%
  select(-escapement, -harvest)

# generate grid of candidate models for one-year-ahead forecasts from 2010 to most recent year, running both boosted regression trees and random forests, with and without "wide" format where every stock is a predictor of all other stocks.
model_grid <-
  tidyr::expand_grid(
    test_year = 2010:max(pinks$year),
    model_type = c("rand_forest"),
    delta_returns = c(FALSE,TRUE),
    run_wide = c(TRUE, FALSE)
  )

pinks$stock <- paste(pinks$subregion, pinks$cycle, sep = "_")

write_rds(pinks, file = file.path(results_path,"pink_data.rds"))
# fit or load forecast models
if (fit_forecasts) {
  fits <- model_grid %>%
    mutate(fit = pmap(
      list(
        test_year = test_year,
        model_type = model_type,
        run_wide = run_wide,
        delta_returns = delta_returns
      ),
      fit_pink_model,
      # map over candidate models
      data = data,
      log_returns = FALSE # log returns
    ))
  
  write_rds(fits, file = file.path(results_path, "fits.rds"))
  
  
  
  # fit GAMS

  sfpp <- safely(fit_parametric_pinks)
  
  
  gams <-
    expand_grid(
      this_stock = unique(pinks$stock),
      delta_returns = c(FALSE,TRUE),
      run_wide = c(TRUE,FALSE),
      test_year = 2010:max(pinks$year)
    ) 
  
  gams <- gams %>%
    mutate(fit = pmap(
      list(
        this_stock = this_stock,
        run_wide = run_wide,
        delta_returns = delta_returns,
        test_year = test_year
      ),
      sfpp,
      data = pinks,
      refresh = 0
    ))
  
  write_rds(gams, file = file.path(results_path, "gam_fits.rds"))
  
  # 
  # i <- sample(1:nrow(gams),1)
  # gams$fit[[i]]$result %>%
  #   ggplot() +
  #   geom_point(aes(year, returns, color = training)) +
  #   geom_line(aes(year, pred))
  # 
  
  juvenile_cpue <-
    expand_grid(this_stock = unique(pinks$stock),
                test_year = 2010:max(pinks$year)) %>%
    mutate(fit = pmap(
      list(this_stock = this_stock, test_year = test_year),
      fit_parametric_pinks,
      data = pinks,
      model = "juvenile_cpue"
    )) %>% 
    select(-this_stock)
  
  write_rds(juvenile_cpue, file = file.path(results_path, "juvenile_cpue_fits.rds"))
  
  
  
} else {
  fits <- read_rds(file = file.path(results_path, "fits.rds"))
  
  gams <- read_rds(file = file.path(results_path, "gam_fits.rds"))
  
  juvenile_cpue <- read_rds(file = file.path(results_path, "juvenile_cpue_fits.rds"))
  
}

region_lookup <- pinks %>% 
  select(stock, region) %>% 
  unique()

gams_worked <- map_lgl(map(gams$fit, "error"), is.null)

gams <- gams %>% 
  filter(gams_worked) %>% 
  mutate(fit = map(fit, "result")) %>% 
  unnest(cols = fit) %>% 
  mutate(model_type = "gam") %>% 
  rename(stock = this_stock) %>% 
  mutate(split = ifelse(year >= test_year, "testing", "training")) %>% 
  left_join(region_lookup, by = "stock") %>% 
  filter(!(run_wide == TRUE & delta_returns == FALSE))

juvenile_cpue <- juvenile_cpue %>% 
  unnest(cols = fit) %>% 
  mutate(split = ifelse(year >= test_year, "testing", "training")) 

# unpack  and prepare results
results <- fits %>%
  mutate(tmp = map(fit, "data")) %>%
  unnest(cols = tmp) %>%
  mutate(year = lubridate::year(year))



# generate benchmark model which is just a lag(1) model
lag_benchmark <- results %>%
  filter(model_type == "rand_forest") %>%
  group_by(test_year, region, cycle, stock, run_wide, delta_returns) %>%
  mutate(pred = lag(returns, 1),
         model_type = "lag(1)")

# add lag(1) model to candidate models
results <- results %>%
  bind_rows(lag_benchmark) %>% 
  bind_rows(gams) %>% 
  bind_rows(juvenile_cpue)


# exploratory graph of forecast performance by model type and testing vs. training split
results %>%
  filter(year == test_year) %>% 
  ggplot(aes(returns, pmin(5 * returns,pred), color = model_type)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(delta_returns~run_wide, labeller = label_both)

results %>%
  unite(col = "model",model_type, run_wide, delta_returns) %>% 
  filter(year == test_year) %>% 
  ggplot(aes(pred, returns, color = model, fill = model)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = 1) +
  geom_smooth(method = "lm", alpha = 0.5) + 
  coord_cartesian(xlim = c(0, 100)) + 
  scale_x_continuous(name = "Predicted") + 
  scale_y_continuous(name  = "Observed")



results %>% 
  filter(run_wide == FALSE,model_type != "lag(1)") %>% 
  separate(stock, sep = "_", into = c("subregion", "cycle")) %>% 
  filter(year == test_year) %>% 
  ggplot() + 
  geom_point(aes(year, returns)) + 
  geom_line(aes(year, pmin(5*returns,pred), color = model_type, linetype = delta_returns)) +
  facet_wrap(~subregion, scales = "free_y")


results %>% 
  filter(year == test_year) %>% 
  group_by(region, model_type, run_wide, year, delta_returns) %>% 
  summarise(returns = sum(returns),
            pred = sum(pred)) %>% 
  ggplot() + 
  geom_point(aes(year, returns)) + 
  geom_line(aes(year, pred, color = delta_returns, linetype = run_wide)) +
  facet_grid(model_type~region, scales = "free_y")


# calculate performance metrics for one-step-ahead forecast
performance <- results %>%
  filter(year == test_year) %>%
  group_by(model_type, stock, region, run_wide, delta_returns) %>%
  summarise(
    rsq = yardstick::rsq_vec(truth = returns, estimate = pred),
    mape =  yardstick::mape_vec(truth = returns, estimate = pred),
    rmse =  yardstick::rmse_vec(truth = returns, estimate = pred)
  )


# quick plot of various performance metrics
performance %>%
  pivot_longer(
    cols = c("rsq", "mape", "rmse"),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(metric == "rmse") %>% 
  ggplot(aes(model_type, value, fill = delta_returns)) +
  geom_boxplot() +
  facet_grid(run_wide ~ metric, scales = "free_y") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) 




# estimate average change in rsq per model type, taking into account clustered nature of stocks per region
rsq_model <-
  stan_glmer(
    rsq ~ model_type + (stock |
                          region) - 1,
    data = performance %>% filter(run_wide == FALSE),
    chains = 4
  )

plot(rsq_model, regex_pars = "model_type")


# estimate average change in mape per model type, taking into account clustered nature of stocks per region

rmse_model <-
  stan_glmer(
    (rsq) ~ model_type:delta_returns + (stock |
                                      region) - 1,
    data = performance %>%  filter(run_wide == FALSE),
    chains = 4,
    cores = 4
  )

rmse_pars <- as.data.frame(rmse_model) %>% 
  as_tibble()

model_type_loc <- grepl("model_type", colnames(rmse_pars))

model_type_results <- rmse_pars[,model_type_loc] %>% 
  pivot_longer(everything())


model_type_results %>% 
  ggplot(aes(y = name, x = value)) + 
  ggdist::stat_halfeye()




