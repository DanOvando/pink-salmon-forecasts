####### make pink salmon forecasting model ########


# libraries -----------------------------------------------------------------
set.seed(42)
library(here)
library(tidyverse)
library(tidymodels)
library(ranger)
library(xgboost)
library(readr)
library(stringr)
library(readxl)
library(forcats)
library(janitor)
library(tidyxl)
library(unpivotr)
library(corrr)
library(RcppRoll)
library(furrr)

functions <- list.files(here::here("functions"))

purrr::walk(functions, ~ source(here::here("functions", .x)))


# OPTIONS -----------------------------------------------------------------

run <- "v0.1"

run_forecast <- TRUE

min_years <- 10 # minimum number of years with returns needed for inclusion

cores <- 1 #parallel::detectCores()-4

# SETUP -------------------------------------------------------------------


future::plan(future::multisession, workers = cores)


results_dir <- here("results",run)

if (!dir.exists(results_dir)) {
  
  dir.create(results_dir, recursive = TRUE)
  
  dir.create(file.path(results_dir,"figs"), recursive = TRUE)
  
}

# load data ---------------------------------------------------------------

# load RAM data

if (file.exists(here("data", "ram.zip")) == FALSE) {
  download.file(
    "https://zenodo.org/record/4824192/files/RAMLDB%20v4.495.zip?download=1",
    destfile = here("data", "ram.zip"),
    mode = "wb"
  )
  
  unzip(here("data", "ram.zip"), exdir = here("data", "ram")) # unzip = 'unzip' needed for windows
}

ram_files <-
  list.files(here("data", "ram", "R Data"), recursive = TRUE)

ram_files <- ram_files[str_detect(ram_files, ".RData")]

load(here("data", "ram", "R Data", ram_files[1]))

stock <- stock %>%
  left_join(area, by = "areaid")

pink_stocks <- stock %>%
  filter(str_detect(tolower(scientificname), "oncorhynchus") &
           str_detect(tolower(commonname), "pink"))

# ss biomass

ram_ss_biomass <- ssb.data %>%
  mutate(year = rownames(.) %>% as.integer()) %>%
  as_tibble() %>%
  gather(stockid, ss_biomass,-year)

pink_ram <- ram_ss_biomass %>%
  filter(stockid %in% unique(pink_stocks$stockid)) %>% 
  mutate(cohort = ifelse((year %% 2) == 0, "Even", "Odd")) %>% 
  rename(returns = ss_biomass)


write_csv(pink_stocks, file = file.path(results_dir,"pink_stocks.csv"))

pink_ram_plot <- pink_ram %>%
  group_by(stockid) %>%
  mutate(
    has = any(!is.na(returns)),
    ssb = scale(returns),
    cohort = ifelse((year %% 2) == 0, "Even", "Odd")
  ) %>%
  ungroup() %>%
  filter(has,
         year > 1960, !is.na(ssb)) %>%
  ggplot(aes(year, ssb, color = stockid)) +
  geom_line(show.legend = FALSE, alpha = 0.5) +
  facet_wrap( ~ cohort) + 
  scale_x_continuous(name = "Year") + 
  scale_y_continuous("Centered and Scaled Spawning Biomass") + labs(caption = "Spawning biomass of pink salmon stocks in RAM over time. Even/Odd refers to stocks spawning in even or odd years")

pink_ram_plot


# load Ruggerone data

ruggerone_data <-
  tidyxl::xlsx_cells(
    here(
      "data",
      "mcf210023-sup-0001-tables1-s24 6.00.07 PM.xlsx"
    ),
    sheet = "ST 1-4 Nat-orig ret (nos)52-15"
  )


# find blank columns to indicate where the spaces are between the tables
# https://nacnudus.github.io/spreadsheet-munging-strategies/small-multiples-with-all-headers-present-for-each-multiple.html

seperators <- ruggerone_data %>%
  group_by(col) %>%
  summarise(is_sep = all(is_blank == TRUE)) %>%
  filter(is_sep == TRUE)

species <- ruggerone_data %>%
  filter(str_detect(character, "Natural-origin"))

exent <-
  ruggerone_data$local_format_id[ruggerone_data$row == unique(species$row)]

subtables <-
  ruggerone_data %>%
  filter(str_detect(character, "Year")) %>%
  select(row, col)

ruggerone_data <- ruggerone_data %>%
  filter(!col %in% seperators$col,
         col < seperators$col[3])


partitions <- unpivotr::partition(ruggerone_data, subtables) %>%
  mutate(species = species$character)


foo <- function(cells) {
  cells %>%
    unpivotr::behead("N", header) %>%
    select(row, data_type, header, numeric) %>%
    unpivotr::spatter(header) %>%
    select(-row, -Total) %>%
    janitor::clean_names() %>%
    gather(region, returns, -year)
}

rugg_salmon = partitions %>%
  mutate(cells = map(cells, foo)) %>%
  unnest(cols = cells) %>%
  select(-corner_row,-corner_col) %>%
  mutate(species = str_trim(tolower(
    str_replace_all(species, "(Natural-origin)|(Salmon)", '')
  ))) %>%
  filter(species != "sockeye")

ping_rugg <- partitions %>%
  mutate(cells = map(cells, foo)) %>%
  unnest(cols = cells) %>%
  select(-corner_row,-corner_col) %>%
  mutate(species = str_trim(tolower(
    str_replace_all(species, "(Natural-origin)|(Salmon)", '')
  ))) %>% 
  mutate(cohort = ifelse((year %% 2) == 0, "Even", "Odd"))



pink_rugg_plot <- ping_rugg %>%
  group_by(species, region) %>%
  mutate(
    returns = scale(returns)
  ) %>%
  ungroup() %>%
  filter(year > 1960) %>%
  ggplot(aes(year, returns, color = region)) +
  geom_line(show.legend = FALSE, alpha = 0.5) +
  facet_wrap( ~ cohort)


# load NPAFC data

# problems with these data, a bunch of characters hidden in the numbers, skipping processing this for now
npafc_data <- readxl::read_xlsx(here('data',"NPAFC_Catch_Stat_Web_6August2021.xlsx"), skip = 1) #%>% 
  # pivot_longer(matches("^\\d"), names_to = "year", values_to = "value")


# plot data ---------------------------------------------------------------

year_coverage_plot <- pink_ram %>% 
  filter(!is.na(returns)) %>% 
  group_by(stockid) %>% 
  summarise(ny = n_distinct(year)) %>% 
  ggplot(aes(ny)) + 
  geom_histogram()
  
year_coverage_plot

# very little cross cohort signal but not none

even_odd_plot <- pink_ram %>% 
  group_by(stockid) %>% 
  mutate(lag_returns = lag(returns)) %>% 
  ggplot(aes(returns, lag_returns, color = cohort)) + 
  geom_point()


# much clearer lag signal within cohort, demeaned correlation breaks down pretty quickly after 1 year, but a lag 2 might be useful. need system specific intercepts

lag_1_plot <- pink_ram %>% 
  group_by(stockid, cohort) %>% 
  mutate(returns = returns - mean(returns, na.rm = TRUE)) %>% 
  mutate(lag_returns = lag(returns,1)) %>% 
  ggplot(aes(returns, lag_returns, color = cohort)) + 
  geom_point() + 
  facet_wrap(~cohort)


# prepare data ------------------------------------------------------------


pink_data <- pink_ram %>% 
  unite("stock",c(stockid,cohort), sep = ":", remove = TRUE) %>% 
  filter(!is.na(returns)) %>% 
  group_by(stock) %>% 
  mutate(year_gap = year - lag(year)) %>% 
  filter(!is.na(year_gap)) %>% # need at least two years of data
  mutate(has_gap = any(year_gap > 2)) %>% 
  ungroup() %>% 
  filter(!has_gap) %>% 
  group_by(stock) %>% 
  mutate(ny = n_distinct(year)) %>% 
  filter(ny >= min_years) %>% 
  ungroup() %>% 
  select(-ny, -has_gap,-year_gap) %>% 
  group_by(stock)


# generate splits

pink_training <- pink_data


# pre-processing


pink_recipe <- recipe(returns ~ ., data = pink_data) %>% 
  step_lag(returns, lag = 1:3)


# fit models --------------------------------------------------------------

if (run_forecast){
  
  pink_forecasts <- expand_grid(
    test_year = 1990:max(pink_data$year),
    model_type = c("boost_tree", "rand_forest")
  )

  pink_forecasts <- pink_forecasts %>%
    mutate(forecast = map2(test_year, model_type, fit_pink_model, data = pink_data))
  
write_rds(pink_forecasts, file = file.path(results_dir,"pink_forecasts.rds"))

} else {
  
  pink_forecasts <- read_rds(file = file.path(results_dir,"pink_forecasts.rds"))
  
}

# evaluate models ---------------------------------------------------------

pink_results <- pink_forecasts %>% 
  mutate(tmp = map(forecast,"data")) %>% 
  select(-forecast) %>% 
  unnest(cols = tmp) %>% 
  arrange(stock) %>% 
  ungroup() %>% 
  mutate(year = lubridate::year(year)) %>% 
  rename(forecast = pred)

null_results <- pink_results %>% 
  filter(model_type == unique(pink_forecasts$model_type)[1]) %>% 
  group_by(stock, test_year) %>% 
  arrange(stock,year) %>% 
  mutate(forecast = lag(returns,1),
         model_type = "null_model")

pink_results <- pink_results %>% 
  bind_rows(null_results)


pink_results %>% 
  ggplot(aes(returns, forecast)) + 
  geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = 0.75) + 
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(model_type~split)


pink_results %>% 
  filter(split == "training" | year == test_year) %>% 
  ggplot(aes(returns / 1e6, forecast / 1e6)) + 
  geom_abline(slope = 1, intercept = 0) +
  geom_jitter(alpha = 0.5) + 
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(model_type~split) + 
  scale_x_continuous("Observed SSB Returns (millions MT)") + 
  scale_y_continuous("Predicted SSB Returns (millions MT)") + 
  labs(caption = "Training are first to the training data, Testing are out-of-sample forcasts.")


pink_performance_summary <- pink_results %>%
  filter(year == test_year | split == "training",!is.na(forecast)) %>%
  group_by(model_type, split) %>%
  summarise(
    rmse =  yardstick::rmse_vec(truth = returns, estimate = forecast),
    r2 = yardstick::rsq_vec(truth = returns, estimate = forecast),
    mape = yardstick::mape_vec(truth = returns, estimate = forecast)
  )

  
pink_performance_summary %>% 
  ggplot(aes(split, rmse, fill = model_type)) + 
  geom_col(position = "dodge") + 
  scale_y_continuous(name = "RMSE")


pink_performance_summary %>% 
  ggplot(aes(split, r2, fill = model_type)) + 
  geom_col(position = "dodge") + 
  scale_y_continuous(name = "R2")

