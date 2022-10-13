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

library(sf)

library(rerddap)

cores <- 8 # cores for parallel processing

functions <- list.files(here::here("functions")) # load functions

purrr::walk(functions, ~ source(here::here("functions", .x)))

future::plan(future::multisession, workers = cores) # set up parallel processing

results_path <- file.path("results", "v0.5") # location for results

if (!dir.exists(results_path)){
  dir.create(results_path, recursive = TRUE)
}


get_environmental_data <- TRUE

fit_forecasts <- TRUE # fit forecasts or load saved run

regionalize <- FALSE

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



# add in environmental data -----------------------------------------------


title <- c(
  "  ICOADS, 1-degree, Enhanced, Monthly, 1960-present
",
"HadISST Sea Ice Component, 1°, Global, Monthly, 1870-present
",
"HadISST Average Sea Surface Temperature, 1°, Global, Monthly, 1870-present",
"TAO/TRITON, RAMA, and PIRATA Buoys, Daily, 1987-present, Salinity
",
"TAO/TRITON, RAMA, and PIRATA Buoys, Monthly, 1987-present, Salinity",
"Aquarius Sea Surface Salinity, L3 SMI, Version 5, 1.0°, Global, 2011-2015, Daily
"
)

type = c("SLP",
         "Sea Ice",
         "SST",
         "Salinity",
         "Salinity",
         "Salinity")

gridded <- c("Y",
             "Y",
             "Y",
             "N",
             "N",
             "Y")

gridded <- gridded == "Y"

erddap_data <- tibble(
  title = title,
  type = type ,
  id =
    c(
      "esrlIcoads1ge",
      "erdHadISSTIce",
      "erdHadISST",
      "pmelTaoDyS",
      "pmelTaoMonS",
      "jplAquariusSSSDailyV5"
    ),
  gridded = gridded
) %>% 
  filter(type == "SST")

if (get_environmental_data == TRUE) {
  
  kodiak_pws_erddap <- map2(
    erddap_data$id,
    erddap_data$gridded,
    query_erddap,
    min_lat = 53.18,
    max_lat = 60,
    min_lon = -159,
    max_lon = -143,
    min_year = min(pinks$year) - 10,
    max_year = max(pinks$year),
    stride = 2
  ) %>%
    set_names(erddap_data$id)
  
  
  seak_erddap <- map2(
    erddap_data$id,
    erddap_data$gridded,
    query_erddap,
    min_lat = 50,
    max_lat = 59,
    min_lon = -142,
    max_lon = -131,
    min_year = min(pinks$year) - 10,
    max_year = max(pinks$year),
    stride = 2
  ) %>%
    set_names(erddap_data$id)
  
  
  
  land = rnaturalearth::ne_download(category = "physical", type = "land",scale = 110, returnclass = "sf") # get borderless land
  
  a = ed_search(query = 'bathymetry', which = "grid")
  
  
  
  bath_dat <- info('srtm30plus_v11_bathy')
  
  (bathy <- griddap(bath_dat,
                    latitude = c(50, 62),
                    longitude = c(-159, -130),
                    stride = 5
  ))
  
  
  
  shallow <- bathy$data %>%
    filter(!is.na(elev) & elev > -500) %>% 
    mutate(latitude = plyr::round_any(latitude, .5),
           longitude = plyr::round_any(longitude, .5)) %>% 
    group_by(latitude, longitude) %>% 
    summarise(elev= mean(elev, na.rm = TRUE)) %>% 
    ungroup()
  
  
  shallow <- st_as_sf(shallow,
                      coords = c("longitude", "latitude"),
                      remove = FALSE, 
                      crs = sf::st_crs(kodiak_pws_erddap$erdHadISST$data))
  
  
  
  
  shallow_kodiak_pws <- sf::st_join(kodiak_pws_erddap$erdHadISST$data %>% filter(!is.na(sst)), shallow) %>% 
    filter(!is.na(elev)) %>% 
    mutate(month = lubridate::month(time),
           year = lubridate::year(time)) %>% 
    filter(!(lat > 55 & lon < -155))
  
  
  
  shallow_seak <- sf::st_join(seak_erddap$erdHadISST$data %>% filter(!is.na(sst)), shallow) %>% 
    filter(!is.na(elev)) %>% 
    mutate(month = lubridate::month(time),
           year = lubridate::year(time))  
  
  
  ggplot() + 
    geom_sf(data =   shallow_kodiak_pws %>% mutate(region = "kod/pws"), aes(color = region)) +
    geom_sf(data =   shallow_seak %>% mutate(region = "seak"), aes(color = region)) +
    geom_sf(data = land) + 
  coord_sf(xlim = c(-165,-125), ylim = c(40, 60))

  kod_pws_sst <- shallow_kodiak_pws  %>%
    filter(month %in% c(9:12)) %>%
    group_by(year) %>%
    summarise(
      mean_sst = mean(sst),
      max_sst = max(sst),
      min_sst = min(sst)
    ) %>%
    sf::st_drop_geometry()
  
  seak_sst <- shallow_seak  %>%
    filter(month %in% c(9:12)) %>%
    group_by(year) %>%
    summarise(
      mean_sst = mean(sst),
      max_sst = max(sst),
      min_sst = min(sst)
    ) %>%
    sf::st_drop_geometry()
  
  seak_sst %>% 
    ggplot(aes(year, min_sst)) + 
    geom_point()
  
  
  pdo <-
    read_table(
      "http://research.jisao.washington.edu/pdo/PDO.latest",
      na = c("-99.99", "99.99", '-99',"-9,90"),
      skip = 29,
      n_max = lubridate::year(Sys.time()) - 1900,
      col_names = c("year", 1:12)
    ) %>%
    slice(-1) %>% 
    gather(month, pdo,-year) %>%
    mutate(month = as.double(month),
           date = lubridate::ymd(paste(year, month, '01', sep = '-'))) %>%
    filter(month %in% (5:8)) %>% 
    mutate(year = str_replace_all(year, "\\D",'') %>%as.numeric(),
           pdo = as.numeric(pdo)) %>% 
    group_by(year) %>%
    summarise(env_pdo = mean(pdo, na.rm = TRUE)) 
  
  write_rds(pdo, here("data", "pdo.rds"))
  
  write_rds(seak_sst, here("data", "seak_sst.rds"))
  
  write_rds(kod_pws_sst, here("data", "kod_pws_sst.rds"))
  
} else {
  
  seak_sst <-  read_rds(here("data", "seak_sst.rds"))
  
  kod_pws_sst <- read_rds(here("data", "kod_pws_sst.rds"))
  
  pdo <- read_rds(here("data", "pdo.rds"))
  
  
}


# for now treat kod and PWS as same region for the purposes of juvenile conditions
sst <- seak_sst %>% 
  mutate(region = "seak") %>% 
  bind_rows(kod_pws_sst %>% mutate(region = "kod")) %>% 
  bind_rows(kod_pws_sst %>% mutate(region = "pws")) %>% 
  select(year, region, mean_sst)

sst %>% 
  ggplot(aes(year, mean_sst, color = region)) + 
  geom_line()
  


# add in hatchery ---------------------------------------------------------

kodiak_hatchery <-
  readxl::read_xlsx(
    here("data", "pink salmon data request C. Boatright.xlsx"),
    sheet = "Kitoi Hatchery",
    skip = 2
  ) %>%
  janitor::clean_names() %>% 
  select(year_2, number) %>%  # release yeear, pink fry released in numbers
  rename(year = year_2, hatchery_numbers = number) %>% 
  mutate(hatchery_numbers = hatchery_numbers / 1e6) %>% 
  mutate(region = "kod")# convert to millions of fish

kodiak_hatchery %>% 
  ggplot(aes(year, hatchery_numbers)) + 
  geom_line()

pws_hatchery <-   readxl::read_xlsx(
  here("data", "Appendix D1.xlsx"),
  sheet = "D1",
  skip = 3,
  n_max = 23
) %>% 
  janitor::clean_names() %>% 
  mutate(pws_hatchery_returns = rowSums(across(!starts_with("year")))) %>% 
  select(year, pws_hatchery_returns)

npafc <- readxl::read_xls(
  here("data", "NPAFC_Hatchery_Rel_Stat_Web_21June2022.xls"),
  sheet = "Hatchery Release",
  skip = 1
) %>% 
  janitor::clean_names() %>% 
  pivot_longer(starts_with("x"), names_to = "year", values_to = "number_released", names_prefix = "x", names_transform = list(year = as.integer)) %>% 
  mutate(number_released = tidyr::replace_na(number_released, 0))

alaska_hatchery_release <- npafc %>% 
  filter(whole_country_province_state == "Alaska") %>% 
  filter(species == "Pink", 
         reporting_area == "Southeast")


alaska_hatchery_release %>% 
  ggplot(aes(year, number_released, color = reporting_area)) + 
  geom_line() + 
  facet_wrap(~species)


hatchery_production <- kodiak_hatchery

# put together final dataset ----------------------------------------------



pinks <- pinks %>% 
  left_join(sst %>% mutate(year = year + 1 ), by = c("region", "year")) %>%  # track sea surface temperature in the year prior to returning
  left_join(hatchery_production %>% mutate(year = year ) %>% select(-region), by = c("year")) %>% 
  left_join(pws_hatchery %>%  mutate(year = year + 1 ), by = "year") %>% 
  ungroup() %>% 
  left_join(alaska_hatchery_release %>% select(year, number_released) %>% mutate(year = year ), by = c("year"))


pinks %>% 
  ggplot(aes(number_released / 1e6, returns, color = region, fill = region)) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~cycle) + 
  scale_x_continuous(name = "Millions of Pink Salmon Released The Year Before") + 
  scale_y_continuous(name = "Millions of Wild Pink Returns")

pinks %>% 
  ggplot(aes(year, returns, color = region, fill = region)) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~cycle)


pinks %>% 
  ggplot(aes(pws_hatchery_returns, returns, color = region, fill = region)) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~cycle)


pinks %>% 
  ggplot(aes(hatchery_numbers, returns, color = region, fill = region)) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~cycle)

pinks %>% 
  ggplot(aes(mean_sst, returns, color = region, fill = region)) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~cycle)

pinks %>% 
  ggplot(aes(year, returns, color = region)) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~cycle)



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
  bind_cols(map_dfc(1:2, foo, y = pinks_tmp)) %>% # bind on columns of lags per stock
  rename_with( ~ str_remove_all(.x, "L"), starts_with("lag_")) %>% # remove L from number indicating how many lags back we are
left_join(seak_index, by = "year") %>% 
  mutate(across(where(is.numeric), ~ replace_na(.x, -999)))  # prepare for machine learning
  

# run a quick exploratory random forest
exploratory_model <-
  ranger(
    returns ~ .,
    data = pinks_tmp %>% select(returns, contains("region"), year, contains("lag_"), contains("juvenile"), mean_sst, hatchery_numbers),
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




