make_wide_pinks <- function(x,stocks){
  
  covs <- x %>% 
    mutate(stock = stocks) %>% 
    select(stock,year, contains("lag")) %>% 
    pivot_longer(contains("lag"), names_to = "lag", values_to = "returns") %>% 
    mutate(lag = str_remove_all(lag, "_returns")) %>% 
    mutate(id = paste(lag, stock, sep = "_")) %>% 
    select(year, id, returns) %>% 
    pivot_wider(names_from = id, values_from = returns,values_fill = -999)
  
  
  z <- x %>% 
    select(-contains("lag_")) %>% 
    left_join(covs, by = "year")
  
  
  return(z)
  
}