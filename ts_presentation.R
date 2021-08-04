library(fpp3)
library(slider)
library(broom)

### Mortality data ###

mortality <- readr::read_csv("mortality.csv") %>%
  mutate(week = yearweek(date), year = year(date)) %>% 
  as_tsibble(index = week, key = c(jurisdiction, cause)) %>%
  select(-date) 
mortality

# Time series plot
mortality %>%
  filter(jurisdiction == "United States") %>%
  autoplot(deaths)

### Population data ###

population <- readr::read_csv("population.csv")

population %>%
  as_tsibble(index = year, key = jurisdiction) %>%
  filter(jurisdiction == "United States") %>%
  autoplot(population)

range(population$year) 
range(mortality$week) 

# Population data covers 2013 - 2019, mortality data covers 2014 - 2021. We'll
# forecast the population using an AR model to get population estimates for
# 2020 and 2021. 

population_ar <- population %>%
  as_tsibble(index = year, key = jurisdiction) %>%
  model(AR(population ~ order(1))) %>%
  forecast(h = 2) %>%
  as_tsibble(index = year, key = jurisdiction) %>%
  select(jurisdiction, year, .mean) %>%
  rename(population = .mean) %>%
  bind_rows(population) %>%
  as_tibble()

# Join the mortality and population data, and use to calculate deaths per 
# 100,000 people (field: deaths_pc). NOTE: The population values are updated 
# each year, so there will be a potential sudden change in deaths_pc at the 
# start of each year. 
mortality <- inner_join(population_ar, mortality, by=c("year", "jurisdiction")) %>%
  select(-year) %>%
  as_tsibble(index=week, key = c(jurisdiction, cause)) %>%
  mutate(deaths_pc = 100000*(deaths / population))

mortality

mortality %>%
  filter(jurisdiction == "United States") %>%
  autoplot(deaths_pc)


### Decomposition ###

# Pre-2020
mortality %>% 
  filter(jurisdiction == "United States" & cause == "all" & week < yearweek("2020")) %>%
  model(stl = STL(deaths_pc)) %>% 
  components() %>%
  autoplot()


### Features ###

mortality %>% 
  features(deaths_pc, quantile) %>%
  arrange(`50%`)

# STL features
stl_feat <- mortality %>% 
  features(deaths_pc, feat_stl)
view(stl_feat)

# Features PCA
mortality_features <- mortality %>%
  features(deaths_pc, feature_set(pkgs = "feasts"))
pcs <- mortality_features %>%
  drop_na() %>%
  select(-jurisdiction, -cause) %>%
  prcomp(scale = FALSE) %>%
  augment(mortality_features %>% drop_na())
pcs %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, col = cause)) +
  geom_point() +
  theme(aspect.ratio = 1)


### Forecasting introduction ###

simple_fit <- mortality %>% 
  filter(week < yearweek("2017")) %>%
  model(
    mean = MEAN(deaths_pc),
    naive = NAIVE(deaths_pc),
    snaive = SNAIVE(deaths_pc)
  ) 

simple_fit

simple_fit %>%
  filter(jurisdiction == "United States" & cause == "all") %>%
  forecast(h = "5 years") %>%
  autoplot() +
  autolayer(mortality %>% 
              filter(jurisdiction == "United States" & cause == "all"), deaths_pc)

simple_fit %>%
  select(naive) %>%
  filter(jurisdiction == "United States" & cause == "all") %>%
  augment()

# simple_fit %>%
#   filter(jurisdiction == "United States" & cause == "all") %>%
#   augment() %>%
#   autoplot(.fitted) +
#   autolayer(mortality %>% 
#               filter(jurisdiction == "United States" & cause == "all" 
#                      & week < yearweek("2017")), deaths_pc)


### Residuals ###

simple_fit %>%
  filter(jurisdiction == "United States" & cause == "all") %>%
  select(naive) %>%
  gg_tsresiduals()

simple_fit %>%
  filter(jurisdiction == "United States" & cause == "all") %>%
  select(naive) %>%
  augment() %>%
  features(.innov, ljung_box ,lag = 156 / 5, dof = 0) # lag = T / 5


### Prediction intervals ###

simple_sim <- simple_fit %>%
  filter(jurisdiction == "United States" & cause == "all") %>%
  select(naive) %>%
  generate(h = 30, times = 10, bootstrap = TRUE)

mortality %>%
  filter(jurisdiction == "United States" & cause == "all") %>%
  ggplot(aes(x = week)) +
  geom_line(aes(y = deaths_pc)) +
  geom_line(aes(y = .sim, colour = as.factor(.rep)), data = simple_sim) + 
  guides(col = FALSE)

simple_fit %>%
  filter(jurisdiction == "United States" & cause == "all") %>%
  select(naive) %>%
  forecast(h = 30, bootstrap = TRUE) %>%
  autoplot() +
  autolayer(mortality %>% filter(jurisdiction == "United States" & cause == "all"), deaths_pc)


### Forecasting with decomposition ###

mortality_all_us <- mortality %>% 
  filter(jurisdiction == "United States" & cause == "all")
mortality_all_us_train <- mortality_all_us %>% 
  filter(week < yearweek("2017"))

dcmp <- mortality_all_us_train %>%
  model(STL(deaths_pc ~ trend(), robust = TRUE)) %>% 
  components() %>%
  select(-.model)

mortality_all_us_train %>% 
  model(STL(deaths_pc ~ trend(), robust = TRUE)) %>% 
  components() %>%
  autoplot()

fit_dcmp <- mortality_all_us_train %>% 
  model(stlf = decomposition_model(
    STL(deaths_pc ~ trend(), robust = TRUE),
    NAIVE(season_adjust)
  )) 
  
fit_dcmp %>%
  forecast(h = "5 years") %>%
  autoplot(mortality_all_us)
  
fit_dcmp %>% gg_tsresiduals()

fit_dcmp %>%
  filter(jurisdiction == "United States" & cause == "all") %>%
  augment() %>%
  features(.innov, ljung_box ,lag = 156 / 5, dof = 0) # Not sure on dof
  

### Evaluation ###

decomp_fit <- mortality %>%
  filter(week < yearweek("2017")) %>%
  model(stlf = decomposition_model(
    STL(deaths_pc ~ trend(), robust = TRUE),
    NAIVE(season_adjust)
  ))

decomp_fc <- decomp_fit %>%
  forecast(h = "5 years")

fabletools::accuracy(decomp_fc, mortality) %>%
  drop_na() %>%
  arrange(MASE)

decomp_fc %>%
  filter(jurisdiction == "South Dakota" & cause == "cancer") %>%
  autoplot(mortality %>% filter(jurisdiction == "South Dakota" & cause == "cancer"))


### Cross validation ###

fc_cv <- mortality_all_us %>%
  stretch_tsibble(.init = 3, .step = 1) %>%
  model(naive = NAIVE(deaths_pc)) %>%
  forecast(h = 4) %>%
  group_by(.id, .model) %>%
  mutate(h = row_number()) %>%
  ungroup()

fc_cv %>%
  filter(h == 4) %>%
  ggplot() +
  geom_line(data = mortality_all_us, aes(x=week, y=deaths_pc), color='grey') + 
  geom_line(aes(x=week, y=.mean), color='blue')


### Stationary models ###

mortality_all_us_diffs <- mortality_all_us %>%
  filter(week < yearweek("2020")) %>%
  mutate(
    diff = difference(deaths_pc), 
    diff2 = difference(diff),
    sdiff = difference(deaths_pc, 52),
    sdiffdiff = difference(sdiff)
  )

mortality_all_us_diffs %>% 
  autoplot(sdiffdiff)

# Different tests, different results. 
mortality_all_us_diffs %>% 
  mutate(diff0 = deaths_pc) %>%
  pivot_longer(matches("diff"), names_to = "diff_type", values_to = "diff") %>%
  group_by(diff_type) %>%
  features(diff, c(unitroot_nsdiffs, unitroot_kpss))


### Autoregressive models ###

mortality_all_us_train %>% 
  model(
    ar = AR(deaths_pc)
  ) %>%
  forecast(h = "5 years") %>%
  autoplot(mortality_all_us)

fit_ar <- mortality_all_us_train %>% 
  model(
    ar = AR(deaths_pc)
  )

fit_ar %>%
  augment() %>%
  autoplot(.fitted) +
  geom_line(data = mortality_all_us_train, aes(x=week, y=deaths_pc), color='grey')

fit_ar %>%
  forecast(h = "5 years") %>%
  autoplot(mortality_all_us)


### ARIMA ###

arima_fit <- mortality_all_us_train %>%
  model(
    arima210 = ARIMA(deaths_pc ~ pdq(2,1,0) + PDQ(0,0,0)),
    arima_stepwise = ARIMA(deaths_pc),
    arima_search = ARIMA(deaths_pc, stepwise=FALSE)
  )

glance(arima_fit) %>% arrange(AICc)

arima_fit$arima_search[[1]]$fit$spec

arima_fit %>% 
  select(arima_search) %>%
  gg_tsresiduals()

arima_fit %>%
  filter(jurisdiction == "United States" & cause == "all") %>%
  augment() %>%
  features(.innov, ljung_box ,lag = 156 / 5, dof = 6)

arima_fit %>%
  select(arima_search) %>%
  forecast(h = "5 years") %>% 
  autoplot() +
  autolayer(mortality_all_us, deaths_pc, linetype = "dashed") + 
  autolayer(mortality_all_us_train, deaths_pc)

