# This is my following along the reservoirnet tutorial posted at https://rdrr.io/cran/reservoirnet/f/vignettes/basic_usage_01.Rmd

library(reservoirnet)
library(ggplot2)
library(dplyr)

data("dfCovid")
dist_forecast = 14
traintest_date = as.Date("2022-01-01")

dfOutcome <- dfCovid %>%
  # outcome at 14 days
  mutate(outcome = lead(x = hosp, n = dist_forecast),
         outcomeDate = date + dist_forecast) %>%
  # rolling average for iptcc and positive_pcr
  mutate_at(.vars = c("Positive","Tested"),
            .funs = function(x) slider::slide_dbl(.x = x,
                                                  .before = 6,
                                                  .f = mean))

dfOutcome %>%
  tidyr::pivot_longer(cols = c("hosp", "Positive", "Tested")) %>%
  ggplot2::ggplot(mapping = aes(x = date, y = value)) + 
  geom_line() +
  facet_grid(name ~ ., scales = "free_y") +
  theme_bw() +
  geom_vline(mapping = aes(color = "train-test sets", xintercept = traintest_date)) +
  labs(color = "") +
  theme(legend.position = "bottom")

reservoir <- reservoirnet::createNode(nodeType = "Reservoir",
                                      seed = 1,
                                      units = 500,
                                      lr = 0.7,
                                      sr = 1,
                                      input_scaling = 1)

## select explanatory and transform it to an array
X <- dfOutcome %>%
  filter(outcomeDate < traintest_date) %>%
  select(hosp, Positive, Tested) %>%
  as.matrix() %>%
  as.array()

reservoir_state <- predict_seq(node = reservoir, X = X, reset = TRUE)

plot(reservoir_state)
