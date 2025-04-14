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

stand_max <- function(x) return(x/max(x))

# scaled features
Xstand <- dfOutcome %>%
  filter(date < traintest_date) %>%
  select(hosp, Positive, Tested) %>%
  mutate_all(.funs = stand_max) %>%
  as.matrix() %>%
  as.array()
# feed them to the reservoir
reservoir_state_stand <- predict_seq(node = reservoir, X = Xstand, reset = TRUE)
# plot the output
plot(reservoir_state_stand)

readout <- reservoirnet::createNode(nodeType = "Ridge", ridge = 0.1)
model <- reservoirnet::link(reservoir, readout)

# train set
yTrain <- dfOutcome %>% filter(outcomeDate <= traintest_date) %>% select(outcome)
xTrain <- dfOutcome %>% filter(outcomeDate <= traintest_date) %>% select(hosp, Positive, Tested)
# test set
xTest <- dfOutcome %>% select(hosp, Positive, Tested)

# standardize based on training set values
ls_fct_stand <- apply(xTrain,
                      MARGIN = 2,
                      FUN = function(x) function(feature) return(feature/(max(x))))

xTrainstand <- xTrain
xTeststand <- xTest
lapply(X = names(ls_fct_stand),
       FUN = function(x){
         xTrainstand[,x] <<- ls_fct_stand[[x]](feature = xTrain[,x])
         xTeststand[,x] <<- ls_fct_stand[[x]](feature = xTest[,x])
         return()
       })

# convert to array
lsdf <- lapply(list(yTrain = yTrain,
                    xTrain = xTrainstand,
                    xTest = xTeststand),
               function(x) as.array(as.matrix(x)))

### train the reservoir ridge output
fit <- reservoirnet::reservoirR_fit(node = model, X = lsdf$xTrain, Y = lsdf$yTrain, warmup = 30, reset = TRUE)

### predict with the reservoir
vec_pred <- reservoirnet::predict_seq(node = fit$fit, X = lsdf$xTest, reset = TRUE)

dfOutcome %>%
  mutate(pred = vec_pred) %>%
  ggplot(mapping = aes(x = outcomeDate)) +
  geom_line(mapping = aes(y = outcome, color = "observed")) +
  geom_line(mapping = aes(y = pred, color = "forecast")) +
  geom_vline(mapping = aes(color = "train-test sets", xintercept = traintest_date)) +
  scale_color_manual(values = c("#3772ff", "#080708", "#df2935")) +
  theme_bw() +
  labs(color = "", x = "Date", y = "Hospitalisations")