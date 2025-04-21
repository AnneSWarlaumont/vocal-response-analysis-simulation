source('VocEventSim.R')

sim_length = 10*60*60
rthresh = 1

a2_minp = .0001
a2_respsensitivity = 1
a2_othersensitivity = 1.5

a1_minp = .000001
a1_respsensitivity = 1.25
a1_othersensitivity = 1

a1_meanlog = 0
a2_meanlog = 0
a1_sdlog = .2
a2_sdlog = .2
a1_maxp = .3
a2_maxp = .3
a1_p_voc = runif(1,min=a1_minp,max=a1_maxp)
a2_p_voc = runif(1,min=a2_minp,max=a2_maxp)

voc_records = two_agent_vocal_sim(sim_length,a1_p_voc,a2_p_voc,a1_meanlog,a2_meanlog,a1_sdlog,a2_sdlog,a1_minp,a2_minp,a1_maxp,a2_maxp,a1_othersensitivity,a2_othersensitivity,a1_respsensitivity,a2_respsensitivity,rthresh)
t = seq(from = 1, to = sim_length)
a1_voc_record = voc_records[[1]]
a2_voc_record = voc_records[[2]]
a2toa1_r_record = c(voc_records[[3]],NA)
a1toa2_r_record = c(voc_records[[4]],NA)
voc_dat = data.frame(t,a1_voc_record,a2_voc_record,a2toa1_r_record,a1toa2_r_record)

library(reservoirnet)
library(ggplot2)
library(dplyr)
library(slider)
library(ggplot2)

a1_dist_forecast = 60
traintest_sec = 5*60*60

voc_predict_dat <- voc_dat %>%
  # outcome at 14 seconds later
  mutate(a1outcome = lead(x = a1_voc_record, a1_dist_forecast),
         a1outcome_t = t + a1_dist_forecast) %>%
  #rolling average variables
  mutate_at(.vars = c("a1_voc_record","a2_voc_record","a2toa1_r_record","a1outcome"),
            .funs = function(x) slider::slide_dbl(.x = x,
                                        .before = 10*60,
                                        .f = sum,
                                        na.rm = TRUE))


voc_predict_dat %>%
  tidyr::pivot_longer(cols = c("a1_voc_record","a2_voc_record","a2toa1_r_record")) %>%
  ggplot2::ggplot(mapping = aes(x = t, y = value)) +
  geom_line() +
  facet_grid(name ~ ., scales = "free_y") +
  theme_bw() +
  geom_vline(mapping = aes(color = "train-test sets", xintercept = traintest_sec)) +
  labs(color = "") +
  theme(legend.position = "bottom")

reservoir <- reservoirnet::createNode(nodeType = "Reservoir",
                                      seed = 1,
                                      units = 500,
                                      lr = 0.7,
                                      sr = 1,
                                      input_scaling = 1)


stand_max <- function(x) return(x/max(x))

# select and scale explanatory and transform it to an array
Xstand <- voc_predict_dat %>%
  filter(t < traintest_sec) %>%
  select(a1_voc_record, a2_voc_record, a2toa1_r_record) %>%
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
yTrain <- voc_predict_dat %>% filter(t <= traintest_sec) %>% select(a1outcome)
xTrain <- voc_predict_dat %>% filter(t <= traintest_sec) %>% select(a1_voc_record, a2_voc_record, a2toa1_r_record)
# test set
xTest <- voc_predict_dat %>% select(a1_voc_record, a2_voc_record, a2toa1_r_record)

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

voc_predict_dat %>%
  mutate(pred = vec_pred) %>%
  ggplot(mapping = aes(x = a1outcome_t)) +
  geom_line(mapping = aes(y = a1outcome, color = "observed")) +
  geom_line(mapping = aes(y = pred, color = "forecast")) +
  geom_vline(mapping = aes(color = "train-test sets", xintercept = traintest_sec)) +
  scale_color_manual(values = c("#3772ff", "#080708", "#df2935")) +
  theme_bw() +
  labs(color = "", x = "Time", y = "a1 vocalizations")
