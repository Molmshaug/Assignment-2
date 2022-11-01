library(tidyverse)
library(exscidata)
library(gt)
data("cyclingstudy")

anaerob <- cyclingstudy %>% 
  select(subject, timepoint, lac.125:lac.375) %>% 
  filter(subject == 10,
         timepoint == "pre") %>% 
  pivot_longer(names_to = "watt",
               values_to = "laktat",
               cols = lac.125:lac.375,
               names_prefix = "lac.",
               names_transform = list(watt = as.numeric)) %>%
  filter(!is.na(laktat))

modanaerob <- lm(laktat ~ poly(watt, 3, raw = TRUE), data = anaerob)


anaerob %>% 
  ggplot(aes(watt, laktat, group = subject)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ poly(x, 3), color = "#4daf4a")

wattdata <- data.frame(watt = seq(from = 125, to = 375, by = 0.1))

wattdata$pred <- predict(modanaerob, newdata = wattdata)

laktat.pre <- wattdata %>% 
  filter(abs(pred - 2) == min(abs(pred - 2)) |
        abs(pred - 4) == min(abs(pred - 4))) %>% 
  mutate(timepoint = "pre",
         laktat = "2mmol")

laktat.pre[2, 4] <- "4mmol"

########################################################3

library(tidyverse)
library(exscidata)
data("cyclingstudy")

anaerob.meso1 <- cyclingstudy %>% 
  select(subject, timepoint, lac.125:lac.375) %>% 
  filter(subject == 10,
         timepoint == "meso1") %>% 
  pivot_longer(names_to = "watt",
               values_to = "laktat",
               cols = lac.125:lac.375,
               names_prefix = "lac.",
               names_transform = list(watt = as.numeric)) %>%
  filter(!is.na(laktat))

modanaerob.meso1 <- lm(laktat ~ poly(watt, 3, raw = TRUE), data = anaerob.meso1)


anaerob.meso1 %>% 
  ggplot(aes(watt, laktat, group = subject)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ poly(x, 3), color = "#4daf4a")

wattdata.meso1 <- data.frame(watt = seq(from = 125, to = 375, by = 0.1))

wattdata.meso1$pred <- predict(modanaerob.meso1, newdata = wattdata.meso1)

laktat.meso1 <- wattdata.meso1 %>% 
  filter(abs(pred - 2) == min(abs(pred - 2)) |
           abs(pred - 4) == min(abs(pred - 4))) %>% 
  mutate(timepoint = "meso1",
         laktat = "2mmol")

laktat.meso1[2, 4] <- "4mmol"

###############################################################

library(tidyverse)
library(exscidata)
data("cyclingstudy")

anaerob.meso2 <- cyclingstudy %>% 
  select(subject, timepoint, lac.125:lac.375) %>% 
  filter(subject == 10,
         timepoint == "meso2") %>% 
  pivot_longer(names_to = "watt",
               values_to = "laktat",
               cols = lac.125:lac.375,
               names_prefix = "lac.",
               names_transform = list(watt = as.numeric)) %>%
  filter(!is.na(laktat))

modanaerob.meso2 <- lm(laktat ~ poly(watt, 3, raw = TRUE), data = anaerob.meso2)


anaerob.meso2 %>% 
  ggplot(aes(watt, laktat, group = subject)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ poly(x, 3), color = "#4daf4a")

wattdata.meso2 <- data.frame(watt = seq(from = 125, to = 375, by = 0.1))

wattdata.meso2$pred <- predict(modanaerob.meso2, newdata = wattdata.meso2)

laktat.meso2 <- wattdata.meso2 %>% 
  filter(abs(pred - 2) == min(abs(pred - 2)) |
           abs(pred - 4) == min(abs(pred - 4))) %>% 
  mutate(timepoint = "meso2",
         laktat = "2mmol")

laktat.meso2[2, 4] <- "4mmol"

########################################################################

library(tidyverse)
library(exscidata)
data("cyclingstudy")

anaerob.meso3 <- cyclingstudy %>% 
  select(subject, timepoint, lac.125:lac.375) %>% 
  filter(subject == 10,
         timepoint == "meso3") %>% 
  pivot_longer(names_to = "watt",
               values_to = "laktat",
               cols = lac.125:lac.375,
               names_prefix = "lac.",
               names_transform = list(watt = as.numeric)) %>%
  filter(!is.na(laktat))

modanaerob.meso3 <- lm(laktat ~ poly(watt, 3, raw = TRUE), data = anaerob.meso3)


anaerob.meso3 %>% 
  ggplot(aes(watt, laktat, group = subject)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ poly(x, 3), color = "#4daf4a")

wattdata.meso3 <- data.frame(watt = seq(from = 125, to = 375, by = 0.1))

wattdata.meso3$pred <- predict(modanaerob.meso3, newdata = wattdata.meso3)

laktat.meso3 <- wattdata.meso3 %>% 
  filter(abs(pred - 2) == min(abs(pred - 2)) |
           abs(pred - 4) == min(abs(pred - 4))) %>% 
  mutate(timepoint = "meso3",
         laktat = "2mmol")

  laktat.meso3[2, 4] <- "4mmol"
  
         

############################################################################


full_data<- laktat.pre %>% 
  full_join(laktat.meso1) %>% 
  full_join(laktat.meso2) %>% 
  full_join(laktat.meso3)
 
    
  full_data_rel <- full_data %>%
     select(-pred) %>% 
     pivot_wider(names_from = timepoint,
                values_from = watt) %>% 
    mutate(diff1 = pre - meso1,
           diff2 = meso1 - meso2,
           diff3 = meso2 - meso3) %>% 
    group_by(laktat) %>% 
    summarise(s = sd(c(diff1, diff2, diff3)),
              m = mean(c(pre, meso1, meso2, meso3)),
              te = s / sqrt(2),
              cv = 100 * (te / m))
  
  gt(full_data_rel)
    

  
  