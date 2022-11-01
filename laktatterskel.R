library(tidyverse)
library(exscidata)
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

laktatterskel <- wattdata %>% 
  filter(abs(pred - 2) == min(abs(pred - 2)) |
        abs(pred - 4) == min(abs(pred - 4)))

        