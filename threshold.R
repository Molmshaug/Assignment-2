library(tidyverse)
library(exscidata)
library(modelr)
library(gt)
library(gtExtras)

# import data
data(cyclingstudy)

# new function containing a model-function for later iteration
subject_model <- function(df) { 
  lm(lactate ~ poly(x, 3, raw = TRUE), data = df)
}

# new function for later iteration
df_fun <- function(x) {
  as.data.frame(x)
}



# Table showing LT
cyclingstudy %>%
  select(subject, group, timepoint, lac.125:lac.375) %>%
  pivot_longer(names_to = "x", 
               values_to = "lactate", 
               names_prefix = "lac.",
               names_transform = list(x = as.numeric),
               cols = lac.125:lac.375) %>% 
  filter(!is.na(lactate)) %>% 
  group_by(subject, timepoint) %>%
  nest() %>% 
  mutate(model = map(data, subject_model),
         watt_range = list(watt = seq(from = 125, to = 375, by = 0.1)),
         df = map(watt_range, df_fun),
         pred = map2(df, model, add_predictions)) %>%
  unnest(pred) %>% 
  filter(abs(pred - 4) == min(abs(pred - 4)) | 
           abs(pred - 2) == min(abs(pred - 2))) %>% 
  select(subject, timepoint, x, pred) %>% 
  mutate(lactate = if_else(pred > 3, "ana", "aer")) %>% 
  select(-pred) %>% 
  rename(watt = x) %>% 
  pivot_wider(names_from = c(timepoint, lactate),
              values_from = watt) %>%
  select(subject, pre_aer, meso1_aer, meso2_aer, meso3_aer,
         pre_ana, meso1_ana, meso2_ana, meso3_ana) %>% 
  filter(subject <= 4) %>% 
  ungroup() %>% 
  gt() %>% 
  tab_spanner(columns = pre_aer:meso3_aer, label = "LT at 2 mmol*L") %>% 
  tab_spanner(columns = pre_ana:meso3_ana, label = "LT at 4 mmol*L") %>% 
  fmt_number(columns = pre_aer:meso3_ana, decimals = 0) %>% 
  gt_add_divider(columns = meso3_aer, color = "lightgrey", style = "dotted") %>%
  cols_label(subject = "Subject",
             pre_aer = "Pre",
             meso1_aer = "Meso1",
             meso2_aer = "Meso2",
             meso3_aer = "Meso3",
             pre_ana = "Pre", 
             meso1_ana = "Meso1",
             meso2_ana = "Meso2",
             meso3_ana = "Meso3") %>% 
  tab_header(title = md("**Calculated Lactate Threshold**"),
             subtitle = "Using fixed blood lactate concentration method (FBLC)") %>% 
  tab_footnote(footnote = "Showing the first 4 rows in the dataset") %>% 
  tab_footnote(footnote = "Abbriviations: LT, lactate threshold")
  
  


# Table showing Reliability
cyclingstudy %>%
  select(subject, group, timepoint, lac.125:lac.375) %>%
  pivot_longer(names_to = "x", 
               values_to = "lactate", 
               names_prefix = "lac.",
               names_transform = list(x = as.numeric),
               cols = lac.125:lac.375) %>% 
  filter(!is.na(lactate)) %>% 
  group_by(subject, timepoint) %>%
  nest() %>% 
  mutate(model = map(data, subject_model),
         watt_range = list(watt = seq(from = 125, to = 375, by = 0.1)),
         df = map(watt_range, df_fun),
         pred = map2(df, model, add_predictions)) %>%
  unnest(pred) %>% 
  filter(abs(pred - 4) == min(abs(pred - 4)) | 
           abs(pred - 2) == min(abs(pred - 2))) %>% 
  select(subject, timepoint, x, pred) %>% 
  mutate(lactate = if_else(pred > 3, "4 mmol*L", "2 mmol*L")) %>% 
  select(-pred) %>% 
  rename(watt = x) %>% 
  pivot_wider(names_from = timepoint,
              values_from = watt) %>% 
  mutate(diff1 = meso1 - pre,
         diff2 = meso2 - meso1,
         diff3 = meso3 - meso2) %>% 
  pivot_longer(names_to = "diff",
               values_to = "value",
               cols = diff1:diff3) %>% 
  group_by(lactate) %>% 
  summarise(m = mean(c(pre, meso1, meso2, meso3), na.rm = T),
            s = round(sd(value, na.rm = T), 1),
            te = round(s / sqrt(2), 1), 
            cv = round(100 * (te / m), 1)) %>% 
  select(-s, -te) %>% 
  gt(rowname_col = "lactate") %>% 
  fmt_number(columns = m:cv, decimals = 1) %>% 
  cols_label(m = "M",
             cv = "CV") %>%
  tab_stubhead(label = "Lactate [BLa-]") %>% 
  tab_header(title = md("**Reliability of lactate threshold**"),
             subtitle = "Using fixed blood lactate concentration method (FBLC)") %>% 
  tab_footnote(footnote = "Abbriviations: M, mean; CV, coefficient of variations")
  
  
  
  





  

