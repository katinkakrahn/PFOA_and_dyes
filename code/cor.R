# necessary packages can be installed as follows:
# install.packages(c("tidyverse", "skimr", "ggcorrplot", "broom"))
library(tidyverse)
library(readxl)
library(skimr)
library(ggcorrplot)
library(broom)

# improting the data and removing standard deviation columns
data <- read_xlsx("data_raw/20230317_BC additional data_KAS.xlsx",
                  n_max = 23, # this specifies the amount of lines to read int
                  na = c("-", "n.a.", "", "n.d.", "Yaxin")) |> # this specifies what should count as missing value
  select(-contains("stdev"))

skim(data)

# scaling of all numeric variables. NB: log kF should _not_ be scaled!
data_scale <- data |>  
  mutate(across(where(is.numeric) & !matches("Y"), ~ as.vector(scale(.x)))) # the Y variable name needs to be replaced with the correct name of the log kF column

skim(data_scale)

# (Pearson) correlation analysis
data_cor <- data_scale |>
  select(-c(New_sample_id, Old_sample_id, C_arom_BPCA)) |> # remove any non-numeric variables. NB: Kf should also be removed here!
  cor(use = "pairwise.complete.obs", # to include the maximum of data even when values are missing
      method = "pearson") # this is the default

data_cor

# visualization of the correlation analysis
ggcorrplot(data_cor, hc.order = TRUE, type = "lower", lab = TRUE)
