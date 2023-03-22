library(tidyverse)
library(devtools)

input_samples <- read_delim("/Users/rochelleyap/Desktop/JgiAllSampleCounts.csv", delim = ",")
head(input_samples)


input_samples.sub <- input_samples %>% 
  dplyr::select(-Chr, -Start, -End, -Strand, -Length)

head(input_samples.sub)
