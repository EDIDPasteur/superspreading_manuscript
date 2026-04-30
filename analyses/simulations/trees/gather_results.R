#!/usr/bin/Rscript
# Jordi Sevilla Fortuny

library(beastio)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)


get_params <- function(file){
  # get attributes
  clean_file <- basename(file)
  prop <- strsplit(clean_file, "_")[[1]][2] %>% as.numeric()
  m <- strsplit(clean_file, "_")[[1]][3] %>% as.numeric()
  tree <- strsplit(clean_file, "_")[[1]][4]
  tiptypes <- !str_detect(clean_file, "tiptypes")
  
  # read log file
  data <- readLog(file)
  hpds <- beastio::getHPDMedian(data)
  
  interesting_params <- c("Re", "ReMultiplier", "SSFrac")
  rownames <- rownames(hpds)
  hpds %>% as.tibble() %>%
    mutate(param = rownames) %>%
    filter(
      param %in% interesting_params
    ) %>%
    mutate(tree = tree,
           tiptypes = ifelse(tiptypes, "yes", "no"),
           prop = prop,
           multiplier = m)
}

# Process data
lapply(
  args, 
  get_params
  ) %>%
  bind_rows() %>%
  write_csv("Simulations_results.csv")
