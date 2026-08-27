#!/usr/bin/Rscript
# Jordi Sevilla Fortuny

library(beastio)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)

# First argument may be the output CSV path; remaining arguments are log files.
# Usage: Rscript gather_results.R output.csv path/to/*.log
if (length(args) == 0) stop("Usage: Rscript gather_results.R output.csv file1.log [file2.log ...]")
if (endsWith(args[1], ".csv")) {
  outfile   <- args[1]
  log_files <- args[-1]
} else {
  outfile   <- "Simulations_results.csv"
  log_files <- args
}
if (length(log_files) == 0) stop("No log files provided.")


get_params <- function(file){
  # get attributes
  clean_file <- basename(file)
  prop <- strsplit(clean_file, "_")[[1]][2] %>% as.numeric()
  m <- strsplit(clean_file, "_")[[1]][3] %>% as.numeric()
  tree <- strsplit(clean_file, "_")[[1]][4]
  has_tiptypes <- !str_detect(clean_file, "tiptypes")

  # read log file
  data <- readLog(file)
  hpds <- beastio::getHPDMedian(data)

  interesting_params <- c("Re", "ReMultiplier", "SSFrac")
  rownames <- rownames(hpds)
  hpds %>% as_tibble() %>%
    mutate(param = rownames) %>%
    filter(
      param %in% interesting_params
    ) %>%
    mutate(tree = tree,
           tiptypes = ifelse(has_tiptypes, "yes", "no"),
           prop = prop,
           multiplier = m)
}

# Process data
lapply(
  log_files,
  get_params
  ) %>%
  bind_rows() %>%
  write_csv(outfile)
