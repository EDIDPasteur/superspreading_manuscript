#!/usr/bin/env Rscript
# Jordi Sevilla Fortuny
# Script to get pairwise coales index between trees

library(tidyverse)
library(ape)
library(phangorn)
library(phyloTop)

# Read trees
trees_0.05_2 <- read.tree("../neisseria_ss_data/sim_0.05_2.trees")
trees_0.05_10 <- read.tree("../neisseria_ss_data/sim_0.05_10.trees")
trees_0.25_2 <- read.tree("../neisseria_ss_data/sim_0.25_2.trees")
trees_0.25_10 <- read.tree("../neisseria_ss_data/sim_0.25_10.trees")

results <- data.frame(
  tree = character(),
  frac = numeric(),
  effect = numeric(), 
  index = numeric()
)

get_indexes <- function(tree_set, f, m){
  
  tmp <- data.frame(
    tree = character(),
    frac = numeric(),
    effect = numeric(), 
    index = numeric()
  )
  
  ntrees <- length(tree_set)
  
  for (x in 1:ntrees) {
    ci <- phyloTop::colless.phylo(tree_set[[x]], normalise = T)
    
    tmp <- tmp %>%
      add_row(
        tree = sprintf("tree_%s", x),
        frac = f,
        effect = m,
        index = ci
      )
  }
  
  tmp
}


results <- bind_rows(
  get_indexes(trees_0.05_10, 0.05, 10),
  get_indexes(trees_0.05_2, 0.05, 2),
  get_indexes(trees_0.25_10, 0.25, 10),
  get_indexes(trees_0.25_2, 0.25, 2)
)


write_csv(results, "../neisseria_ss_data/trees_ci.csv")
