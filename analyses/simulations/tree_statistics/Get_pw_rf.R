#!/usr/bin/env Rscript
# Jordi Sevilla Fortuny
# Script to get pairwise RFdists between trees

library(tidyverse)
library(ape)
library(phangorn)
library(treespace)

# Read trees
trees_0.05_2 <- read.tree("../neisseria_ss_data/sim_0.05_2.trees")
trees_0.05_10 <- read.tree("../neisseria_ss_data/sim_0.05_10.trees")
trees_0.25_2 <- read.tree("../neisseria_ss_data/sim_0.25_2.trees")
trees_0.25_10 <- read.tree("../neisseria_ss_data/sim_0.25_10.trees")

results <- data.frame(
  frac_1 = numeric(),
  effect_1 = numeric(),
  frac_2 = numeric(),
  effect_2 = numeric(),
  tree1_label = character(),
  tree2_label = character(),
  dist = numeric()
)


# get within group distances
get_distances <- function(tree_set, f, m){
  
  tmp <- data.frame(
    frac_1 = numeric(),
    effect_1 = numeric(),
    frac_2 = numeric(),
    effect_2 = numeric(),
    tree1_label = character(),
    tree2_label = character(),
    dist = numeric()
  )
  
  n_trees <- length(tree_set)
  
  for (tree1 in 1:n_trees) {
    for (tree2 in tree1:n_trees) {
      rf_dist = treespace::treeDist(tree_set[[tree1]], tree_set[[tree2]])
      
      tmp <- tmp %>%
        add_row(frac_1 = f, effect_1 = m, frac_2 = f, effect_2 = m, tree1_label = sprintf("tree_%s", tree1),
                tree2_label = sprintf("tree_%s", tree2), dist = rf_dist
                )
      
    }
  }
  
  tmp
}


df_0.05_10 <- get_distances(trees_0.05_10, 0.05, 10)
df_0.05_2 <- get_distances(trees_0.05_2, 0.05, 2)
df_0.25_10 <- get_distances(trees_0.25_10, 0.25, 10)
df_0.25_2 <- get_distances(trees_0.25_2, 0.25, 2)

results <- bind_rows(
  df_0.05_10,
  df_0.05_2,
  df_0.25_10,
  df_0.25_2
)


# Between groups dists
get_distances_bg <- function(tree_set_1, tree_set_2, f1, f2, m1, m2){
  
  tmp <- data.frame(
    frac_1 = numeric(),
    effect_1 = numeric(),
    frac_2 = numeric(),
    effect_2 = numeric(),
    tree1_label = character(),
    tree2_label = character(),
    dist = numeric()
  )
  
  n_trees_1 <- length(tree_set_1)
  n_trees_2 <- length(tree_set_2)
  
  for (tree1 in 1:n_trees_1) {
    for (tree2 in 1:n_trees_2) {
      rf_dist = treespace::treeDist(tree_set_1[[tree1]], tree_set_2[[tree2]])
      
      tmp <- tmp %>%
        add_row(frac_1 = f1, effect_1 = m1, frac_2 = f2, effect_2 = m2, tree1_label = sprintf("tree_%s", tree1),
                tree2_label = sprintf("tree_%s", tree2), dist = rf_dist
        )
      
    }
  }
  
  tmp
}


# get bg distances
a_b <- get_distances_bg(trees_0.05_10, trees_0.05_2, 0.05, 0.05, 10, 2)
a_c <- get_distances_bg(trees_0.05_10, trees_0.25_10, 0.05, 0.25, 10, 10)
a_d <- get_distances_bg(trees_0.05_10, trees_0.25_2, 0.05, 0.25, 10, 2)
b_c <- get_distances_bg(trees_0.05_2, trees_0.25_10, 0.05, 0.25, 2, 10)
b_d <- get_distances_bg(trees_0.05_2, trees_0.25_2, 0.05, 0.25, 2, 2)
c_d <- get_distances_bg(trees_0.25_10, trees_0.25_2, 0.25, 0.25, 10, 2)

results <- bind_rows(
  results,
  a_b, a_c, a_d,
  b_c, b_d,
  c_d
)

write_csv(results, file = "../neisseria_ss_data/pairwise_rf.csv")
