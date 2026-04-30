# Jordi Sevilla Fortuny - 20260430

# Script to generate the figures for the manuscript:
#   -Quantifying superspreading in bacterial STI outbreaks using phylodynamics

# libraries
library(tidyverse)
library("latex2exp")
library(ggpubr)
library(beastio)
library(coda)

# Set wd
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Aesthetics:
cluster_order <- c(
  "180", "146", "557", "502", "230", "26", "818", "879", "898", "980",
  "33", "1003", "1094", "15", "1067", "844", "1177", "560"
)

data_types_colors <- c(
  seqs = "#5EBCB2",
  trees = "#AB5A7D"
)

data_types_names <- c(
  seqs = "Sequences",
  trees = "Trees"
)

conditions_names <- c(
  "f=0.05; M=10" = "A",
  "f=0.05; M=2" = "B",
  "f=0.25; M=10" = "C",
  "f=0.25; M=2" = "D"
)

skyline_colors <- c(
  "1" = "#289BFF",
  "2" = "#F28E2B"
)

skyline_names <- c(
  "1" = "Before lockdown",
  "2" = "After lockdown"
)

# FIGURE 1 ####

# simulated trees
data_trees <- read_csv(
  "data/Simulation_results_trees.csv"
)

## HPDs
accuracy_trees <- data_trees %>%
  filter(tiptypes == "no") %>%
  mutate(group = paste("f=", prop, "; ", "M=", multiplier, sep = "")) %>%
  pivot_wider(names_from = param, values_from = c(med, upper, lower)) %>%
  mutate(data = "trees")

# simulated seqs
data_seqs <- read_csv(
  "data/Simulation_results_seqs.csv"
)

## HPDs
accuracy_seqs <- data_seqs %>%
  filter(tiptypes == "no") %>%
  mutate(group = paste("f=", prop, "; ", "M=", multiplier, sep = "")) %>%
  pivot_wider(names_from = param, values_from = c(med, upper, lower)) %>%
  mutate(data = "seqs")


# subplots
## comparisions

comparisions <- bind_rows(accuracy_seqs, accuracy_trees) %>%
  ggplot() +
  aes(x = med_SSFrac, y = med_ReMultiplier, color = data) +
  geom_segment(aes(y = lower_ReMultiplier, yend = upper_ReMultiplier), alpha = 0.3) +
  geom_segment(aes(x = lower_SSFrac, xend = upper_SSFrac), alpha = 0.3) +
  geom_point(alpha = 0.9, size = 2.5) +
  geom_hline(aes(yintercept = multiplier), size = 2, color = "darkgreen", alpha = 0.5) +
  geom_vline(aes(xintercept = prop), size = 2, color = "darkgreen", alpha = 0.5) +
  scale_color_manual(
    values = data_types_colors,
    labels = data_types_names, guide = "none"
  ) +
  facet_grid(~group,
    labeller = labeller(group = conditions_names)
  ) +
  theme_bw() +
  labs(
    x = "Fraction of superspreaders (95% HPD)",
    y = "Relative effect of superspreading (95% HPD)",
    color = "Type of data"
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.title = element_text(size = 9),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_line(colour = "gray80", linewidth = 1),
    legend.position = "bottom"
  )

## stats multiplier
multiplier <- data_seqs %>%
  filter(param == "ReMultiplier", tiptypes == "no") %>%
  mutate(group = paste("f=", prop, "; ", "M=", multiplier, sep = "")) %>%
  group_by(group, tiptypes) %>%
  summarise(
    accuracy = sum(multiplier >= lower & multiplier <= upper) / n(),
    error = abs(mean(med - multiplier)),
    precission = mean((upper - lower) / med)
  ) %>%
  pivot_longer(c(accuracy, error, precission),
    names_to = "stat",
    values_to = "value"
  ) %>%
  mutate(data = "seqs") %>%
  bind_rows(
    data_trees %>%
      filter(param == "ReMultiplier", tiptypes == "no") %>%
      mutate(group = paste("f=", prop, "; ", "M=", multiplier, sep = "")) %>%
      group_by(group, tiptypes) %>%
      summarise(
        accuracy = sum(multiplier >= lower & multiplier <= upper) / n(),
        error = abs(mean(med - multiplier)),
        precission = mean((upper - lower) / med)
      ) %>%
      pivot_longer(c(accuracy, error, precission),
        names_to = "stat", values_to = "value"
      ) %>%
      mutate(data = "trees")
  ) %>%
  ggplot() +
  aes(x = group, fill = data, y = value) +
  geom_col(position = "dodge") +
  scale_fill_manual(
    values = data_types_colors,
    labels = data_types_names, guide = "none"
  ) +
  theme_bw() +
  facet_wrap(~stat, scales = "free_y", labeller = labeller(
    stat = c(
      "accuracy" = "Accuracy",
      "error" = "Error",
      "precission" = "Uncertainty"
    )
  )) +
  scale_x_discrete(labels = c("A", "B", "C", "D")) +
  labs(x = "Simulation setting", y = "Value", data = "Type of data")

## stats ssfrac
ssfrac <- data_seqs %>%
  filter(param == "SSFrac", tiptypes == "no") %>%
  mutate(group = paste("f=", prop, "; ", "M=", multiplier, sep = "")) %>%
  group_by(group, tiptypes) %>%
  summarise(
    accuracy = sum(prop >= lower & prop <= upper) / n(),
    error = abs(mean(med - prop)),
    precission = mean((upper - lower) / med)
  ) %>%
  pivot_longer(c(accuracy, error, precission),
    names_to = "stat",
    values_to = "value"
  ) %>%
  mutate(data = "seqs") %>%
  bind_rows(
    data_trees %>%
      filter(param == "SSFrac", tiptypes == "no") %>%
      mutate(group = paste("f=", prop, "; ", "M=", multiplier, sep = "")) %>%
      group_by(group, tiptypes) %>%
      summarise(
        accuracy = sum(prop >= lower & prop <= upper) / n(),
        error = abs(mean(med - prop)),
        precission = mean((upper - lower) / med)
      ) %>%
      pivot_longer(c(accuracy, error, precission),
        names_to = "stat",
        values_to = "value"
      ) %>%
      mutate(data = "trees")
  ) %>%
  ggplot() +
  aes(x = group, fill = data, y = value) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = data_types_colors, labels = data_types_names) +
  theme_bw() +
  facet_wrap(~stat, scales = "free_y", labeller = labeller(
    stat = c(
      "accuracy" = "Accuracy",
      "error" = "Error",
      "precission" = "Uncertainty"
    )
  )) +
  scale_x_discrete(labels = c("A", "B", "C", "D")) +
  labs(x = "Simulation setting", y = "Value", fill = "Type of data") +
  theme(legend.position = "bottom")


# MAIN FIGURE
ggarrange(
  comparisions,
  multiplier,
  ssfrac,
  nrow = 3,
  labels = c("A)", "B)", "C)"),
  heights = c(0.5, 0.25, 0.3)
)


# FIGURE 2 ####
data_no_ss <- read_csv(
  "data/Simulation_results_no_ss.csv"
) %>%
  mutate(tree = gsub(".log", "", tree))

# exponential prior
exponential_prior <- data.frame(sample = rexp(1000000, rate = 0.2)) %>%
  ggplot() +
  aes(y = sample, x = ..scaled..) +
  geom_density(fill = "#BC8357", alpha = 0.3, color = "#BC8357", linewidth = 2) +
  theme_bw() +
  ylim(c(NA, 20)) +
  scale_x_reverse() +
  theme(
    axis.text.y = element_blank(), axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(y = "Prior Exp(5)", x = "")

# Beta prior
beta_prior <- data.frame(sample = rbeta(1000000, 1, 1)) %>%
  ggplot() +
  aes(x = sample, y = ..scaled..) +
  geom_density(fill = "#BC8357", alpha = 0.3, color = "#BC8357", linewidth = 2) +
  theme_bw() +
  scale_y_reverse() +
  theme(
    axis.text.x = element_blank(), axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(x = "Prior Beta(1,1)", y = "")

# data figure
figure_presentation <- data_no_ss %>%
  filter(param %in% c("ReMultiplier", "SSFrac")) %>%
  pivot_wider(names_from = param, values_from = c("lower", "med", "upper")) %>%
  ggplot() +
  aes(x = med_SSFrac, y = med_ReMultiplier) +
  geom_errorbar(aes(ymin = lower_ReMultiplier, ymax = upper_ReMultiplier),
    alpha = 0.3, color = "#BC8357"
  ) +
  geom_errorbar(aes(xmax = upper_SSFrac, xmin = lower_SSFrac),
    alpha = 0.3, color = "#BC8357"
  ) +
  geom_point(size = 3, color = "#BC8357") +
  theme_bw() +
  geom_hline(aes(yintercept = 1), linewidth = 3, color = "#5CBEB1", alpha = 0.7) +
  geom_vline(aes(xintercept = 0), linewidth = 3, color = "#5CBEB1", alpha = 0.7) +
  labs(
    x = TeX("$f_{ss}$ (95% HPD)"),
    y = TeX("$T_a$ (95% HPD)")
  )

# MAIN FIGURE
ggarrange(exponential_prior,
  figure_presentation,
  NA,
  beta_prior,
  ncol = 2, nrow = 2, align = "hv",
  widths = c(0.3, 0.9), heights = c(0.9, 0.3)
)


# FIGURE 3 ####
# read the log file

log <- readLog("data/BDMMPrime_neisseria_superspreading.thinned2x.log.gz", as.mcmc = F)

stats <- lapply(names(log), function(param_name) {
  x <- log[[param_name]]

  m <- mean(x)
  hpd <- HPDinterval(as.mcmc(x), prob = 0.95)

  data.frame(
    param = param_name,
    med = m,
    lower = hpd[1],
    upper = hpd[2]
  )
}) %>%
  bind_rows()

# Origin figure
data_dates <- read_csv("data/tip_dates.txt")
data_dates_last <- read_csv("data/last_sample_date.txt")

all_data_origin <- stats %>%
  extract(param,
    into = c("nombre", "cluster_id"),
    regex = "([A-Za-z]+)(\\d+)",
    convert = TRUE
  ) %>%
  filter(nombre == "originBDMMPrime") %>%
  left_join(data_dates, by = c("cluster_id" = "cluster")) %>%
  left_join(data_dates_last) %>%
  summarise(
    cluster = cluster_id,
    end = date,
    origin = date - med,
    upper_o = date - lower,
    lower_o = date - upper,
    first = last_samp_date
  )

## figure
plot_origin <- ggplot(all_data_origin) +
  aes(y = reorder(factor(cluster), -origin), x = end) +
  geom_rect(aes(xmin = 2020.216, xmax = 2020.243, ymin = -Inf, ymax = Inf),
    fill = "darkgreen", alpha = 1
  ) +
  geom_point(size = 3) +
  geom_segment(aes(x = end, xend = origin)) +
  geom_point(aes(x = origin), size = 3) +
  geom_segment(aes(x = lower_o, xend = upper_o),
    lineend = "round",
    linewidth = 3, alpha = 0.5
  ) +
  geom_point(
    size = 2, color = "darkred", aes(x = first),
    shape = 4, stroke = 1.5
  ) +
  theme_bw() +
  labs(y = "Cluster", x = "Time")

# subplot clockrate

plot_cr <- stats %>%
  extract(param,
    into = c("nombre", "cluster_id"),
    regex = "([A-Za-z]+)(\\d+)",
    convert = TRUE
  ) %>%
  filter(nombre == "clockRate") %>%
  ggplot() +
  aes(x = factor(cluster_id, cluster_order), y = med) +
  geom_point(size = 3, colour = "#BC8357") +
  geom_segment(aes(y = lower, yend = upper),
    lineend = "round",
    linewidth = 3, alpha = 0.5, colour = "#BC8357"
  ) +
  theme_bw() +
  labs(x = "Cluster", y = "Clock rate (95% HPD)") +
  geom_density(
    inherit.aes = F, aes(y = estimate, x = ..scaled.. * 3),
    data = data.frame(
      sample = c(1:1000000),
      estimate = rlnorm(
        n = 1000000, meanlog = log$ClockM,
        sdlog = log$ClockS
      )
    ),
    linewidth = 2, color = "#BC8357",
    fill = "#BC8357", alpha = 0.5,
    position = position_nudge(x = 19)
  ) +
  geom_hline(aes(yintercept = exp(mean(log$ClockM) + ((mean(log$ClockS))**2 / 2))),
    linetype = 2, linewidth = 1, color = "#BC8357"
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(~., name = "Prior distribution"),
    limits = c(NA, 0.000022)
  ) +
  theme(
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )


# subplot sampling prop

plot_sp <- stats %>%
  extract(param,
    into = c("nombre", "cluster_id", "time_slice"),
    regex = "([A-Za-z]+)(\\d+)\\.(\\d+)",
    convert = TRUE
  ) %>%
  filter(nombre == "samplingProportionEpi", time_slice == 1) %>%
  ggplot() +
  aes(x = factor(cluster_id, cluster_order), y = med) +
  geom_point(size = 3, colour = "#BC8357") +
  geom_segment(
    aes(y = lower, yend = upper),
    lineend = "round", linewidth = 3,
    alpha = 0.5, colour = "#BC8357"
  ) +
  theme_bw() +
  labs(
    x = "Cluster",
    y = "Sampling proportion (95% HPD)"
  ) +
  geom_density(
    inherit.aes = F,
    trim = TRUE,
    aes(y = estimate, x = ..scaled.. * 3),
    data = data.frame(
      sample = c(1:1000000),
      estimate = rbeta(
        n = 1000000, shape1 = log$SamplingPropAlpha,
        shape2 = log$SamplingPropBeta
      )
    ),
    linewidth = 2, color = "#BC8357",
    fill = "#BC8357", alpha = 0.5,
    position = position_nudge(x = 19)
  ) +
  geom_hline(
    aes(
      yintercept = mean(log$SamplingPropAlpha) / (mean(log$SamplingPropAlpha) + mean(log$SamplingPropBeta))
    ),
    linetype = 2, linewidth = 1, color = "#BC8357"
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(~., name = "Prior distribution"),
    limits = c(NA, 0.4)
  ) +
  theme(
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )


## FINAL FIGURE
ggarrange(
  ggarrange(plot_origin),
  ggarrange(
    plot_sp,
    plot_cr,
    ncol = 2,
    align = "h",
    labels = c("B)", "C)"),
    label.x = -0.005
  ),
  ncol = 1,
  align = "v",
  labels = c("A)", ""),
  label.x = -0.005
)

# FIGURE 4 ####

# posterior of change after the lockdown
change_posterior <- log %>%
  select(c(contains("Re"), contains("SSFrac"))) %>%
  mutate(sample = seq_len(n())) %>%
  pivot_longer(!sample, names_to = "nombre", values_to = "value") %>%
  extract(nombre,
    into = c("nombre", "cluster_id", "time_slice"),
    regex = "([A-Za-z]+)(\\d+).(\\d+)",
    convert = TRUE
  ) %>%
  filter(!is.na(time_slice)) %>%
  group_by(sample, nombre, cluster_id) %>%
  summarise(gt = value[2] < value[1]) %>%
  ungroup() %>%
  group_by(nombre, cluster_id) %>%
  summarise(p = mean(gt)) %>%
  ungroup() %>%
  ggplot() +
  aes(x = p, y = nombre) +
  geom_boxplot(outlier.alpha = 0, color = "#CC79A7", size = 1, fill = "#CC79A7", alpha = 0.1, width = 0.2) +
  geom_jitter(size = 2, color = "#CC79A7", position = position_jitter(height = 0.1), size = 2) +
  theme_bw() +
  scale_y_discrete(labels = c(TeX("$R_e$"), TeX("$T_a$"), TeX("$f_{ss}$"))) +
  labs(x = "P(after lockdown < before lockdown)", y = "", title = "") +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.text.y = element_text(size = 15)
  ) +
  geom_vline(aes(xintercept = 0.5), linetype = 2, color = "#CC79A7", linewidth = 1.5, alpha = 0.6)

# posterior distributions for the ss params

stats_sep <- stats %>%
  extract(param,
    into = c("nombre", "cluster_id", "time_slice"),
    regex = "([A-Za-z]+)(\\d+)\\.(\\d+)",
    convert = TRUE
  ) %>%
  filter(!is.na(time_slice))

## Re
Re <- log %>%
  select(c(contains("Re"), contains("SSFrac"))) %>%
  mutate(sample = seq_len(n())) %>%
  pivot_longer(!sample, names_to = "nombre", values_to = "value") %>%
  extract(nombre,
    into = c("nombre", "cluster_id", "time_slice"),
    regex = "([A-Za-z]+)(\\d+).(\\d+)",
    convert = TRUE
  ) %>%
  mutate(cluster_id = factor(cluster_id)) %>%
  left_join(stats_sep %>% mutate(cluster_id = as.factor(cluster_id))) %>%
  filter((value > lower) & (value < upper)) %>%
  filter(nombre == "Re") %>%
  ggplot() +
  aes(x = factor(cluster_id, cluster_order), y = value, color = as.factor(time_slice), fill = as.factor(time_slice)) +
  geom_violin(scale = "width", width = 0.6, position = position_dodge(width = 0.6), alpha = 0.1) +
  geom_point(aes(y = med),
    data = stats_sep %>% filter(nombre != "samplingProportionEpi") %>% filter(nombre == "Re"),
    position = position_dodge(width = 0.6), size = 2.5
  ) +
  scale_color_manual("Time slice",
    labels = skyline_names,
    values = skyline_colors, guide = "none"
  ) +
  scale_fill_manual("Time slice",
    labels = skyline_names,
    values = skyline_colors, guide = "none"
  ) +
  theme_bw() +
  labs(y = TeX("Estimate and 95% HPD"), x = "Cluster", title = TeX("$R_e$")) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )


## Ta
ReM <- log %>%
  select(c(contains("Re"), contains("SSFrac"))) %>%
  mutate(sample = seq_len(n())) %>%
  pivot_longer(!sample, names_to = "nombre", values_to = "value") %>%
  extract(nombre,
    into = c("nombre", "cluster_id", "time_slice"),
    regex = "([A-Za-z]+)(\\d+).(\\d+)",
    convert = TRUE
  ) %>%
  mutate(cluster_id = factor(cluster_id)) %>%
  left_join(stats_sep %>% mutate(cluster_id = as.factor(cluster_id))) %>%
  filter((value > lower) & (value < upper)) %>%
  filter(nombre == "ReMultiplier") %>%
  ggplot() +
  aes(x = factor(cluster_id, cluster_order), y = value, color = as.factor(time_slice), fill = as.factor(time_slice)) +
  geom_violin(scale = "width", width = 0.6, position = position_dodge(width = 0.6), alpha = 0.1) +
  geom_point(aes(y = med),
    data = stats_sep %>% filter(nombre != "samplingProportionEpi") %>% filter(nombre == "ReMultiplier"),
    position = position_dodge(width = 0.6), size = 2.5
  ) +
  scale_color_manual("Time slice",
    labels = skyline_names,
    values = skyline_colors, guide = "none"
  ) +
  scale_fill_manual("Time slice",
    labels = skyline_names,
    values = skyline_colors, guide = "none"
  ) +
  theme_bw() +
  labs(y = TeX("Estimate and 95% HPD"), x = "Cluster", title = TeX("$T_a$")) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )

## SSFrac

SS <- log %>%
  select(contains("SSFrac")) %>%
  mutate(sample = seq_len(n())) %>%
  pivot_longer(!sample, names_to = "nombre", values_to = "value") %>%
  extract(nombre,
    into = c("nombre", "cluster_id", "time_slice"),
    regex = "([A-Za-z]+)(\\d+).(\\d+)",
    convert = TRUE
  ) %>%
  mutate(cluster_id = factor(cluster_id)) %>%
  left_join(stats_sep %>% mutate(cluster_id = as.factor(cluster_id))) %>%
  filter((value > lower) & (value < upper)) %>%
  ggplot() +
  aes(x = factor(cluster_id, cluster_order), y = value, color = as.factor(time_slice), fill = as.factor(time_slice)) +
  geom_violin(scale = "width", width = 0.6, position = position_dodge(width = 0.6), alpha = 0.1) +
  geom_point(aes(y = med),
    data = stats_sep %>% filter(nombre != "samplingProportionEpi") %>% filter(nombre == "SSFrac"),
    position = position_dodge(width = 0.6), size = 2.5
  ) +
  scale_color_manual("Time slice",
    labels = skyline_names,
    values = skyline_colors, guide = "none"
  ) +
  scale_fill_manual("Time slice",
    labels = skyline_names,
    values = skyline_colors, guide = "none"
  ) +
  theme_bw() +
  labs(y = TeX("Estimate and 95% HPD"), x = "Cluster", title = TeX("$f_{ss}$")) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )

# bayes factors

ss_frac_alpha <- mean(log$SSfracAlpha)
ss_frac_beta <- 5
epsilon <- 0.001

den <- pbeta(epsilon, shape1 = ss_frac_alpha, shape2 = ss_frac_beta) / (1 - pbeta(epsilon, shape1 = ss_frac_alpha, shape2 = ss_frac_beta))

cond_probs <- lapply(names(log), function(param_name) {
  x <- log[[param_name]]


  m <- mean(x < epsilon)
  num <- mean(x < epsilon) / mean(x >= epsilon)

  n <- mean(x)

  data.frame(
    param = param_name,
    n = n,
    med = mean(x),
    bf = den / num
  )
}) %>%
  bind_rows() %>%
  extract(param,
    into = c("nombre", "cluster_id", "time_slice"),
    regex = "([A-Za-z]+)(\\d+).(\\d+)",
    convert = TRUE
  ) %>%
  filter(nombre == "SSFrac") %>%
  mutate(
    ss = case_when(bf > 1 ~ "No", T ~ "Yes")
  )

bf <- cond_probs %>%
  mutate(bf2 = base::log(ifelse(bf < 85, bf, exp(5)))) %>%
  ggplot() +
  aes(x = factor(cluster_id, cluster_order), yend = bf2, color = as.factor(time_slice)) +
  geom_segment(aes(y = 0), position = position_dodge(width = 0.5), linewidth = 4, lineend = "round", alpha = 0.9) +
  theme_bw() +
  scale_color_manual("Time slice",
    labels = skyline_names,
    values = skyline_colors
  ) +
  labs(
    x = "Cluster",
    y = "Log(BF)"
  )

# plot to grob the legend
parameters <- stats_sep %>%
  filter(nombre %in% c("ReMultiplier", "Re", "SSFrac")) %>%
  ggplot() +
  aes(x = factor(cluster_id, cluster_order), y = med, color = as.factor(time_slice)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.5) +
  geom_segment(aes(y = upper, yend = lower),
    alpha = 0.5, position = position_dodge(width = 0.6),
    linewidth = 2.5, lineend = "round"
  ) +
  scale_color_manual("Time slice",
    labels = c("Before lockdown", "After lockdown"),
    values = c("#289BFF", "#F28E2B")
  ) +
  theme_bw() +
  labs(y = TeX("Estimate and 95% HPD"), x = "Cluster") +
  facet_wrap(~nombre, ncol = 1, scales = "free_y", labeller = labeller(
    nombre = c(
      "Re" = TeX("$R_e$", output = "character"),
      "ReMultiplier" = TeX("$T_a$", output = "character"),
      "SSFrac" = TeX("$f_{ss}$", output = "character")
    )
  )) +
  theme(legend.position = "bottom")

# MAIN FIGURE
ggarrange(
  ggarrange(
    ggarrange(change_posterior + scale_y_discrete(expand = expansion(add = c(0.35, 0.35)), labels = c(TeX("$R_e$"), TeX("$T_a$"), TeX("$f_{ss}$")))),
    ggarrange(SS, ReM, Re, ncol = 1),
    common.legend = T,
    legend = "bottom",
    widths = c(1, 2),
    align = "hv",
    labels = c("A)", "B)")
  ),
  ggarrange(bf + theme(legend.position = "none")),
  labels = c("", "C)"),
  ncol = 1,
  heights = c(0.9, 0.3),
  align = "v",
  common.legend = T,
  legend.grob = get_legend(parameters),
  legend = "bottom"
)


# FIGURE 5 ####
# Get reproductive numbers
log_m <- log

for (x in cluster_order) {
  # first interval
  mult <- sprintf("ReMultiplier%s.1", x)
  frac <- sprintf("SSFrac%s.1", x)
  re <- sprintf("Re%s.1", x)

  # New names
  ress_ss <- sprintf("ReSS_SS%s.1", x)
  rens_ns <- sprintf("ReNS_NS%s.1", x)
  ress_ns <- sprintf("ReSS_NS%s.1", x)
  rens_ss <- sprintf("ReNS_SS%s.1", x)

  # Operations
  log_m[[rens_ns]] <- ((log_m[[re]] / (1 + log_m[[frac]] * (log_m[[mult]] - 1)))) * (1 - log[[frac]])
  log_m[[ress_ss]] <- ((log_m[[re]] / (1 + log_m[[frac]] * (log_m[[mult]] - 1)))) * log_m[[mult]] * log[[frac]]
  log_m[[rens_ss]] <- ((log_m[[re]] / (1 + log_m[[frac]] * (log_m[[mult]] - 1)))) * (log[[frac]])
  log_m[[ress_ns]] <- ((log_m[[re]] / (1 + log_m[[frac]] * (log_m[[mult]] - 1)))) * log_m[[mult]] * (1 - log[[frac]])

  # Second_slice
  mult2 <- sprintf("ReMultiplier%s.2", x)
  frac2 <- sprintf("SSFrac%s.2", x)
  re2 <- sprintf("Re%s.2", x)

  # New names
  ress_ss2 <- sprintf("ReSS_SS%s.2", x)
  rens_ns2 <- sprintf("ReNS_NS%s.2", x)
  ress_ns2 <- sprintf("ReSS_NS%s.2", x)
  rens_ss2 <- sprintf("ReNS_SS%s.2", x)


  # Operations
  log_m[[rens_ns2]] <- ((log_m[[re2]] / (1 + log_m[[frac2]] * (log_m[[mult2]] - 1)))) * (1 - log[[frac2]])
  log_m[[ress_ss2]] <- ((log_m[[re2]] / (1 + log_m[[frac2]] * (log_m[[mult2]] - 1)))) * log_m[[mult2]] * log[[frac2]]
  log_m[[rens_ss2]] <- ((log_m[[re2]] / (1 + log_m[[frac2]] * (log_m[[mult2]] - 1)))) * (log[[frac2]])
  log_m[[ress_ns2]] <- ((log_m[[re2]] / (1 + log_m[[frac2]] * (log_m[[mult2]] - 1)))) * log_m[[mult2]] * (1 - log[[frac2]])
}

## HPDs data
stats_res <- lapply(names(log_m), function(param_name) {
  x <- log_m[[param_name]]

  m <- mean(x)
  hpd <- HPDinterval(as.mcmc(x), prob = 0.95)

  data.frame(
    param = param_name,
    med = m,
    lower = hpd[1],
    upper = hpd[2]
  )
}) %>%
  bind_rows() %>%
  extract(param,
    into = c("nombre", "cluster_id", "time_slice"),
    regex = "([A-Za-z_]+)(\\d+).(\\d+)",
    convert = TRUE
  ) %>%
  filter(nombre %in% c("ReSS_SS", "ReSS_NS", "ReNS_SS", "ReNS_NS", "Re")) %>%
  mutate(cluster_id = factor(cluster_id, levels = cluster_order))


df_long <- log_m %>%
  select(contains("Re")) %>%
  pivot_longer(everything(), names_to = "id", values_to = "estimate") %>%
  extract(id,
    into = c("nombre", "cluster_id", "time_slice"),
    regex = "([A-Za-z_]+)(\\d+).(\\d+)",
    convert = TRUE
  ) %>%
  mutate(cluster_id = factor(cluster_id)) %>%
  left_join(stats_res) %>%
  filter((estimate > lower) & (estimate < upper))

## nsns
plot_nsns <- df_long %>%
  filter(nombre != "Re") %>%
  filter(!(cluster_id == "1177" & time_slice == "1")) %>%
  filter(nombre == "ReNS_NS") %>%
  ggplot() +
  aes(x = factor(cluster_id, cluster_order), y = estimate, color = factor(time_slice)) +
  geom_violin(aes(fill = factor(time_slice)), position = position_dodge(width = 0.6), alpha = 0.1, scale = "width", width = 0.6) +
  geom_point(aes(y = med),
    size = 2.5, position = position_dodge(width = 0.6),
    data = stats_res %>% filter(nombre != "Re") %>%
      filter(!(cluster_id == "1177" & time_slice == "1")) %>% filter(nombre == "ReNS_NS")
  ) +
  scale_color_manual("Time slice",
    labels = skyline_names,
    values = skyline_colors, guide = "none"
  ) +
  scale_fill_manual("Time slice",
    labels = skyline_names,
    values = skyline_colors, guide = "none"
  ) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(
    x = "Cluster",
    y = "Estimate and 95% HPD",
    title = TeX("$Re_{ns,ns}$")
  ) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )

## ssss
plot_ssss <- df_long %>%
  filter(nombre != "Re") %>%
  filter(!(cluster_id == "1177" & time_slice == "1")) %>%
  filter(nombre == "ReSS_SS") %>%
  ggplot() +
  aes(x = factor(cluster_id, cluster_order), y = estimate, color = factor(time_slice)) +
  geom_violin(aes(fill = factor(time_slice)), position = position_dodge(width = 0.6), alpha = 0.1, scale = "width", width = 0.6) +
  geom_point(aes(y = med),
    size = 2.5, position = position_dodge(width = 0.6),
    data = stats_res %>% filter(nombre != "Re") %>%
      filter(!(cluster_id == "1177" & time_slice == "1")) %>% filter(nombre == "ReSS_SS")
  ) +
  scale_color_manual("Time slice",
    labels = skyline_names,
    values = skyline_colors, guide = "none"
  ) +
  scale_fill_manual("Time slice",
    labels = skyline_names,
    values = skyline_colors, guide = "none"
  ) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(
    x = "Cluster",
    y = "Estimate and 95% HPD",
    title = TeX("$Re_{ss,ss}$")
  ) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )


## nsss
plot_nsss <- df_long %>%
  filter(nombre != "Re") %>%
  filter(!(cluster_id == "1177" & time_slice == "1")) %>%
  filter(nombre == "ReNS_SS") %>%
  ggplot() +
  aes(x = factor(cluster_id, cluster_order), y = estimate, color = factor(time_slice)) +
  geom_violin(aes(fill = factor(time_slice)), position = position_dodge(width = 0.6), alpha = 0.1, scale = "width", width = 0.6) +
  geom_point(aes(y = med),
    size = 2.5, position = position_dodge(width = 0.6),
    data = stats_res %>% filter(nombre != "Re") %>%
      filter(!(cluster_id == "1177" & time_slice == "1")) %>% filter(nombre == "ReNS_SS")
  ) +
  scale_color_manual("Time slice",
    labels = skyline_names,
    values = skyline_colors, guide = "none"
  ) +
  scale_fill_manual("Time slice",
    labels = skyline_names,
    values = skyline_colors, guide = "none"
  ) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(
    x = "Cluster",
    y = "Estimate and 95% HPD",
    title = TeX("$Re_{ns,ss}$")
  ) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )

## ssns
plot_ssns <- df_long %>%
  filter(nombre != "Re") %>%
  filter(!(cluster_id == "1177" & time_slice == "1")) %>%
  filter(nombre == "ReSS_NS") %>%
  ggplot() +
  aes(x = factor(cluster_id, cluster_order), y = estimate, color = factor(time_slice)) +
  geom_violin(aes(fill = factor(time_slice)), position = position_dodge(width = 0.6), alpha = 0.1, scale = "width", width = 0.6) +
  geom_point(aes(y = med),
    size = 2.5, position = position_dodge(width = 0.6),
    data = stats_res %>% filter(nombre != "Re") %>%
      filter(!(cluster_id == "1177" & time_slice == "1")) %>% filter(nombre == "ReSS_NS")
  ) +
  scale_color_manual("Time slice",
    labels = skyline_names,
    values = skyline_colors, guide = "none"
  ) +
  scale_fill_manual("Time slice",
    labels = skyline_names,
    values = skyline_colors, guide = "none"
  ) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(
    x = "Cluster",
    y = "Estimate and 95% HPD",
    title = TeX("$Re_{ss,ns}$")
  ) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )

# mook plot for legend
legend_plot <- df_long %>%
  filter(nombre != "Re") %>%
  filter(!(cluster_id == "1177" & time_slice == "1")) %>%
  filter(nombre == "ReSS_NS") %>%
  ggplot() +
  aes(x = factor(cluster_id, cluster_order), y = estimate, color = factor(time_slice)) +
  geom_violin(aes(fill = factor(time_slice)), position = position_dodge(width = 0.6), alpha = 0.1, scale = "width", width = 0.6) +
  geom_point(aes(y = med),
    size = 2.5, position = position_dodge(width = 0.6),
    data = stats_res %>% filter(nombre != "Re") %>%
      filter(!(cluster_id == "1177" & time_slice == "1")) %>% filter(nombre == "ReSS_NS")
  ) +
  scale_color_manual("Time slice",
    labels = skyline_names,
    values = skyline_colors
  ) +
  scale_fill_manual("Time slice",
    labels = skyline_names,
    values = skyline_colors
  ) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(
    x = "Cluster",
    y = "Estimate and 95% HPD",
    title = TeX("$Re_{ss,ns}$")
  ) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )
# subplot reproductive numbers
subplot_res <- ggarrange(
  plot_nsns, plot_nsss,
  plot_ssns, plot_ssss,
  ncol = 2,
  nrow = 2,
  common.legend = TRUE,
  legend.grob = get_legend(legend_plot),
  legend = "bottom"
)

# subplot growts
probs_dynamics <- log_m %>%
  select(contains("Re")) %>%
  pivot_longer(everything(), names_to = "id", values_to = "estimate") %>%
  extract(id,
    into = c("nombre", "cluster_id", "time_slice"),
    regex = "([A-Za-z_]+)(\\d+).(\\d+)",
    convert = TRUE
  ) %>%
  filter(nombre != "ReMultiplier") %>%
  group_by(cluster_id, nombre, time_slice) %>%
  mutate(sample = seq_len(n())) %>%
  ungroup() %>%
  pivot_wider(names_from = nombre, values_from = estimate) %>%
  group_by(cluster_id, time_slice) %>%
  filter(mean(Re) > 1) %>%
  summarise(
    growth_ns = mean(Re > 1 & ReNS_NS > 1 & ReSS_SS < 1) / mean(Re > 1),
    growth_ss = mean(Re > 1 & ReNS_NS < 1 & ReSS_SS > 1) / mean(Re > 1),
    growth_inter = mean(Re > 1 & ReNS_NS < 1 & ReSS_SS < 1) / mean(Re > 1),
    growth_both = mean(Re > 1 & ReNS_NS > 1 & ReSS_SS > 1) / mean(Re > 1),
    dec_ns = mean(Re < 1 & ReNS_NS > 1 & ReSS_SS < 1) / mean(Re < 1),
    dec_ss = mean(Re < 1 & ReNS_NS < 1 & ReSS_SS > 1) / mean(Re < 1),
    dec_inter = mean(Re < 1 & ReNS_NS < 1 & ReSS_SS < 1) / mean(Re < 1),
    dec_both = mean(Re < 1 & ReNS_NS > 1 & ReSS_SS > 1) / mean(Re < 1)
  )


subplot_growths <- probs_dynamics %>%
  select(!contains("dec")) %>%
  left_join(cond_probs %>% select(cluster_id, time_slice, bf)) %>%
  # filter(bf < 1) %>%
  select(!bf) %>%
  pivot_longer(c(contains("growth"), contains("dec")), names_to = "type", values_to = "prob") %>%
  filter(!(cluster_id == "1177" & time_slice == 1)) %>%
  mutate(
    prob = case_when(time_slice == 1 ~ -prob, T ~ prob)
  ) %>%
  ggplot() +
  aes(x = prob, y = factor(cluster_id, cluster_order), fill = factor(type, c("growth_inter", "growth_ns", "growth_ss", "growth_both"))) +
  geom_col(position = "stack") +
  theme_bw() +
  geom_vline(aes(xintercept = 0), linetype = 5, color = "black", linewidth = 1.5) +
  labs(
    fill = "Dynamic", y = "Cluster",
    x = "Pr before the lockdown                                                                       Pr after the lockdown"
  ) +
  scale_fill_manual(
    values = c("growth_both" = "#615B9A", "growth_ns" = "#C1667D", "growth_ss" = "#5F2644", "growth_inter" = "#E9A6A9"),
    labels = c("Synergic growth", "NS driven", "SS driven", "Complete growth")
  ) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1), labels = abs) +
  theme(legend.position = "bottom")


# MAIN FIGURE
ggarrange(
  subplot_res,
  subplot_growths,
  ncol = 1
)
