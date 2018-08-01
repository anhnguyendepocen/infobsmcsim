### ---
# Packages, options, and stuff
### ---

# Packages
library(tidyverse)
library(formattable)
library(haven)
library(labelled)
library(rsimsum)
library(cowplot)
library(knitr)
library(kableExtra)
library(pracma)

# Number of significant digits
sigd <- 4

# comma2 function for labelling
comma2 <- function(x, ...) formattable::comma(x = x, digits = sigd, ...)

### ---
# Some summary statistics of the simulated data
### ---

# Summarise average number of observation and average gap time per DGM
# sim1_dl.RDS is a file containing all simulated dataset
if (!file.exists("sim1_dl.RDS")) {
  df_list <- lapply(list.files(path = "Data/", full.names = TRUE), read_dta)
  df_list <- bind_rows(df_list)
  df_list <- split(df_list, f = df_list[["dgm"]])
  saveRDS(df_list, file = "sim1_dl.RDS")
  rm(df_list)
  gc()
}
df_list <- readRDS("sim1_dl.RDS")
df <- bind_rows(df_list)
dgm <- distinct(df, dgm, gamma, lambda) %>%
  mutate(
    dgm = as.character(dgm),
    gamma = as.character(comma(gamma, 2)),
    lambda = as.character(comma(lambda, 2))
  )

df %>%
  group_by(dgm, i) %>%
  mutate(ss = n()) %>%
  ungroup() %>%
  distinct(dgm, gamma, lambda, i, .keep_all = TRUE) %>%
  group_by(lambda) %>%
  summarise(
    min_ss = min(ss),
    q1_ss = quantile(ss, probs = 0.25),
    mean_ss = mean(ss),
    median_ss = median(ss),
    q3_ss = quantile(ss, probs = 0.75),
    max_ss = max(ss)
  ) %>%
  View(title = "ss")

df %>%
  distinct(dgm, gamma, lambda, i, id, .keep_all = TRUE) %>%
  group_by(lambda) %>%
  summarise(
    min_nn = min(nn),
    q1_nn = quantile(nn, probs = 0.25),
    mean_nn = mean(nn),
    median_nn = median(nn),
    q3_nn = quantile(nn, probs = 0.75),
    max_nn = max(nn)
  ) %>%
  View(title = "nn")

df %>%
  filter(ind == 1) %>%
  group_by(lambda) %>%
  summarise(
    min_obtime = min(obtime),
    q1_obtime = quantile(obtime, probs = 0.25),
    mean_obtime = mean(obtime),
    median_obtime = median(obtime),
    q3_obtime = quantile(obtime, probs = 0.75),
    max_obtime = max(obtime)
  ) %>%
  View(title = "obtime")

# Summarise marginal survival for given values of lambda
ms <- function(t, lambda, v, trt, beta = 1, sigmau = 1) {
  S <- exp(-lambda * exp(beta * trt) * t^v)
  fn <- function(eta) S^exp(eta) * dnorm(eta, mean = 0, sd = sigmau)
  ms <- quadinf(f = fn, xa = -Inf, xb = Inf)$Q
  return(ms)
}
msdf <- crossing(
  t = seq(0.01, 15, by = 0.01),
  lambda = c(0.10, 0.30, 1.00),
  v = 1.05,
  trt = 0:1
)
msdf$ms <- vapply(1:nrow(msdf), function(i) with(msdf[i, ], ms(t, lambda, v, trt)), numeric(1))
msdf <- split(x = msdf, f = lapply(c("lambda", "trt"), function(f) msdf[[f]]), sep = "#")
msfun <- lapply(msdf, function(x) with(x, splinefun(x = ms, y = t, method = "natural")))
tmedian <- lapply(seq_along(msfun), function(i) data.frame(
    t = msfun[[i]](0.5),
    what = names(msfun)[i]
  ))
tmedian <- do.call(rbind.data.frame, tmedian) %>%
  separate(what, into = c("lambda", "trt"), sep = "#") %>%
  spread(key = trt, value = t) %>%
  rename(t0 = `0`, t1 = `1`)
View(tmedian)

### ---
# Simulation results
### ---
file <- "sim1.dta"
sim1 <- read_dta(file)
attr(sim1, "spec") <- NULL
sim1 <- labelled::remove_labels(sim1, user_na_to_na = TRUE)
sim1 <- filter(sim1, est %in% 1:3, model %in% c(1:3, 5:6)) %>%
  mutate(
    model = factor(model, levels = c(1:3, 5:6), labels = c("True", "ME (Total)", "ME (Cumulative)", "ME (No adjustment)", "GEE (IIVW)")),
    est = factor(est, levels = 1:3, labels = c("a0", "a1", "a2"))
  ) %>%
  rename(par = est)

# Make multisimsum object
ms <- multisimsum(sim1, par = "par", true = c("a0" = 0, "a1" = 1, "a2" = 0.2), estvarname = "b", se = "se", methodvar = "model", by = "dgm", ref = "True")

# DF with results
res <- get_data(summary(ms)) %>%
  filter(stat %in% c("thetamean", "se2mean", "empse", "bias", "cover", "mse")) %>%
  mutate(est = ifelse(stat == "se2mean", sqrt(est), est))

# Summary tables
tres <- res %>%
  mutate(value = ifelse(stat %in% c("thetamean", "se2mean"), paste0(comma(est, sigd)), paste0(comma(est, sigd), " (", comma(mcse, sigd), ")"))) %>%
  select(-est, -mcse, -lower, -upper) %>%
  arrange(dgm) %>%
  spread(key = stat, value = value) %>%
  mutate(
    par = case_when(
      par == "a0" ~ "\\(\\alpha_0\\)",
      par == "a1" ~ "\\(\\alpha_1\\)",
      par == "a2" ~ "\\(\\alpha_2\\)"
    ),
    model = factor(model, levels = c("True", "ME (Total)", "ME (Cumulative)", "ME (No adjustment)", "GEE (IIVW)")),
  ) %>%
  select(dgm, model, par, thetamean, se2mean, empse, bias, cover, mse) %>%
  arrange(dgm, model, par)

for (i in as.character(4:6)) {
  sink(paste0("sim1_files/tab-dgm-", i, ".tex"))
  filter(tres, dgm == i) %>%
    select(-dgm, -model) %>%
    kable(
      format = "latex",
      booktabs = TRUE,
      align = rep("r", ncol(res) - 2),
      col.names = c("", "Estimate", "SE", "Empirical SE", "Bias", "Coverage", "MSE"),
      caption = paste0("\\label{tab:dgm-", i, "}Results for data-generating mechanism with \\(\\gamma\\) = ", dgm[["gamma"]][dgm[["dgm"]] == i], " and \\(\\lambda\\) = ", dgm[["lambda"]][dgm[["dgm"]] == i]),
      escape = FALSE
    ) %>%
    kable_styling(latex_options = c("scale_down", "hold_position")) %>%
    group_rows(paste0("Model: ", levels(tres$model)[1]), 1, 3) %>%
    group_rows(paste0("Model: ", levels(tres$model)[2]), 4, 6) %>%
    group_rows(paste0("Model: ", levels(tres$model)[3]), 7, 9) %>%
    group_rows(paste0("Model: ", levels(tres$model)[4]), 10, 12) %>%
    group_rows(paste0("Model: ", levels(tres$model)[5]), 13, 15) %>%
    print()
  sink()
  tmpf <- read_lines(file = paste0("sim1_files/tab-dgm-", i, ".tex"))
  tmpf <- str_replace_all(string = tmpf, pattern = "\\\\label\\{tab:\\}", replacement = "")
  write_lines(x = tmpf, path = paste0("sim1_files/tab-dgm-", i, ".tex"))
}

# Make a plot for bias, coverage, MSE
pres <- res %>%
  left_join(dgm, "dgm") %>%
  mutate(
    gamma = factor(gamma),
    lambda = factor(lambda),
    model = factor(model, levels = c("True", "ME (Total)", "ME (Cumulative)", "ME (No adjustment)", "GEE (IIVW)"), labels = LETTERS[1:5])
  ) %>%
  mutate(sig = case_when(
    stat == "bias" & (lower > 0 | upper < 0) ~ "*",
    stat == "cover" & (lower > 0.95 | upper < 0.95) ~ "*",
    TRUE ~ ""
  ))

bias_a0 <- pres %>%
  filter(stat == "bias" & par == "a0") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  geom_text(aes(y = min(est) - 0.1, label = sig), position = position_dodge(width = 1 / 3), show.legend = FALSE) +
  scale_color_grey() +
  scale_y_continuous(labels = comma2) +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = expression(Bias ~ of ~ alpha[0]), shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
bias_a1 <- pres %>%
  filter(stat == "bias" & par == "a1") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  geom_text(aes(y = max(est) + 0.05, label = sig), position = position_dodge(width = 1 / 3), show.legend = FALSE) +
  scale_color_grey() +
  scale_y_continuous(labels = comma2) +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = expression(Bias ~ of ~ alpha[1]), shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
bias_a2 <- pres %>%
  filter(stat == "bias" & par == "a2") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  geom_text(aes(y = min(est) - 0.01, label = sig), position = position_dodge(width = 1 / 3), show.legend = FALSE) +
  scale_color_grey() +
  scale_y_continuous(labels = comma2) +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = expression(Bias ~ of ~ alpha[2]), shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
bias <- plot_grid(
  plot_grid(bias_a0 + theme(legend.position = "none"),
    bias_a1 + theme(legend.position = "none"),
    bias_a2 + theme(legend.position = "none"),
    ncol = 1,
    labels = LETTERS,
    align = "vh", axis = "ltrb"
  ),
  get_legend(bias_a2),
  ncol = 1,
  rel_heights = c(9, 1)
)

ggsave(bias, filename = "sim1_files/bias.pdf", dpi = 1200, height = 9, width = 6)

cover_a0 <- pres %>%
  filter(stat == "cover" & par == "a0") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0.95, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  geom_text(aes(y = 1, label = sig), position = position_dodge(width = 1 / 3), show.legend = FALSE) +
  scale_color_grey() +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = expression(Coverage ~ of ~ alpha[0]), shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
cover_a1 <- pres %>%
  filter(stat == "cover" & par == "a1") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0.95, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  geom_text(aes(y = 1, label = sig), position = position_dodge(width = 1 / 3), show.legend = FALSE) +
  scale_color_grey() +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = expression(Coverage ~ of ~ alpha[1]), shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
cover_a2 <- pres %>%
  filter(stat == "cover" & par == "a2") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0.95, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  geom_text(aes(y = 1, label = sig), position = position_dodge(width = 1 / 3), show.legend = FALSE) +
  scale_color_grey() +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = expression(Coverage ~ of ~ alpha[2]), shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
cover <- plot_grid(
  plot_grid(cover_a0 + theme(legend.position = "none"),
    cover_a1 + theme(legend.position = "none"),
    cover_a2 + theme(legend.position = "none"),
    ncol = 1,
    labels = LETTERS,
    align = "vh", axis = "ltrb"
  ),
  get_legend(cover_a2),
  ncol = 1,
  rel_heights = c(9, 1)
)

ggsave(cover, filename = "sim1_files/cover.pdf", dpi = 1200, height = 9, width = 6)

mse_a0 <- pres %>%
  filter(stat == "mse" & par == "a0") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  scale_color_grey() +
  scale_y_continuous(labels = comma2) +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = expression(MSE ~ of ~ alpha[0]), shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
mse_a1 <- pres %>%
  filter(stat == "mse" & par == "a1") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  scale_color_grey() +
  scale_y_continuous(labels = comma2) +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = expression(MSE ~ of ~ alpha[1]), shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
mse_a2 <- pres %>%
  filter(stat == "mse" & par == "a2") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  scale_color_grey() +
  scale_y_continuous(labels = comma2) +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = expression(MSE ~ of ~ alpha[2]), shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
mse <- plot_grid(
  plot_grid(mse_a0 + theme(legend.position = "none"),
    mse_a1 + theme(legend.position = "none"),
    mse_a2 + theme(legend.position = "none"),
    ncol = 1,
    labels = LETTERS,
    align = "vh", axis = "ltrb"
  ),
  get_legend(mse_a2),
  ncol = 1,
  rel_heights = c(9, 1)
)

ggsave(mse, filename = "sim1_files/mse.pdf", dpi = 1200, height = 9, width = 6)

### ---
# Plots for the application
### ---

# ex-res.dta contains predicted values from the joint model
ex_data <- haven::read_dta(file = "sim1_files/ex-res.dta")
ex_data <- dplyr::mutate(ex_data, sex = factor(sex, levels = 1:2, labels = c("Males", "Females")))
ex_p1 <- ex_data %>%
  ggplot(aes(x = tY, fill = sex)) +
  geom_ribbon(aes(ymin = egfr_me_lci, ymax = egfr_me_uci, lty = "Mixed Effects Model"), alpha = 1 / 4) +
  geom_ribbon(aes(ymin = egfr_jm_lci, ymax = egfr_jm_uci, lty = "Joint Model"), alpha = 1 / 4) +
  geom_line(aes(y = egfr_me, color = sex, lty = "Mixed Effects Model")) +
  geom_line(aes(y = egfr_jm, color = sex, lty = "Joint Model")) +
  scale_color_grey() +
  scale_fill_grey() +
  guides(fill = FALSE) +
  theme_bw() +
  coord_cartesian(ylim = c(53, 63)) +
  theme(legend.position = "bottom") +
  labs(x = "Follow-up time (Years)", y = "Predicted eGFR", colour = "", fill = "", lty = "")
ggsave(ex_p1, filename = "sim1_files/ex_p1.pdf", dpi = 1200, height = 4, width = 6)

### ---
#   Supplementary material
### ---

# Tabs for gamma = 0
for (i in as.character(1:3)) {
  sink(paste0("sim1_files/tab-dgm-", i, ".tex"))
  filter(tres, dgm == i) %>%
    select(-dgm, -model) %>%
    kable(
      format = "latex",
      booktabs = TRUE,
      align = rep("r", ncol(res) - 2),
      col.names = c("", "Estimate", "SE", "Empirical SE", "Bias", "Coverage", "MSE"),
      caption = paste0("\\label{tab:dgm-", i, "}Results for data-generating mechanism with \\(\\gamma\\) = ", dgm[["gamma"]][dgm[["dgm"]] == i], " and \\(\\lambda\\) = ", dgm[["lambda"]][dgm[["dgm"]] == i]),
      escape = FALSE
    ) %>%
    kable_styling(latex_options = c("scale_down", "hold_position")) %>%
    group_rows(paste0("Model: ", levels(tres$model)[1]), 1, 3) %>%
    group_rows(paste0("Model: ", levels(tres$model)[2]), 4, 6) %>%
    group_rows(paste0("Model: ", levels(tres$model)[3]), 7, 9) %>%
    group_rows(paste0("Model: ", levels(tres$model)[4]), 10, 12) %>%
    group_rows(paste0("Model: ", levels(tres$model)[5]), 13, 15) %>%
    print()
  sink()
  tmpf <- read_lines(file = paste0("sim1_files/tab-dgm-", i, ".tex"))
  tmpf <- str_replace_all(string = tmpf, pattern = "\\\\label\\{tab:\\}", replacement = "")
  write_lines(x = tmpf, path = paste0("sim1_files/tab-dgm-", i, ".tex"))
}

# Bias and coverage of variances (residual, random effects)
sim1_var <- read_dta(file)
attr(sim1_var, "spec") <- NULL
sim1_var <- labelled::remove_labels(sim1_var, user_na_to_na = TRUE)
sim1_var <- filter(sim1_var, est %in% 4:6 & model %in% c(1:3, 5:6)) %>%
  mutate(
    model = factor(model, levels = c(1:3, 5:6), labels = c("True", "ME (Total)", "ME (Cumulative)", "ME (No adjustment)", "GEE (IIVW)")),
    est = factor(est, levels = 4:6, labels = c("V(u)", "V(v)", "V(e)"))
  ) %>%
  rename(par = est)
ms_var <- multisimsum(sim1_var, par = "par", true = c("V(u)" = 1, "V(v)" = 0.5, "V(e)" = 1), estvarname = "b", se = "se", methodvar = "model", by = "dgm", ref = "True")
pres_var <- get_data(summary(ms_var)) %>%
  filter(stat %in% c("bias", "cover", "mse")) %>%
  left_join(dgm, "dgm") %>%
  mutate(
    gamma = factor(gamma, levels = c("0.00", "1.50")),
    lambda = factor(lambda),
    model = factor(model, levels = c("True", "ME (Total)", "ME (Cumulative)", "ME (No adjustment)", "GEE (IIVW)"), labels = LETTERS[1:5]),
  ) %>%
  mutate(sig = case_when(
    stat == "bias" & (lower > 0 | upper < 0) ~ "*",
    stat == "cover" & (lower > 0.95 | upper < 0.95) ~ "*",
    TRUE ~ ""
  ))
bias_vu <- pres_var %>%
  filter(stat == "bias" & par == "V(u)") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  geom_text(aes(y = -0.01, label = sig), position = position_dodge(width = 1 / 3), show.legend = FALSE) +
  scale_y_continuous(labels = comma2) +
  scale_color_grey() +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = "Bias of V(u)", shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
bias_vv <- pres_var %>%
  filter(stat == "bias" & par == "V(v)") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  geom_text(aes(y = min(est, na.rm = TRUE) - 0.1, label = sig), position = position_dodge(width = 1 / 3), show.legend = FALSE) +
  scale_y_continuous(labels = comma2) +
  scale_color_grey() +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = "Bias of V(v)", shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
bias_ve <- pres_var %>%
  filter(stat == "bias" & par == "V(e)") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  geom_text(aes(y = max(est, na.rm = TRUE) + 0.01, label = sig), position = position_dodge(width = 1 / 3), show.legend = FALSE) +
  scale_y_continuous(labels = comma2) +
  scale_color_grey() +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = expression(Bias ~ of ~ V(epsilon)), shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
bias_var <- plot_grid(
  plot_grid(bias_vu + theme(legend.position = "none"),
    bias_vv + theme(legend.position = "none"),
    bias_ve + theme(legend.position = "none"),
    ncol = 1,
    labels = LETTERS,
    align = "vh", axis = "ltrb"
  ),
  get_legend(bias_ve),
  ncol = 1,
  rel_heights = c(9, 1)
)

ggsave(bias_var, filename = "sim1_files/bias_var.pdf", dpi = 1200, height = 9, width = 6)

cover_vu <- pres_var %>%
  filter(stat == "cover" & par == "V(u)") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0.95, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  geom_text(aes(y = 1, label = sig), position = position_dodge(width = 1 / 3), show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_grey() +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = "Coverage of V(u)", shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
cover_vv <- pres_var %>%
  filter(stat == "cover" & par == "V(v)") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0.95, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  geom_text(aes(y = 1, label = sig), position = position_dodge(width = 1 / 3), show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_grey() +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = "Coverage of V(v)", shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
cover_ve <- pres_var %>%
  filter(stat == "cover" & par == "V(e)") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0.95, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  geom_text(aes(y = 1, label = sig), position = position_dodge(width = 1 / 3), show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_grey() +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = expression(Coverage ~ of ~ V(epsilon)), shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
cover_var <- plot_grid(
  plot_grid(cover_vu + theme(legend.position = "none"),
    cover_vv + theme(legend.position = "none"),
    cover_ve + theme(legend.position = "none"),
    ncol = 1,
    labels = LETTERS,
    align = "vh", axis = "ltrb"
  ),
  get_legend(cover_ve),
  ncol = 1,
  rel_heights = c(9, 1)
)

ggsave(cover_var, filename = "sim1_files/cover_var.pdf", dpi = 1200, height = 9, width = 6)

mse_vu <- pres_var %>%
  filter(stat == "bias" & par == "V(u)") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  scale_y_continuous(labels = comma2) +
  scale_color_grey() +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = "MSE of V(u)", shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
mse_vv <- pres_var %>%
  filter(stat == "bias" & par == "V(v)") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  scale_y_continuous(labels = comma2) +
  scale_color_grey() +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = "MSE of V(v)", shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
mse_ve <- pres_var %>%
  filter(stat == "mse" & par == "V(e)") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  scale_y_continuous(labels = comma2) +
  scale_color_grey() +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = expression(MSE ~ of ~ V(epsilon)), shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
mse_var <- plot_grid(
  plot_grid(mse_vu + theme(legend.position = "none"),
    mse_vv + theme(legend.position = "none"),
    mse_ve + theme(legend.position = "none"),
    ncol = 1,
    labels = LETTERS,
    align = "vh", axis = "ltrb"
  ),
  get_legend(mse_ve),
  ncol = 1,
  rel_heights = c(9, 1)
)

ggsave(mse_var, filename = "sim1_files/mse_var.pdf", dpi = 1200, height = 9, width = 6)

# Summary statistics for gamma
sim1_gamma <- read_dta(file)
attr(sim1_gamma, "spec") <- NULL
sim1_gamma <- labelled::remove_labels(sim1_gamma, user_na_to_na = TRUE)
sim1_gamma <- filter(sim1_gamma, est == 7 & model %in% c(1:3, 5:6)) %>%
  mutate(
    model = factor(model, levels = c(1:3, 5:6), labels = c("True", "ME (Total)", "ME (Cumulative)", "ME (No adjustment)", "GEE (IIVW)")),
    est = ifelse(dgm %in% 1:3, "gamma.0", "gamma.1.5")
  ) %>%
  rename(par = est)
ms_gamma <- multisimsum(sim1_gamma, par = "par", true = c("gamma.0" = 0, "gamma.1.5" = 1.5), estvarname = "b", se = "se", methodvar = "model", by = "dgm", ref = "True")
pres_gamma <- get_data(summary(ms_gamma)) %>%
  filter(stat %in% c("bias", "cover", "mse")) %>%
  left_join(dgm, "dgm") %>%
  mutate(
    gamma = factor(gamma, levels = c("0.00", "1.50")),
    lambda = factor(lambda),
    model = factor(model, levels = c("True", "ME (Total)", "ME (Cumulative)", "ME (No adjustment)", "GEE (IIVW)"), labels = LETTERS[1:5]),
  ) %>%
  mutate(sig = case_when(
    stat == "bias" & (lower > 0 | upper < 0) ~ "*",
    stat == "cover" & (lower > 0.95 | upper < 0.95) ~ "*",
    TRUE ~ ""
  )) %>%
  mutate(par = "gamma")
bias_gamma <- pres_gamma %>%
  filter(stat == "bias" & par == "gamma") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  geom_text(aes(y = min(est) - 0.01, label = sig), position = position_dodge(width = 1 / 3), show.legend = FALSE) +
  scale_y_continuous(labels = comma2) +
  scale_color_grey() +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = bquote("Bias" ~ "of" ~ gamma), shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
cover_gamma <- pres_gamma %>%
  filter(stat == "cover" & par == "gamma") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0.95, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  geom_text(aes(y = 1, label = sig), position = position_dodge(width = 1 / 3), show.legend = FALSE) +
  scale_y_continuous(labels = comma2) +
  scale_color_grey() +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = bquote("Coverage" ~ "of" ~ gamma), shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
mse_gamma <- pres_gamma %>%
  filter(stat == "mse" & par == "gamma") %>%
  ggplot(aes(x = model, y = est, ymin = lower, ymax = upper, shape = gamma, linetype = gamma, colour = lambda)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dotted") +
  geom_point(position = position_dodge(width = 1 / 3)) +
  geom_text(aes(y = -0.01, label = sig), position = position_dodge(width = 1 / 3), show.legend = FALSE) +
  scale_y_continuous(labels = comma2) +
  scale_color_grey() +
  theme_bw() +
  theme(plot.margin = unit(c(1, 0, 1, 0), "pt"), legend.position = "bottom") +
  labs(x = "Model", y = bquote("MSE" ~ "of" ~ gamma), shape = bquote(lambda ~ ":"), linetype = bquote(lambda ~ ":"), colour = bquote(gamma ~ ":"))
summary_gamma <- plot_grid(
  plot_grid(bias_gamma + theme(legend.position = "none"),
    cover_gamma + theme(legend.position = "none"),
    mse_gamma + theme(legend.position = "none"),
    ncol = 1,
    labels = LETTERS,
    align = "vh", axis = "ltrb"
  ),
  get_legend(bias_gamma),
  ncol = 1,
  rel_heights = c(9, 1)
)

ggsave(summary_gamma, filename = "sim1_files/summary_gamma.pdf", dpi = 1200, height = 9, width = 6)

res_gamma <- get_data(summary(ms_gamma)) %>%
  filter(stat %in% c("thetamean", "se2mean", "empse", "bias", "cover", "mse")) %>%
  mutate(est = ifelse(stat == "se2mean", sqrt(est), est)) %>%
  filter(model == "True") %>%
  select(-model, -par)

tres_gamma <- res_gamma %>%
  mutate(value = ifelse(stat %in% c("thetamean", "se2mean"), paste0(comma(est, sigd)), paste0(comma(est, sigd), " (", comma(mcse, sigd), ")"))) %>%
  select(-est, -mcse, -lower, -upper) %>%
  arrange(dgm) %>%
  spread(key = stat, value = value) %>%
  left_join(dgm, by = "dgm") %>%
  select(gamma, lambda, thetamean, se2mean, empse, bias, cover, mse) %>%
  arrange(gamma, lambda)

sink("sim1_files/tab-gamma.tex")
kable(tres_gamma,
  format = "latex",
  booktabs = TRUE,
  align = rep("r", ncol(tres_gamma)),
  col.names = c("True \\(\\gamma\\)", "True \\(\\lambda\\)", "Estimate", "SE", "Empirical SE", "Bias", "Coverage", "MSE"),
  caption = "\\label{tab:gamma}Results for the association parameter \\(\\gamma\\)",
  escape = FALSE
) %>%
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  print()
sink()
tmpf <- read_lines(file = "sim1_files/tab-gamma.tex")
tmpf <- str_replace_all(string = tmpf, pattern = "\\\\label\\{tab:\\}", replacement = "")
write_lines(x = tmpf, path = "sim1_files/tab-gamma.tex")
