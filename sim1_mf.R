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
library(ggthemes)
library(ggridges)
library(readxl)

# Number of significant digits
sigd <- 4

# comma2 function for labelling
comma2 <- function(x, ...) formattable::comma(x = x, digits = sigd, ...)
comma.2sigd <- function(x, ...) formattable::comma(x = x, digits = 2, ...)
comma.3sigd <- function(x, ...) formattable::comma(x = x, digits = 3, ...)

### ---
# Some summary statistics of the simulated data
### ---

# Summarise average number of observation and average gap time per DGM
# sim1_dl.RDS is a file containing all simulated dataset
if (!file.exists("sim1_dl.RDS")) {
  df <- lapply(list.files(path = "Data/", full.names = TRUE), read_dta)
  df <- dplyr::bind_rows(df)
  df$u <- NULL
  df$v <- NULL
  df$ge <- NULL
  df$iivw <- NULL
  df$xb <- NULL
  df$cumnn <- NULL
  df$cnn <- NULL
  saveRDS(df, file = "sim1_dl.RDS")
  rm(df)
  gc()
}
df <- readRDS("sim1_dl.RDS")
dgm <- data.frame(
  dgm = seq(10),
  lambda = c(rep(c(0.10, 0.30, 1.00), times = 2), rep(NA, times = 3), 0.05),
  gamma = c(rep(c(0.00, 1.50, NA), each = 3), 3.00),
  shape = c(rep(NA, times = 6), rep(2, times = 3), NA),
  theta = c(rep(NA, times = 6), 0, 2, 2, NA)
) %>%
  dplyr::mutate(
    gamma = as.character(formattable::comma(gamma, 2)),
    lambda = as.character(formattable::comma(lambda, 2)),
    shape = as.character(formattable::comma(shape, 2)),
    theta = as.character(formattable::comma(theta, 2))
  ) %>%
  dplyr::mutate(label = dplyr::case_when(
    dgm <= 6 ~ paste0("JM (\\(\\gamma = ", gamma, ", \\lambda = ", lambda, "\\))"),
    dgm %in% 7:8 ~ paste0("\\(\\Gamma\\) distribution ", ifelse(dgm == 7, "not ", ""), "depending on treatment"),
    dgm == 9 ~ paste0("\\(\\Gamma\\) distribution depending on treatment and previous Y"),
    dgm == 10 ~ paste0("JM (\\(\\gamma = ", gamma, ", \\lambda = ", lambda, "\\)) with regular visits")
  )) %>%
  dplyr::mutate(dgm = as.character(dgm))

ss <- df %>%
  dplyr::group_by(dgm, i) %>%
  dplyr::mutate(ss = n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(dgm, i, .keep_all = TRUE) %>%
  dplyr::group_by(dgm) %>%
  dplyr::summarise(
    ss = paste0(as.character(formattable::comma(median(ss), 0)), " (", as.character(formattable::comma(quantile(ss, probs = 0.25), 0)), " - ", as.character(formattable::comma(quantile(ss, probs = 0.75), 0)), ")")
  ) %>%
  dplyr::ungroup()

nn <- df %>%
  dplyr::distinct(dgm, i, id, .keep_all = TRUE) %>%
  dplyr::group_by(dgm) %>%
  dplyr::summarise(
    nn = paste0(as.character(formattable::comma(median(nn), 0)), " (", as.character(formattable::comma(quantile(nn, probs = 0.25), 0)), " - ", as.character(formattable::comma(quantile(nn, probs = 0.75), 0)), ")")
  ) %>%
  dplyr::ungroup()

obtime <- df %>%
  dplyr::filter(ind == 1) %>%
  dplyr::group_by(dgm) %>%
  dplyr::summarise(
    obtime = paste0(as.character(formattable::comma(median(obtime), 2)), " (", as.character(formattable::comma(quantile(obtime, probs = 0.25), 2)), " - ", as.character(formattable::comma(quantile(obtime, probs = 0.75), 2)), ")")
  ) %>%
  dplyr::ungroup()

dgmsumm <- dplyr::left_join(ss, nn, by = "dgm") %>%
  dplyr::left_join(obtime, by = "dgm") %>%
  dplyr::mutate(dgm = as.character(dgm)) %>%
  dplyr::left_join(select(dgm, dgm, label), by = "dgm") %>%
  dplyr::select(label, ss, nn, obtime)

# Table with summaries of simulated data
sink("sim1_files/tab-dgmsumm.tex")
dgmsumm %>%
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    align = rep("r", ncol(dgmsumm)),
    col.names = c("Data-Generating Mechanism", "Sample Size", "N. of Measurements", "Gap Time"),
    caption = "\\label{tab:dgmsumm}Summary characteristics of simulated data under each data-generating mechanism. Values are median with inter-quartile interval [IQI].",
    escape = FALSE,
    linesep = ""
  ) %>%
  kableExtra::kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  print()
sink()
tmpf <- readr::read_lines(file = "sim1_files/tab-dgmsumm.tex")
tmpf <- stringr::str_replace_all(string = tmpf, pattern = "\\\\label\\{tab:\\}", replacement = "")
readr::write_lines(x = tmpf, path = "sim1_files/tab-dgmsumm.tex")

rm(dgmsumm, nn, obtime, ss)
gc()

# Ridgelines plot comparing the distribution of gap times between treated and untreated individuals across DGMs

obtime_rp <- df %>%
  dplyr::filter(ind == 1) %>%
  dplyr::mutate(
    trt = factor(trt, levels = 0:1, labels = c("Untreated Individuals", "Treated Individuals")),
    dgm = factor(dgm)
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = obtime, y = dgm, color = trt, fill = trt)) +
  ggridges::geom_density_ridges(alpha = 0.50) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = c(1, 1), legend.justification = c(1, 1), legend.background = ggplot2::element_blank(), legend.key = ggplot2::element_blank()) +
  ggplot2::labs(x = "Gap time between observations", y = "Data-generating mechanism", fill = "", color = "")
ggplot2::ggsave(obtime_rp, filename = "sim1_files/obtime_rp.pdf", dpi = 1200, height = 6, width = 6)

rm(obtime_rp)
gc()

# Summarise marginal survival for given values of lambda
# ms <- function(t, lambda, v, trt, beta = 1, sigmau = 1) {
#   S <- exp(-lambda * exp(beta * trt) * t^v)
#   fn <- function(eta) S^exp(eta) * dnorm(eta, mean = 0, sd = sigmau)
#   ms <- quadinf(f = fn, xa = -Inf, xb = Inf)$Q
#   return(ms)
# }
# msdf <- crossing(
#   t = seq(0.01, 15, by = 0.01),
#   lambda = c(0.10, 0.30, 1.00),
#   v = 1.05,
#   trt = 0:1
# )
# msdf$ms <- vapply(1:nrow(msdf), function(i) with(msdf[i, ], ms(t, lambda, v, trt)), numeric(1))
# msdf <- split(x = msdf, f = lapply(c("lambda", "trt"), function(f) msdf[[f]]), sep = "#")
# msfun <- lapply(msdf, function(x) with(x, splinefun(x = ms, y = t, method = "natural")))
# tmedian <- lapply(seq_along(msfun), function(i) data.frame(
#     t = msfun[[i]](0.5),
#     what = names(msfun)[i]
#   ))
# tmedian <- do.call(rbind.data.frame, tmedian) %>%
#   separate(what, into = c("lambda", "trt"), sep = "#") %>%
#   spread(key = trt, value = t) %>%
#   rename(t0 = `0`, t1 = `1`)
# View(tmedian)

### ---
# Simulation results
### ---
file <- "sim1.dta"
sim1 <- haven::read_dta(file)
attr(sim1, "spec") <- NULL
sim1 <- labelled::remove_labels(sim1, user_na_to_na = TRUE)
sim1 <- dplyr::filter(sim1, est %in% 1:3, model %in% c(1:3, 5:6)) %>%
  dplyr::mutate(
    model = factor(model, levels = c(1:3, 5:6), labels = c("JM", "ME (Total)", "ME (Cumulative)", "ME (No adjustment)", "GEE (IIVW)")),
    est = factor(est, levels = 1:3, labels = c("a0", "a1", "a2"))
  ) %>%
  dplyr::rename(par = est)

# Make multisimsum object
ms <- rsimsum::multisimsum(sim1, par = "par", true = c("a0" = 0, "a1" = 1, "a2" = 0.2), estvarname = "b", se = "se", methodvar = "model", by = "dgm", ref = "JM")

# DF with results
res <- rsimsum::get_data(summary(ms)) %>%
  dplyr::filter(stat %in% c("thetamean", "se2mean", "empse", "bias", "cover", "mse")) %>%
  dplyr::mutate(est = ifelse(stat == "se2mean", sqrt(est), est))

# Summary tables
tres <- res %>%
  dplyr::mutate(sig = dplyr::case_when(
    stat == "bias" & (lower > 0 | upper < 0) ~ "*",
    stat == "cover" & (lower > 0.95 | upper < 0.95) ~ "*",
    TRUE ~ ""
  )) %>%
  dplyr::mutate(est = paste0(formattable::comma(est, sigd))) %>%
  dplyr::mutate(est = ifelse(sig == "*", kableExtra::cell_spec(est, "latex", color = "red"), est)) %>%
  dplyr::mutate(value = ifelse(stat %in% c("thetamean", "se2mean"), paste0(est), paste0(est, " (", formattable::comma(mcse, sigd), ")"))) %>%
  dplyr::select(-est, -mcse, -lower, -upper, -sig) %>%
  dplyr::arrange(dgm) %>%
  tidyr::spread(key = stat, value = value) %>%
  dplyr::mutate(
    par = dplyr::case_when(
      par == "a0" ~ "\\(\\alpha_0\\)",
      par == "a1" ~ "\\(\\alpha_1\\)",
      par == "a2" ~ "\\(\\alpha_2\\)"
    ),
    model = factor(model, levels = c("JM", "ME (Total)", "ME (Cumulative)", "ME (No adjustment)", "GEE (IIVW)")),
  ) %>%
  dplyr::select(dgm, model, par, thetamean, se2mean, empse, bias, cover, mse) %>%
  dplyr::arrange(dgm, model, par)

for (i in as.character(1:10)) {
  sink(paste0("sim1_files/tab-dgm-", i, ".tex"))
  dplyr::filter(tres, dgm == i) %>%
    dplyr::select(-dgm, -model) %>%
    knitr::kable(
      format = "latex",
      booktabs = TRUE,
      align = rep("r", ncol(res) - 2),
      col.names = c("", "Estimate", "SE", "Empirical SE", "Bias", "Coverage", "MSE"),
      caption = paste0("\\label{tab:dgm-", i, "}Results for DGM: ", dgm[["label"]][dgm[["dgm"]] == i]),
      escape = FALSE
    ) %>%
    kableExtra::kable_styling(latex_options = c("scale_down", "hold_position")) %>%
    kableExtra::group_rows(paste0("Model: ", levels(tres$model)[1]), 1, 3) %>%
    kableExtra::group_rows(paste0("Model: ", levels(tres$model)[2]), 4, 6) %>%
    kableExtra::group_rows(paste0("Model: ", levels(tres$model)[3]), 7, 9) %>%
    kableExtra::group_rows(paste0("Model: ", levels(tres$model)[4]), 10, 12) %>%
    kableExtra::group_rows(paste0("Model: ", levels(tres$model)[5]), 13, 15) %>%
    print()
  sink()
  tmpf <- readr::read_lines(file = paste0("sim1_files/tab-dgm-", i, ".tex"))
  tmpf <- stringr::str_replace_all(string = tmpf, pattern = "\\\\label\\{tab:\\}", replacement = "")
  readr::write_lines(x = tmpf, path = paste0("sim1_files/tab-dgm-", i, ".tex"))
}

# Make a plot for bias, coverage, MSE
pres <- res %>%
  dplyr::left_join(dgm, "dgm") %>%
  dplyr::mutate(
    model = factor(model, levels = c("JM", "ME (Total)", "ME (Cumulative)", "ME (No adjustment)", "GEE (IIVW)"))
  ) %>%
  dplyr::mutate(sig = dplyr::case_when(
    stat == "bias" & (lower > 0 | upper < 0) ~ "*",
    stat == "cover" & (lower > 0.95 | upper < 0.95) ~ "*",
    TRUE ~ ""
  )) %>%
  dplyr::mutate(
    xx = dplyr::case_when(
      dgm == 1 ~ "JM (gamma = 0.00, lambda = 0.10)",
      dgm == 2 ~ "JM (gamma = 0.00, lambda = 0.30)",
      dgm == 3 ~ "JM (gamma = 0.00, lambda = 1.00)",
      dgm == 4 ~ "JM (gamma = 1.50, lambda = 0.10)",
      dgm == 5 ~ "JM (gamma = 1.50, lambda = 0.30)",
      dgm == 6 ~ "JM (gamma = 1.50, lambda = 1.00)",
      dgm == 7 ~ "Gamma not depending on treatment",
      dgm == 8 ~ "Gamma depending on treatment",
      dgm == 9 ~ "Gamma depending on treatment and previous Y",
      dgm == 10 ~ "JM (gamma = 3.00, lambda = 0.05) with regular visits"
    ),
    order = dplyr::case_when(
      dgm == 1 ~ 2,
      dgm == 2 ~ 3,
      dgm == 3 ~ 4,
      dgm == 4 ~ 7,
      dgm == 5 ~ 8,
      dgm == 6 ~ 9,
      dgm == 7 ~ 1,
      dgm == 8 ~ 5,
      dgm == 9 ~ 6,
      dgm == 10 ~ 10
    )
  ) %>%
  dplyr::mutate(facet = as.numeric(order > 4)) %>%
  dplyr::mutate(facet = factor(facet, levels = 0:1, labels = c("Non-Informative Observation Process", "Informative Observation Process"))) %>%
  dplyr::mutate(xx = stringr::str_wrap(xx, width = 12)) %>%
  dplyr::mutate(xx = factor(xx)) %>%
  dplyr::mutate(xx = forcats::fct_reorder(xx, order))

bias_a0 <- pres %>%
  dplyr::filter(stat == "bias" & par == "a0") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = c(0, 1), legend.justification = c(0, 1), legend.key = ggplot2::element_blank(), legend.background = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(Bias ~ of ~ alpha[0]), shape = "Model")
bias_a1 <- pres %>%
  dplyr::filter(stat == "bias" & par == "a1") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 5)) +
  labs(x = "", y = expression(Bias ~ of ~ alpha[1]), shape = "Model")
bias_a2 <- pres %>%
  dplyr::filter(stat == "bias" & par == "a2") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(Bias ~ of ~ alpha[2]), shape = "Model")
bias <- cowplot::plot_grid(bias_a0,
  bias_a1,
  bias_a2,
  ncol = 1,
  labels = LETTERS,
  align = "vh", axis = "ltrb"
)

ggplot2::ggsave(bias, filename = "sim1_files/bias.pdf", dpi = 1200, height = 10, width = 10 / 3 * 2)
rm(bias_a0, bias_a1, bias_a2, bias)

cover_a0 <- pres %>%
  dplyr::filter(stat == "cover" & par == "a0") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0.95, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = percent) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = c(0, 0), legend.justification = c(0, 0), legend.key = ggplot2::element_blank(), legend.background = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(Coverage ~ of ~ alpha[0]), shape = "Model")
cover_a1 <- pres %>%
  dplyr::filter(stat == "cover" & par == "a1") %>%
  ggplot2::ggplot(aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0.95, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = percent) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(Coverage ~ of ~ alpha[1]), shape = "Model")
cover_a2 <- pres %>%
  dplyr::filter(stat == "cover" & par == "a2") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0.95, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = percent) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(Coverage ~ of ~ alpha[2]), shape = "Model")
cover <- cowplot::plot_grid(cover_a0,
  cover_a1,
  cover_a2,
  ncol = 1,
  labels = LETTERS,
  align = "vh", axis = "ltrb"
)

ggplot2::ggsave(cover, filename = "sim1_files/cover.pdf", dpi = 1200, height = 10, width = 10 / 3 * 2)
rm(cover_a0, cover_a1, cover_a2, cover)

mse_a0 <- pres %>%
  dplyr::filter(stat == "mse" & par == "a0") %>%
  ggplot2::ggplot(aes(x = xx, y = est, shape = model, group = model)) +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = c(0, 1), legend.justification = c(0, 1), legend.key = ggplot2::element_blank(), legend.background = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(MSE ~ of ~ alpha[0]), shape = "Model")
mse_a1 <- pres %>%
  dplyr::filter(stat == "mse" & par == "a1") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model)) +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(MSE ~ of ~ alpha[1]), shape = "Model")
mse_a2 <- pres %>%
  dplyr::filter(stat == "mse" & par == "a2") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model)) +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(MSE ~ of ~ alpha[2]), shape = "Model")
mse <- cowplot::plot_grid(mse_a0,
  mse_a1,
  mse_a2,
  ncol = 1,
  labels = LETTERS,
  align = "vh", axis = "ltrb"
)

ggplot2::ggsave(mse, filename = "sim1_files/mse.pdf", dpi = 1200, height = 10, width = 10 / 3 * 2)
rm(mse_a0, mse_a1, mse_a2, mse)

### ---
# Plots with bias, coverage, mse for a given estimand
### ---
a0_bias <- pres %>%
  dplyr::filter(stat == "bias" & par == "a0") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = c(0, 1), legend.justification = c(0, 1), legend.key = ggplot2::element_blank(), legend.background = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(Bias ~ of ~ alpha[0]), shape = "Model")
a0_cover <- pres %>%
  dplyr::filter(stat == "cover" & par == "a0") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0.95, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(Coverage ~ of ~ alpha[0]), shape = "Model")
a0_mse <- pres %>%
  dplyr::filter(stat == "mse" & par == "a0") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(MSE ~ of ~ alpha[0]), shape = "Model")
a0 <- cowplot::plot_grid(a0_bias,
  a0_cover,
  a0_mse,
  ncol = 1,
  labels = LETTERS,
  align = "vh", axis = "ltrb"
)
ggplot2::ggsave(a0, filename = "sim1_files/a0.pdf", dpi = 1200, height = 10, width = 10 / 3 * 2)
rm(a0, a0_bias, a0_cover, a0_mse)

a1_bias <- pres %>%
  dplyr::filter(stat == "bias" & par == "a1") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(Bias ~ of ~ alpha[1]), shape = "Model")
a1_cover <- pres %>%
  dplyr::filter(stat == "cover" & par == "a1") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0.95, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(Coverage ~ of ~ alpha[1]), shape = "Model")
a1_mse <- pres %>%
  dplyr::filter(stat == "mse" & par == "a1") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = c(0, 1), legend.justification = c(0, 1), legend.key = ggplot2::element_blank(), legend.background = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(MSE ~ of ~ alpha[1]), shape = "Model")
a1 <- cowplot::plot_grid(a1_bias,
  a1_cover,
  a1_mse,
  ncol = 1,
  labels = LETTERS,
  align = "vh", axis = "ltrb"
)
ggplot2::ggsave(a1, filename = "sim1_files/a1.pdf", dpi = 1200, height = 10, width = 10 / 3 * 2)
rm(a1, a1_bias, a1_cover, a1_mse)

a2_bias <- pres %>%
  dplyr::filter(stat == "bias" & par == "a2") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = c(0, 1), legend.justification = c(0, 1), legend.key = ggplot2::element_blank(), legend.background = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(Bias ~ of ~ alpha[2]), shape = "Model")
a2_cover <- pres %>%
  dplyr::filter(stat == "cover" & par == "a2") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0.95, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(Coverage ~ of ~ alpha[2]), shape = "Model")
a2_mse <- pres %>%
  dplyr::filter(stat == "mse" & par == "a2") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(MSE ~ of ~ alpha[2]), shape = "Model")
a2 <- cowplot::plot_grid(a2_bias,
  a2_cover,
  a2_mse,
  ncol = 1,
  labels = LETTERS,
  align = "vh", axis = "ltrb"
)
ggplot2::ggsave(a2, filename = "sim1_files/a2.pdf", dpi = 1200, height = 10, width = 10 / 3 * 2)
rm(a2, a2_bias, a2_cover, a2_mse)

### ---
# Plots for the application
### ---

# ex-res.dta contains predicted values from the joint model
ex_data <- haven::read_dta(file = "sim1_files/ex-res.dta")
ex_data <- lapply(c("jm", "me", "mecnn", "mecumnn", "gee"), function(m) {
  tmp <- ex_data[, c("sex", "tY")]
  tmp[["est"]] <- ex_data[[paste0("egfr_", m)]]
  tmp[["lower"]] <- ex_data[[paste0("egfr_", m, "_lci")]]
  tmp[["upper"]] <- ex_data[[paste0("egfr_", m, "_uci")]]
  tmp[["model"]] <- m
  tmp
})
ex_data <- dplyr::bind_rows(ex_data) %>%
  dplyr::mutate(sex = factor(sex, levels = 1:2, labels = c("Males", "Females"))) %>%
  dplyr::mutate(model = factor(model, levels = c("jm", "mecnn", "mecumnn", "me", "gee"), labels = c("JM", "ME (Total)", "ME (Cumulative)", "ME (No adjustment)", "GEE (IIVW)"))) %>%
  dplyr::mutate(tY2 = round(tY, digits = 2)) %>%
  dplyr::group_by(model, sex, tY) %>%
  dplyr::summarise(
    est = mean(est),
    lower = mean(lower),
    upper = mean(upper)
  ) %>%
  dplyr::ungroup()

ex_p1 <- ex_data %>%
  ggplot2::ggplot(ggplot2::aes(x = tY, y = est, ymin = lower, ymax = upper, fill = model, lty = sex)) +
  ggplot2::geom_ribbon(alpha = 0.1) +
  ggplot2::geom_line(aes(color = model)) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  ggplot2::scale_linetype_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom", legend.box = "vertical") +
  ggplot2::scale_x_continuous(breaks = seq(-5, 0, by = 1)) +
  ggplot2::labs(x = "Follow-up time (Years)", y = "Predicted eGFR", colour = "Model:", fill = "Model:", lty = "Sex:")
ggplot2::ggsave(ex_p1, filename = "sim1_files/ex_p1.pdf", dpi = 1200, height = 6.5 * 2 / 3, width = 6.5)

# Forest plot for model coefficients
modcoefs <- readxl::read_excel(path = "sim1_ex-model-results.xlsx") %>%
  dplyr::mutate(model = factor(model, levels = c("jm", "mixed_cnn", "mixed_cumnn", "mixed", "gee"), labels = c("JM", "ME (Total)", "ME (Cumulative)", "ME (No adjustment)", "GEE (IIVW)"))) %>%
  dplyr::filter(grepl("^egfr_", coef) & coef != "egfr_gamma") %>%
  dplyr::mutate(coef = factor(coef, levels = c("egfr_cons", "egfr_sex", "egfr_tY", "egfr_sex#tY"), labels = c("Intercept", "Gender", "Time", "Gender x Time")))

ex_p2 <- modcoefs %>%
  dplyr::group_by(coef) %>%
  dplyr::mutate(xtext = min(value - 6 * se)) %>%
  dplyr::ungroup() %>%
  ggplot2::ggplot(ggplot2::aes(x = value, y = model)) +
  ggplot2::geom_errorbarh(aes(xmin = value - 1.96 * se, xmax = value + 1.96 * se), height = 0.5) +
  ggplot2::geom_point() +
  ggplot2::geom_text(aes(label = formattable::comma(value, digits = 4), x = xtext), hjust = 0, size = 3, fontface = "italic") +
  ggplot2::scale_x_continuous(labels = formattable::comma) +
  ggplot2::facet_wrap(~coef, scales = "free_x") +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Estimated regression coefficient", y = "")
ggplot2::ggsave(ex_p2, filename = "sim1_files/ex_p2.pdf", dpi = 1200, height = 5.0, width = 6.5)

# Clean up
rm(ex_data, ex_p1, modcoefs, ex_p2)

### ---
#   Supplementary material
### ---

# Bias and coverage of variances (residual, random effects)
sim1_var <- haven::read_dta(file)
attr(sim1_var, "spec") <- NULL
sim1_var <- labelled::remove_labels(sim1_var, user_na_to_na = TRUE)
sim1_var <- dplyr::filter(sim1_var, est %in% 4:6 & model %in% c(1:3, 5:6)) %>%
  dplyr::mutate(
    model = factor(model, levels = c(1:3, 5:6), labels = c("JM", "ME (Total)", "ME (Cumulative)", "ME (No adjustment)", "GEE (IIVW)")),
    est = factor(est, levels = 4:6, labels = c("V(u)", "V(v)", "V(e)"))
  ) %>%
  dplyr::rename(par = est)
ms_var <- rsimsum::multisimsum(sim1_var, par = "par", true = c("V(u)" = 1, "V(v)" = 0.5, "V(e)" = 1), estvarname = "b", se = "se", methodvar = "model", by = "dgm", ref = "JM")
res_var <- rsimsum::get_data(summary(ms_var)) %>%
  dplyr::filter(stat %in% c("thetamean", "se2mean", "empse", "bias", "cover", "mse")) %>%
  dplyr::mutate(est = ifelse(stat == "se2mean", sqrt(est), est))
pres_var <- rsimsum::get_data(summary(ms_var)) %>%
  dplyr::filter(stat %in% c("bias", "cover", "mse")) %>%
  dplyr::left_join(dgm, "dgm") %>%
  dplyr::mutate(
    model = factor(model, levels = c("JM", "ME (Total)", "ME (Cumulative)", "ME (No adjustment)", "GEE (IIVW)"))
  ) %>%
  dplyr::mutate(sig = dplyr::case_when(
    stat == "bias" & (lower > 0 | upper < 0) ~ "*",
    stat == "cover" & (lower > 0.95 | upper < 0.95) ~ "*",
    TRUE ~ ""
  )) %>%
  dplyr::mutate(
    xx = dplyr::case_when(
      dgm == 1 ~ "JM (gamma = 0.00, lambda = 0.10)",
      dgm == 2 ~ "JM (gamma = 0.00, lambda = 0.30)",
      dgm == 3 ~ "JM (gamma = 0.00, lambda = 1.00)",
      dgm == 4 ~ "JM (gamma = 1.50, lambda = 0.10)",
      dgm == 5 ~ "JM (gamma = 1.50, lambda = 0.30)",
      dgm == 6 ~ "JM (gamma = 1.50, lambda = 1.00)",
      dgm == 7 ~ "Gamma not depending on treatment",
      dgm == 8 ~ "Gamma depending on treatment",
      dgm == 9 ~ "Gamma depending on treatment and previous Y",
      dgm == 10 ~ "JM (gamma = 3.00, lambda = 0.05) with regular visits"
    ),
    order = dplyr::case_when(
      dgm == 1 ~ 2,
      dgm == 2 ~ 3,
      dgm == 3 ~ 4,
      dgm == 4 ~ 7,
      dgm == 5 ~ 8,
      dgm == 6 ~ 9,
      dgm == 7 ~ 1,
      dgm == 8 ~ 5,
      dgm == 9 ~ 6,
      dgm == 10 ~ 10
    )
  ) %>%
  dplyr::mutate(facet = as.numeric(order > 4)) %>%
  dplyr::mutate(facet = factor(facet, levels = 0:1, labels = c("Non-Informative Observation Process", "Informative Observation Process"))) %>%
  dplyr::mutate(xx = stringr::str_wrap(xx, width = 15)) %>%
  dplyr::mutate(xx = factor(xx)) %>%
  dplyr::mutate(xx = forcats::fct_reorder(xx, order)) %>%
  dplyr::mutate(est = ifelse(dgm > 6 & par == "V(u)", NA, est))

bias_vu <- pres_var %>%
  dplyr::filter(stat == "bias" & par == "V(u)") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = c(0, 1), legend.justification = c(0, 1), legend.key = ggplot2::element_blank(), legend.background = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = "Bias of V(u)", shape = "Model")
bias_vv <- pres_var %>%
  dplyr::filter(stat == "bias" & par == "V(v)") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = "Bias of V(v)", shape = "Model")
bias_ve <- pres_var %>%
  dplyr::filter(stat == "bias" & par == "V(e)") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(Bias ~ of ~ V(epsilon)), shape = "Model")
bias_var <- cowplot::plot_grid(bias_vu,
  bias_vv,
  bias_ve,
  ncol = 1,
  labels = LETTERS,
  align = "vh", axis = "ltrb"
)
ggplot2::ggsave(bias_var, filename = "sim1_files/bias_var.pdf", dpi = 1200, height = 10, width = 10 / 3 * 2)
rm(bias_vu, bias_vv, bias_ve, bias_var)

cover_vu <- pres_var %>%
  dplyr::filter(stat == "cover" & par == "V(u)") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0.95, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = percent) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = c(0, 1), legend.justification = c(0, 1), legend.key = ggplot2::element_blank(), legend.background = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = "Coverage of V(u)", shape = "Model")
cover_vv <- pres_var %>%
  dplyr::filter(stat == "cover" & par == "V(v)") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0.95, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = percent) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = "Coverage of V(v)", shape = "Model")
cover_ve <- pres_var %>%
  dplyr::filter(stat == "cover" & par == "V(e)") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0.95, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = percent) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(Coverage ~ of ~ V(epsilon)), shape = "Model")
cover_var <- cowplot::plot_grid(cover_vu,
  cover_vv,
  cover_ve,
  ncol = 1,
  labels = LETTERS,
  align = "vh", axis = "ltrb"
)
ggplot2::ggsave(cover_var, filename = "sim1_files/cover_var.pdf", dpi = 1200, height = 10, width = 10 / 3 * 2)
rm(cover_vu, cover_vv, cover_ve, cover_var)

mse_vu <- pres_var %>%
  dplyr::filter(stat == "mse" & par == "V(u)") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = c(1, 1), legend.justification = c(1, 1), legend.key = ggplot2::element_blank(), legend.background = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = "MSE of V(u)", shape = "Model")
mse_vv <- pres_var %>%
  dplyr::filter(stat == "mse" & par == "V(v)") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = "MSE of V(v)", shape = "Model")
mse_ve <- pres_var %>%
  dplyr::filter(stat == "mse" & par == "V(e)") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 5)) +
  ggplot2::labs(x = "", y = expression(MSE ~ of ~ V(epsilon)), shape = "Model")
mse_var <- cowplot::plot_grid(mse_vu,
  mse_vv,
  mse_ve,
  ncol = 1,
  labels = LETTERS,
  align = "vh", axis = "ltrb"
)
ggplot2::ggsave(mse_var, filename = "sim1_files/mse_var.pdf", dpi = 1200, height = 10, width = 10 / 3 * 2)
rm(mse_vu, mse_vv, mse_ve, mse_var)
rm(sim1_var, pres_var, ms_var)

# Tables for variances
tres_var <- res_var %>%
  dplyr::mutate(sig = dplyr::case_when(
    stat == "bias" & (lower > 0 | upper < 0) ~ "*",
    stat == "cover" & (lower > 0.95 | upper < 0.95) ~ "*",
    TRUE ~ ""
  )) %>%
  dplyr::mutate(est = ifelse(dgm > 6 & par == "V(u)", NA, est)) %>%
  dplyr::mutate(est = paste0(formattable::comma(est, sigd))) %>%
  dplyr::mutate(est = ifelse(sig == "*", kableExtra::cell_spec(est, "latex", color = "red"), est)) %>%
  dplyr::mutate(value = ifelse(stat %in% c("thetamean", "se2mean"), paste0(est), paste0(est, " (", formattable::comma(mcse, sigd), ")"))) %>%
  dplyr::select(-est, -mcse, -lower, -upper, -sig) %>%
  dplyr::arrange(dgm) %>%
  tidyr::spread(key = stat, value = value) %>%
  dplyr::mutate(
    model = factor(model, levels = c("JM", "ME (Total)", "ME (Cumulative)", "ME (No adjustment)", "GEE (IIVW)"))
  ) %>%
  dplyr::select(dgm, model, par, thetamean, se2mean, empse, bias, cover, mse) %>%
  dplyr::mutate(par = ifelse(par == "V(e)", "V(\\(\\epsilon\\))", par)) %>%
  dplyr::arrange(dgm, model, par) %>%
  dplyr::mutate(
    thetamean = ifelse((par == "V(u)" & model != "JM") | model == "GEE (IIVW)" | (dgm > 6 & par == "V(u)"), "---", thetamean),
    se2mean = ifelse((par == "V(u)" & model != "JM") | model == "GEE (IIVW)" | (dgm > 6 & par == "V(u)"), "---", se2mean),
    empse = ifelse((par == "V(u)" & model != "JM") | model == "GEE (IIVW)" | (dgm > 6 & par == "V(u)"), "---", empse),
    bias = ifelse((par == "V(u)" & model != "JM") | model == "GEE (IIVW)" | (dgm > 6 & par == "V(u)"), "---", bias),
    cover = ifelse((par == "V(u)" & model != "JM") | model == "GEE (IIVW)" | (dgm > 6 & par == "V(u)"), "---", cover),
    mse = ifelse((par == "V(u)" & model != "JM") | model == "GEE (IIVW)" | (dgm > 6 & par == "V(u)"), "---", mse)
  )

for (i in as.character(1:10)) {
  sink(paste0("sim1_files/tab-var-dgm-", i, ".tex"))
  dplyr::filter(tres_var, dgm == i) %>%
    dplyr::select(-dgm, -model) %>%
    knitr::kable(
      format = "latex",
      booktabs = TRUE,
      align = rep("r", ncol(res_var) - 2),
      col.names = c("", "Estimate", "SE", "Empirical SE", "Bias", "Coverage", "MSE"),
      caption = paste0("\\label{tab:var-dgm-", i, "}Results for DGM: ", dgm[["label"]][dgm[["dgm"]] == i]),
      escape = FALSE
    ) %>%
    kableExtra::kable_styling(latex_options = c("scale_down", "hold_position")) %>%
    kableExtra::group_rows(paste0("Model: ", levels(tres_var$model)[1]), 1, 3) %>%
    kableExtra::group_rows(paste0("Model: ", levels(tres_var$model)[2]), 4, 6) %>%
    kableExtra::group_rows(paste0("Model: ", levels(tres_var$model)[3]), 7, 9) %>%
    kableExtra::group_rows(paste0("Model: ", levels(tres_var$model)[4]), 10, 12) %>%
    kableExtra::group_rows(paste0("Model: ", levels(tres_var$model)[5]), 13, 15) %>%
    print()
  sink()
  tmpf <- readr::read_lines(file = paste0("sim1_files/tab-var-dgm-", i, ".tex"))
  tmpf <- stringr::str_replace_all(string = tmpf, pattern = "\\\\label\\{tab:\\}", replacement = "")
  readr::write_lines(x = tmpf, path = paste0("sim1_files/tab-var-dgm-", i, ".tex"))
}

# Summary statistics for gamma
sim1_gamma <- haven::read_dta(file)
attr(sim1_gamma, "spec") <- NULL
sim1_gamma <- labelled::remove_labels(sim1_gamma, user_na_to_na = TRUE)
sim1_gamma <- dplyr::filter(sim1_gamma, est == 7 & model %in% c(1:3, 5:6) & dgm %in% c(1:6, 10)) %>%
  dplyr::mutate(
    model = factor(model, levels = c(1:3, 5:6), labels = c("JM", "ME (Total)", "ME (Cumulative)", "ME (No adjustment)", "GEE (IIVW)")),
    est = ifelse(dgm %in% 1:3, "gamma.0", ifelse(dgm == 10, "gamma.3.0", "gamma.1.5"))
  ) %>%
  dplyr::rename(par = est)
ms_gamma <- rsimsum::multisimsum(sim1_gamma, par = "par", true = c("gamma.0" = 0, "gamma.1.5" = 1.5, "gamma.3.0" = 3), estvarname = "b", se = "se", methodvar = "model", by = "dgm", ref = "JM")
pres_gamma <- rsimsum::get_data(summary(ms_gamma)) %>%
  dplyr::filter(stat %in% c("bias", "cover", "mse")) %>%
  dplyr::left_join(dgm, "dgm") %>%
  dplyr::mutate(
    model = factor(model, levels = c("JM", "ME (Total)", "ME (Cumulative)", "ME (No adjustment)", "GEE (IIVW)"))
  ) %>%
  dplyr::mutate(sig = dplyr::case_when(
    stat == "bias" & (lower > 0 | upper < 0) ~ "*",
    stat == "cover" & (lower > 0.95 | upper < 0.95) ~ "*",
    TRUE ~ ""
  )) %>%
  dplyr::mutate(
    xx = dplyr::case_when(
      dgm == 1 ~ "JM (gamma = 0.00, lambda = 0.10)",
      dgm == 2 ~ "JM (gamma = 0.00, lambda = 0.30)",
      dgm == 3 ~ "JM (gamma = 0.00, lambda = 1.00)",
      dgm == 4 ~ "JM (gamma = 1.50, lambda = 0.10)",
      dgm == 5 ~ "JM (gamma = 1.50, lambda = 0.30)",
      dgm == 6 ~ "JM (gamma = 1.50, lambda = 1.00)",
      dgm == 10 ~ "JM (gamma = 3.00, lambda = 0.05) with regular visits"
    )
  ) %>%
  dplyr::mutate(facet = as.numeric(as.numeric(dgm) > 3)) %>%
  dplyr::mutate(facet = factor(facet, levels = 0:1, labels = c("Non-Informative Observation Process", "Informative Observation Process"))) %>%
  dplyr::mutate(xx = stringr::str_wrap(xx, width = 15)) %>%
  dplyr::mutate(xx = factor(xx)) %>%
  dplyr::mutate(xx = forcats::fct_reorder(xx, as.numeric(dgm))) %>%
  dplyr::mutate(par = dplyr::case_when(
    par == "gamma.0" ~ "gamma",
    par == "gamma.1.5" ~ "gamma",
    par == "gamma.3.0" ~ "gamma",
    TRUE ~ par
  ))

res_gamma <- rsimsum::get_data(summary(ms_gamma)) %>%
  dplyr::filter(stat %in% c("thetamean", "se2mean", "empse", "bias", "cover", "mse")) %>%
  dplyr::mutate(est = ifelse(stat == "se2mean", sqrt(est), est)) %>%
  dplyr::filter(model == "JM") %>%
  dplyr::select(-model, -par) %>%
  dplyr::mutate(sig = dplyr::case_when(
    stat == "bias" & (lower > 0 | upper < 0) ~ "*",
    stat == "cover" & (lower > 0.95 | upper < 0.95) ~ "*",
    TRUE ~ ""
  ))

tres_gamma <- res_gamma %>%
  dplyr::mutate(est = paste0(formattable::comma(est, sigd))) %>%
  dplyr::mutate(est = ifelse(sig == "*", kableExtra::cell_spec(est, "latex", color = "red"), est)) %>%
  dplyr::mutate(value = ifelse(stat %in% c("thetamean", "se2mean"), paste0(est), paste0(est, " (", formattable::comma(mcse, sigd), ")"))) %>%
  dplyr::select(-est, -mcse, -lower, -upper, -sig) %>%
  dplyr::arrange(dgm) %>%
  tidyr::spread(key = stat, value = value) %>%
  dplyr::left_join(dgm, by = "dgm") %>%
  dplyr::select(gamma, lambda, thetamean, se2mean, empse, bias, cover, mse) %>%
  dplyr::arrange(gamma, lambda)

sink("sim1_files/tab-gamma.tex")
knitr::kable(tres_gamma,
  format = "latex",
  booktabs = TRUE,
  align = rep("r", ncol(tres_gamma)),
  col.names = c("True \\(\\gamma\\)", "True \\(\\lambda\\)", "Estimate", "SE", "Empirical SE", "Bias", "Coverage", "MSE"),
  caption = "\\label{tab:gamma}Results for the association parameter \\(\\gamma\\)",
  escape = FALSE
) %>%
  kableExtra::kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  print()
sink()
tmpf <- readr::read_lines(file = "sim1_files/tab-gamma.tex")
tmpf <- stringr::str_replace_all(string = tmpf, pattern = "\\\\label\\{tab:\\}", replacement = "")
readr::write_lines(x = tmpf, path = "sim1_files/tab-gamma.tex")

bias_gamma <- pres_gamma %>%
  dplyr::filter(stat == "bias" & par == "gamma") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = c(0, 0), legend.justification = c(0, 0), legend.key = ggplot2::element_blank(), legend.background = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 7)) +
  ggplot2::labs(x = "", y = expression(Bias ~ of ~ association ~ parameter ~ gamma), shape = "Model")

cover_gamma <- pres_gamma %>%
  dplyr::filter(stat == "cover" & par == "gamma") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model, color = sig)) +
  ggplot2::geom_hline(yintercept = 0.95, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggthemes::scale_color_colorblind(guide = "none") +
  ggplot2::scale_y_continuous(labels = percent) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 7)) +
  ggplot2::labs(x = "", y = expression(Coverage ~ of ~ association ~ parameter ~ gamma), shape = "Model")

mse_gamma <- pres_gamma %>%
  dplyr::filter(stat == "mse" & par == "gamma") %>%
  ggplot2::ggplot(ggplot2::aes(x = xx, y = est, shape = model, group = model)) +
  ggplot2::geom_hline(yintercept = 0.95, color = "grey50", linetype = "dotted") +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
  ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_x") +
  ggplot2::scale_y_continuous(labels = comma.3sigd) +
  ggplot2::scale_shape_discrete() +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.margin = ggplot2::unit(c(1, 0, 1, 0), "pt"), panel.spacing = ggplot2::unit(0, "pt"), legend.position = "none", axis.text.x = ggplot2::element_text(size = 7)) +
  ggplot2::labs(x = "", y = expression(MSE ~ of ~ association ~ parameter ~ gamma), shape = "Model")

summary_gamma <- cowplot::plot_grid(bias_gamma,
  cover_gamma,
  mse_gamma,
  ncol = 1,
  labels = LETTERS,
  align = "vh", axis = "ltrb"
)

ggplot2::ggsave(summary_gamma, filename = "sim1_files/summary_gamma.pdf", dpi = 1200, height = 10, width = 10 / 3 * 2)

rm(summary_gamma, mse_gamma, cover_gamma, bias_gamma, tmpf, tres_gamma, res_gamma, pres_gamma, ms_gamma, sim1_gamma)
gc()

### ---
# Convergence summaries
### ---
sim1_convergence <- read_dta(file)
attr(sim1_convergence, "spec") <- NULL
sim1_convergence <- labelled::remove_labels(sim1_convergence, user_na_to_na = TRUE)
sim1_convergence <- dplyr::filter(sim1_convergence, model %in% c(1:3, 5:6)) %>%
  dplyr::distinct(dgm, model, i, .keep_all = TRUE) %>%
  dplyr::group_by(dgm, model) %>%
  dplyr::summarise(
    complete_n = sum(!is.na(b) & !is.na(se)),
    complete_p = mean(!is.na(b) & !is.na(se))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(value = paste0(formattable::comma(complete_n, digits = 0), " (", formattable::comma(complete_p * 100, digits = 0), "\\%)")) %>%
  dplyr::select(-dplyr::starts_with("complete_")) %>%
  dplyr::mutate(
    model = factor(model, levels = c(1:3, 5:6), labels = c("JM", "ME (Total)", "ME (Cumulative)", "ME (No adjustment)", "GEE (IIVW)"))
  ) %>%
  tidyr::spread(key = model, value = value) %>%
  dplyr::mutate(dgm = as.character(dgm)) %>%
  dplyr::left_join(dgm, by = "dgm") %>%
  dplyr::select(11, 2:6)

sink("sim1_files/tab-convergence.tex")
sim1_convergence %>%
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    align = rep("r", ncol(sim1_convergence)),
    col.names = c("Data-Generating Mechanism", paste("Model", LETTERS[1:5])),
    caption = "\\label{tab:convergence}Number and percentage of models converging under each data-generating mechanism. Model A is the joint model, model B is the mixed model adjusting for the total number of measurements, model C is the mixed model adjusting for the cumulative number of measurements, model D is the mixed model with no further adjustment, and model E is the marginal model fitted using GEE and IIVW.",
    escape = FALSE,
    linesep = ""
  ) %>%
  kableExtra::kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  print()
sink()
tmpf <- readr::read_lines(file = "sim1_files/tab-convergence.tex")
tmpf <- stringr::str_replace_all(string = tmpf, pattern = "\\\\label\\{tab:\\}", replacement = "")
readr::write_lines(x = tmpf, path = "sim1_files/tab-convergence.tex")
