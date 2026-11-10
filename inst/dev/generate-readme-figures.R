# generate-readme-figures.R
# Run this script from the package root to regenerate man/figures/ PNGs.
# Requires: TemporalHazard (installed), ggplot2, survival
#
# Both figures come from the CABGKUL dataset shipped with the package. The
# layout follows the "three phases add up, one survival curve" pair: hazard in
# percent per year against years, phases labelled on the curves, and the
# survival curve shaded where the early and late phases dominate.

library(TemporalHazard)
library(ggplot2)
library(survival)

# -- Load the CABGKUL dataset (lazy-loaded with the package) -----------------
data(cabgkul)

# -- Fit a 3-phase multiphase model on CABG death (time in months) -----------
fit <- hazard(
  Surv(int_dead, dead) ~ 1,
  data   = cabgkul,
  dist   = "multiphase",
  phases = list(
    early    = hzr_phase("cdf", t_half = 0.2, nu = 1, m = 1,
                          fixed = "shapes"),
    constant = hzr_phase("constant"),
    late     = hzr_phase("g3",  tau = 1, gamma = 3, alpha = 1, eta = 1,
                          fixed = "shapes")
  ),
  fit     = TRUE,
  control = list(n_starts = 10, maxit = 2000, reltol = 1e-10)
)

# -- Time grid: dense near zero, where the early phase changes fastest -------
max_years <- 15
t_grid <- c(seq(0.005, 2, length.out = 400),
            seq(2.05, max_years * 12, length.out = 1600))

# -- Per-phase hazard from the decomposed cumulative hazard ------------------
# predict(decompose = TRUE) returns cumulative hazard by phase; the hazard is
# its slope, taken between grid points and placed at the midpoints.
decomp_H <- predict(fit, newdata = data.frame(time = t_grid),
                    type = "cumulative_hazard", decompose = TRUE)
t_mid  <- (head(t_grid, -1) + tail(t_grid, -1)) / 2
phases <- c("early", "constant", "late", "total")
h_mid  <- sapply(phases, function(p) diff(decomp_H[[p]]) / diff(t_grid))

# The figure claims the phases add to the total. Check that they do, and that
# the slope matches the model's own hazard, before drawing anything.
h_model <- predict(fit, newdata = data.frame(time = t_mid), type = "hazard")
stopifnot(
  max(abs(rowSums(h_mid[, 1:3]) / h_mid[, "total"] - 1)) < 1e-8,
  max(abs(h_mid[, "total"] / h_model - 1)) < 1e-2,
  diff(range(h_mid[, "early"])) > 0,
  diff(range(h_mid[, "late"])) > 0
)

# Months to years; hazard per month to percent per year.
years   <- t_mid / 12
pct_yr  <- h_mid * 12 * 100

h_long <- data.frame(
  years  = rep(years, 4),
  hazard = as.vector(pct_yr),
  Phase  = factor(rep(c("Early", "Constant", "Late", "Total"),
                      each = length(years)),
                  levels = c("Early", "Constant", "Late", "Total"))
)

phase_colours <- c(
  "Early"    = "#0072B2",
  "Constant" = "#009E73",
  "Late"     = "#D55E00",
  "Total"    = "#111111"
)

# -- Figure 1: three phases add up -------------------------------------------
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)

y_top <- 10
at_year <- function(y, phase) pct_yr[which.min(abs(years - y)), phase]

p1 <- ggplot(h_long, aes(years, hazard, colour = Phase, linetype = Phase)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = phase_colours, guide = "none") +
  scale_linetype_manual(values = c(Early = "solid", Constant = "solid",
                                   Late = "solid", Total = "22"),
                        guide = "none") +
  annotate("text", x = 0.6, y = 8.6, label = "Early", hjust = 0,
           colour = phase_colours[["Early"]], fontface = "bold", size = 4.5) +
  annotate("text", x = 1, y = at_year(1, "constant") - 0.45,
           label = "Constant", hjust = 0,
           colour = phase_colours[["Constant"]], fontface = "bold",
           size = 4.5) +
  annotate("text", x = 13.2, y = at_year(13, "late") - 0.9, label = "Late",
           hjust = 0, colour = phase_colours[["Late"]], fontface = "bold",
           size = 4.5) +
  annotate("text", x = 7, y = at_year(7, "total") + 1.1,
           label = "Total = sum", hjust = 0.5,
           colour = phase_colours[["Total"]], fontface = "bold", size = 4.5) +
  coord_cartesian(xlim = c(0, max_years), ylim = c(0, y_top),
                  expand = FALSE, clip = "on") +
  scale_x_continuous(breaks = seq(0, max_years, 5)) +
  labs(x = "Years after CABG",
       y = "Risk of death per year (%)",
       title = "Three phases add up",
       subtitle = "Multiphase hazard fitted to CABGKUL (n = 5,880)") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.margin = margin(8, 14, 8, 8))

ggsave("man/figures/readme-hazard-phases.png", p1,
       width = 7, height = 4.5, dpi = 150, bg = "white")

# -- Figure 2: one survival curve, with Kaplan-Meier -------------------------
# Shade where a phase carries more risk than the constant phase: the early
# phase until it drops below it, the late phase once it rises above it.
early_end  <- years[which(pct_yr[, "early"] < pct_yr[, "constant"])[1]]
late_start <- years[which(pct_yr[, "late"] > pct_yr[, "constant"])[1]]
stopifnot(is.finite(early_end), is.finite(late_start),
          early_end < late_start)

surv_df <- data.frame(
  years    = t_grid / 12,
  survival = predict(fit, newdata = data.frame(time = t_grid),
                     type = "survival") * 100
)
km    <- survfit(Surv(int_dead, dead) ~ 1, data = cabgkul)
km_df <- data.frame(years = c(0, km$time / 12), survival = c(100, km$surv * 100))
km_df <- km_df[km_df$years <= max_years, ]

regions <- data.frame(
  xmin  = c(0, late_start),
  xmax  = c(early_end, max_years),
  Phase = c("Early", "Late")
)

p2 <- ggplot() +
  geom_rect(data = regions,
            aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 100, fill = Phase),
            alpha = 0.12) +
  scale_fill_manual(values = phase_colours, guide = "none") +
  # Kaplan-Meier goes on top: it sits on the model curve for the first
  # decade and would otherwise be hidden under it.
  geom_line(data = surv_df,
            aes(years, survival, colour = "Multiphase model"),
            linewidth = 1.6) +
  geom_step(data = km_df, aes(years, survival, colour = "Kaplan-Meier"),
            linewidth = 0.5) +
  scale_colour_manual(values = c("Multiphase model" = "#111111",
                                 "Kaplan-Meier" = "#E69F00"),
                      breaks = c("Multiphase model", "Kaplan-Meier")) +
  annotate("text", x = (late_start + max_years) / 2, y = 15,
           label = "late\nacceleration", colour = phase_colours[["Late"]],
           fontface = "bold", size = 4, lineheight = 0.9) +
  annotate("text", x = (early_end + late_start) / 2, y = 15,
           label = "steady decline", colour = phase_colours[["Constant"]],
           fontface = "bold", size = 4) +
  coord_cartesian(xlim = c(0, max_years), ylim = c(0, 100), expand = FALSE) +
  scale_x_continuous(breaks = seq(0, max_years, 5)) +
  labs(x = "Years after CABG",
       y = "Freedom from death (%)",
       title = "One survival curve",
       subtitle = "Early phase shaded at left, too brief to label",
       colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.margin = margin(8, 14, 8, 8))

ggsave("man/figures/readme-survival.png", p2,
       width = 7, height = 4.5, dpi = 150, bg = "white")

message(sprintf("Early phase dominates to %.2f years; late from %.1f years.",
                early_end, late_start))
message("Figures written to man/figures/")
