# Script version of figures/03_drift.ipynb (run from the figures/ folder: Rscript scripts/03_drift.R)
source("fig_setup.R")
options(repr.plot.width = 12, repr.plot.height = 5, repr.plot.res = 120)
dc <- read_csv(file.path(DRIFT, "drift_correction_final.csv"), show_col_types = FALSE)
yr <- read_csv(file.path(DRIFT, "ground_drift_500", "ground_drift_within_tile_by_year.csv"), show_col_types = FALSE)
order_bins <- c("0-10%", "10-30%", "30-50%", "50-70%", "70-80%", "80-90%", "90-100%")

fa <- ggplot(dc, aes(x = cover_mid, y = excess_rate)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_ribbon(aes(ymin = excess_lo, ymax = excess_hi), fill = "#92c5de", alpha = 0.5) +
  geom_line(colour = "#2166ac", linewidth = 0.8) + geom_point(colour = "#2166ac", size = 2) +
  labs(x = "Landsat tree cover (%)", y = "Excess ground drift (m / year)",
       title = "a  Ground-finding drift beyond open ground (biases RH)") + theme_fig()
yb <- yr %>% filter(cover_bin %in% c("0-10%", "50-70%", "90-100%")) %>%
  mutate(cover_bin = factor(cover_bin, levels = order_bins))
fb <- ggplot(yb, aes(x = yr, y = estimate, colour = cover_bin)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  scale_colour_manual(values = c("#bdbdbd", "#6baed6", "#08306b")) +
  labs(x = NULL, y = "Ground residual vs TanDEM-X,\nrelative to tile T1 mean (m)", colour = "Tree cover",
       title = "b  Year effect within tiles") + theme_fig() + theme(legend.position = "top")
fs <- fa + fb
fs
save_fig(fs, "figS1_ground_drift", 12, 5)