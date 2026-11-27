# Script version of figures/02_summaries.ipynb (run from the figures/ folder: Rscript scripts/02_summaries.R)
source("fig_setup.R")
options(repr.plot.width = 12, repr.plot.height = 8, repr.plot.res = 120)
d <- load_inputs()
cont <- d$pa %>% distinct(iso3, continent)
cm <- d$country %>% filter(metric == MAIN) %>% left_join(d$boot_country %>% filter(metric == MAIN), by = c("metric", "iso3")) %>%
  left_join(cont, by = "iso3")

top <- cm %>% filter(pct_pas_qualifying >= MIN_PCT_COUNTERFACTUAL) %>% slice_max(abs(avoided_loss_Tg), n = 20)
f4a <- ggplot(top, aes(x = avoided_loss_Tg, y = fct_reorder(iso3, avoided_loss_Tg), colour = continent)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_errorbarh(aes(xmin = avoided_loss_Tg_lo, xmax = avoided_loss_Tg_hi), height = 0, linewidth = 0.6) +
  geom_point(size = 2.6) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Avoided loss: DiD x PA area (Tg AGB)", y = NULL, colour = NULL,
       title = "a  Biomass retained by protection, T1 to T2") + theme_fig()
top2 <- cm %>% slice_max(abs(stock_change_Tg), n = 20)
f4b <- ggplot(top2, aes(x = stock_change_Tg, y = fct_reorder(iso3, stock_change_Tg), colour = continent)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_errorbarh(aes(xmin = stock_change_Tg_lo, xmax = stock_change_Tg_hi), height = 0, linewidth = 0.6) +
  geom_point(size = 2.6) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Change in PA biomass stock (Tg AGB)", y = NULL, colour = NULL,
       title = "b  Change in protected-area stocks (drift-corrected)") + theme_fig()
f4 <- f4a + f4b + plot_layout(guides = "collect") & theme(legend.position = "bottom")
f4
save_fig(f4, "fig4_top_countries", 12, 7)

cb <- d$cb %>% filter(metric == MAIN, n_qualifying >= 5) %>%
  select(continent, biome_name, n_qualifying, change_pa_mean, change_ctl_mean, did_mean) %>%
  pivot_longer(c(change_pa_mean, change_ctl_mean), names_to = "group", values_to = "change") %>%
  mutate(group = recode(group, change_pa_mean = "Protected", change_ctl_mean = "Counterfactual"))
f5 <- ggplot(cb, aes(x = change, y = biome_name)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_line(aes(group = biome_name), colour = "grey70", linewidth = 1) +
  geom_point(aes(colour = group), size = 2.6) +
  facet_wrap(~ continent, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = c(Protected = "#01665e", Counterfactual = "#bf812d")) +
  labs(x = "Mean AGBD change, T1 to T2 (Mg/ha, drift-corrected)", y = NULL, colour = NULL,
       title = "Biomass change in protected areas and their counterfactuals") +
  theme_fig() + theme(legend.position = "top")
f5
save_fig(f5, "fig5_continent_biome_change", 12, 9)

g <- d$global %>% filter(metric %in% names(METRIC_LABEL)) %>%
  left_join(d$boot_global, by = "metric") %>%
  mutate(label = factor(METRIC_LABEL[metric], levels = METRIC_LABEL)) %>%
  transmute(label, `Stock change` = stock_change_Tg / 1e3, `Avoided loss` = avoided_loss_Tg / 1e3,
            sc_lo = stock_change_Pg_lo, sc_hi = stock_change_Pg_hi, al_lo = avoided_loss_Pg_lo, al_hi = avoided_loss_Pg_hi)
gl <- g %>% pivot_longer(c(`Stock change`, `Avoided loss`), names_to = "quantity", values_to = "Pg") %>%
  mutate(lo = ifelse(quantity == "Stock change", sc_lo, al_lo), hi = ifelse(quantity == "Stock change", sc_hi, al_hi))
f6 <- ggplot(gl, aes(x = label, y = Pg, fill = label)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, na.rm = TRUE) +
  facet_wrap(~ quantity, scales = "free_y") +
  scale_fill_manual(values = c("grey60", "#2166ac", "#92c5de", "#053061"), guide = "none") +
  labs(x = NULL, y = "Pg AGB, T1 to T2", title = "Global protected-area biomass: stock change and avoided loss",
       caption = "Error bars: bootstrap 95% CI (uncorrected and central only). Low/high: 95% bounds of the drift rate.") +
  theme_fig() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
f6
save_fig(f6, "fig6_global_sensitivity", 10, 5.5)
d$global %>% filter(metric %in% names(METRIC_LABEL)) %>%
  select(metric, stock_t1_Tg, stock_change_Tg, avoided_loss_Tg, did_mean, pct_area_qualifying) %>% print()

pa <- d$pa %>% filter(metric == MAIN, qualifies)
f7 <- ggplot(pa, aes(x = pa_t1, y = did)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_hex(bins = 45) +
  geom_smooth(method = "loess", se = FALSE, colour = "black", linewidth = 0.7, formula = y ~ x) +
  scale_fill_viridis_c(option = "F", trans = "log10", direction = -1, name = "PAs") +
  facet_wrap(~ continent) +
  labs(x = "PA mean AGBD in T1 (Mg/ha)", y = "DiD (Mg/ha)", title = "Protection effect by starting biomass") +
  theme_fig()
f7
save_fig(f7, "fig7_did_vs_t1", 11, 7)