# Script version of figures/01_maps.ipynb (run from the figures/ folder: Rscript scripts/01_maps.R)
source("fig_setup.R")
options(repr.plot.width = 14, repr.plot.height = 8, repr.plot.res = 120)
d <- load_inputs()
world <- get_countries()
globe <- globe_outline()
base_map <- list(geom_sf(data = globe, fill = "#eef3f7", colour = "grey75", linewidth = 0.2),
                 geom_sf(data = world, fill = "grey92", colour = "white", linewidth = 0.1))
cm <- d$country %>% filter(metric == MAIN) %>%
  mutate(no_cf = pct_pas_qualifying < MIN_PCT_COUNTERFACTUAL,
         stock_change_pct = 100 * stock_change_Tg / stock_t1_Tg)
wd <- world %>% left_join(cm, by = "iso3")
cat(nrow(cm), "countries with results;", sum(cm$no_cf), "flagged as lacking counterfactuals\n")

lim <- quantile(abs(cm$did_mean[!cm$no_cf]), 0.98, na.rm = TRUE)
f1 <- ggplot() +
  base_map +
  geom_sf(data = filter(wd, !is.na(did_mean) & !no_cf), aes(fill = did_mean), colour = "white", linewidth = 0.1) +
  geom_sf(data = filter(wd, no_cf), fill = "grey70", colour = "grey40", linewidth = 0.25) +
  scale_fill_gradient2(low = "#8c510a", mid = "#f5f5f5", high = "#01665e", midpoint = 0,
                       limits = c(-lim, lim), oob = squish,
                       name = "DiD in AGBD change, T1 to T2 (Mg/ha)", guide = guide_colourbar(title.position = "top")) +
  labs(title = "Protected-area effect on biomass change",
       caption = "Grey with dark outline: too few PAs with counterfactuals. Light grey: no PA results.") +
  theme_map()
f1
save_fig(f1, "fig1_did_map", 12, 6.5)

f2a <- ggplot() +
  base_map +
  geom_sf(data = filter(wd, !is.na(stock_t1_Tg)), aes(fill = stock_t1_Tg), colour = "white", linewidth = 0.1) +
  scale_fill_viridis_c(option = "D", trans = "log10", labels = label_number(big.mark = ","),
                       name = "AGB stock in PAs, 2019-2021 (Tg)", guide = guide_colourbar(title.position = "top")) +
  labs(title = "a  Biomass stock in protected areas") + theme_map()
clim <- quantile(abs(cm$stock_change_pct), 0.98, na.rm = TRUE)
f2b <- ggplot() +
  base_map +
  geom_sf(data = filter(wd, !is.na(stock_change_pct)), aes(fill = stock_change_pct), colour = "white", linewidth = 0.1) +
  scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac", midpoint = 0,
                       limits = c(-clim, clim), oob = squish,
                       name = "Change in PA biomass stock, T1 to T2 (%)", guide = guide_colourbar(title.position = "top")) +
  labs(title = "b  Change in protected-area biomass stock (drift-corrected)") + theme_map()
f2 <- f2a / f2b
f2
save_fig(f2, "fig2_stock_maps", 11, 11)

f3 <- ggplot() +
  base_map +
  geom_sf(data = filter(wd, !is.na(pct_pas_qualifying)), aes(fill = pct_pas_qualifying), colour = "white", linewidth = 0.1) +
  geom_sf(data = filter(wd, no_cf), fill = NA, colour = "#d7301f", linewidth = 0.5) +
  scale_fill_viridis_c(option = "G", direction = -1, limits = c(0, 100),
                       name = "PAs with their own matched counterfactual (%)", guide = guide_colourbar(title.position = "top")) +
  labs(title = "Counterfactual availability",
       caption = paste0("Red outline: < ", MIN_PCT_COUNTERFACTUAL, "% of PAs with counterfactuals.")) +
  theme_map()
f3
save_fig(f3, "fig3_matching_success_map", 12, 6.5)
cm %>% filter(no_cf) %>% select(iso3, n_pas, n_qualifying, pct_pas_qualifying, pct_area_qualifying) %>%
  arrange(pct_pas_qualifying) %>% print(n = 50)