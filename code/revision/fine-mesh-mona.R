### fine mesh


tmp8 <- fx.total(res_acd_few[[13]]$mod)
tmp88 <- data.frame(x = tmp8$x, y = tmp8$y, value = covariates(tmp8)$D.sum)

newmask <- expand.grid(x = seq(0.5+.25/2,50.5-.25/2,0.25), y = seq(0.5+.25/2,50.5-.25/2,0.25))
newmask <- read.mask(data = newmask, spacing = 0.25)

tmp9 <- fx.total(res_acd_few[[13]]$mod, mask = newmask)
tmp99 <- data.frame(x = tmp9$x, y = tmp9$y, value = covariates(tmp9)$D.sum)

p2c <- tmp88 %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  #scale_fill_continuous_diverging(palette = "Blue-Red 3", l1 = 30, l2 = 100, p1 = .9, p2 = 1.2, mid = 0.225) +
  # scale_fill_gradientn(colours = c(hcl.colors(16, palette = "Blues"),hcl.colors(16, palette = "Reds", rev = TRUE)),
  #                      values = rescale(c(0,0.08,.45)),
  #                      guide = "colorbar",
  #                      limits=c(0,0.45)) +
  scale_fill_gradientn(colours = hcl.colors(16, palette = "Blues")) +
  #scale_fill_gradientn(colours = hcl.colors(16, palette = "hawaii")) +
  #scale_fill_viridis(direction = 1, option = "viridis") +
  # facet_grid(covtype2 ~ occasions2) +
  # geom_point(data = detectors_df_all %>% filter(occasions %in% occ[1:nn], array_size %in% asz), inherit.aes = T,
  #            colour = "gray80", pch = 4, size = 2) +
  # geom_point(data = simulated_points_few, inherit.aes = F, aes(x=x,y=y),
  #            colour = "darkorange", pch = 16, size = 1, alpha = 0.5) +
  # geom_point(data = r2, inherit.aes = F, aes(x=x,y=y), 
  #            colour = "red", size = .2) + 
  coord_equal() +
  theme_classic(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.spacing=unit(-1, "lines"),
        strip.background = element_rect(fill=NA, colour = NA), 
        legend.position="none", legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.7,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p2d <- tmp99 %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  #scale_fill_continuous_diverging(palette = "Blue-Red 3", l1 = 30, l2 = 100, p1 = .9, p2 = 1.2, mid = 0.225) +
  # scale_fill_gradientn(colours = c(hcl.colors(16, palette = "Blues"),hcl.colors(16, palette = "Reds", rev = TRUE)),
  #                      values = rescale(c(0,0.08,.45)),
  #                      guide = "colorbar",
  #                      limits=c(0,0.45)) +
  scale_fill_gradientn(colours = hcl.colors(16, palette = "Blues")) +
  #scale_fill_gradientn(colours = hcl.colors(16, palette = "hawaii")) +
  #scale_fill_viridis(direction = 1, option = "viridis") +
  # facet_grid(covtype2 ~ occasions2) +
  # geom_point(data = detectors_df_all %>% filter(occasions %in% occ[1:nn], array_size %in% asz), inherit.aes = T,
  #            colour = "gray80", pch = 4, size = 2) +
  # geom_point(data = simulated_points, inherit.aes = F, aes(x=x,y=y),
  #            colour = "darkorange", pch = 16, size = 1, alpha = 0.5) +
  # geom_point(data = r2, inherit.aes = F, aes(x=x,y=y), 
  #            colour = "red", size = .2) + 
  coord_equal() +
  theme_classic(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.spacing=unit(-1, "lines"),
        strip.background = element_rect(fill=NA, colour = NA), 
        legend.position="none", legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.7,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p2c + p2d

ggsave("paper/revision/mona_rd_in_hres.png", p2c+p2d, width=8, height=6, dpi = 600)
