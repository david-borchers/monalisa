library(ggplot2)
library(dplyr)
library(secr)
library(stringr)
library(patchwork)

source("code/predicted_densities_for_D0.R")
source("code/run_secr.R")

load("output/simulated_densities.Rdata")
load("output/mona_inputs.RData")

simulated_points <- simulated_points_few
mlmesh <- read.mask(data = mona_df)

n_pts <- nrow(simulated_points) 

# make a 7x7 grid of detectors
detectors7x7 <- make.grid(nx = 7, ny = 7, spacex = 8, spacey = 8,
                          originxy = c(1.5, 1.5), detector = "count")

sigma = 4
lambda0 = 0.1
max_nocc <- 1201
my.seed = 112233
capture_history_max_occasions <- sim.capthist(detectors7x7, popn = simulated_points, detectfn = "HHN",
                                              detectpar = list(lambda0 = lambda0, sigma = sigma),
                                              noccasions = max_nocc,
                                              nsessions = 1,
                                              seed = my.seed)


parlist33 <- tibble(i=1:1,
                    secr.fitformula = c("D~1"),
                    dx = 8, dy = 8, nx = 7, ny = 7, 
                    xorig = 1.5, yorig = 1.5, 
                    sigma = 4, lambda0 = 0.1, 
                    noccasions = max_nocc,
                    stringsAsFactors = FALSE)
parlist33 <- parlist33 %>% mutate(capthist = list(capture_history_max_occasions))

i <- 1
m1 <- run_secr(simulated_points = simulated_points_few,
               secr.fitformula = parlist33$secr.fitformula[i], 
               dx = parlist33$dx[i], dy = parlist33$dy[i], 
               nx = parlist33$nx[i], ny = parlist33$ny[i], 
               xorig = parlist33$xorig[i], yorig = parlist33$yorig[i], 
               sigma = parlist33$sigma[i], lambda0 = parlist33$lambda0[i], 
               noccasions = parlist33$noccasions[i], 
               capthist = parlist33$capthist[[i]],
               return.model = TRUE)

tmp8 <- fx.total(m1$mod)
plot(tmp8, covariate = 'D.sum', col = hcl.colors(13, "Blues"), plottype = 'shaded', legend = FALSE)

tmp88 <- data.frame(x = tmp8$x, y = tmp8$y, value = covariates(tmp8)$D.sum)

newmask <- expand.grid(x = seq(from = 1, to = 50, length.out = 200),
                       y = seq(from = 1, to = 50, length.out = 200))
newmask <- read.mask(data = newmask, spacing = 49/500)

tmp9 <- fx.total(m1$mod, mask = newmask)
tmp99 <- data.frame(x = tmp9$x, y = tmp9$y, value = covariates(tmp9)$D.sum)

p2a <- tmp88 %>%
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

p2a

p2b <- tmp99 %>%
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
  #geom_point(data = simulated_points_few, inherit.aes = F, aes(x=x,y=y),
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

p2b

p2a + p2b

ggsave("paper/mona_rd_in_hres2.png", p2a+p2b, width=8, height=6, dpi = 600)
