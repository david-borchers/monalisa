library(secr)
library(dplyr)

load("output/mona_inputs.RData")
load("output/capthists.RData")
load("output/simulated_densities.Rdata")

dx = 8
dy = 8
nx = 7
ny = 7
xorig = 1.5
yorig = 1.5
sigma = 4 
lambda0 = 0.1
noccasions = 100

mlmesh <- read.mask(data = mona_df)

simulated_points <- simulated_points_lots
n_pts <- nrow(simulated_points) 

dx <- dx
dy <- dy

# make a grid of detectors
detectors <- make.grid(nx = nx, ny = ny, spacex = dx, spacey = dy,
                       originxy = c(xorig, yorig), detector = "count")

## Generate capture histories
lambda0 <- lambda0
sigma <- sigma
D <- n_pts * 4 #(as 50x50 area is 1/4 of a hectare)

# # if loading existing capthist
# ct <- capthists_realised_and_expected_acd_few %>%
#   dplyr::filter(secr.fitformula == "D~1", noccasions == max(noccasions))
# capthist <- ct$capthist[[1]]
# detectors <- traps(capthist)

capthist <- NULL
if(!is.null(capthist)){
  capture_history_max_occasions <- capthist
} else { 
  capture_history_max_occasions <- sim.capthist(detectors, popn = simulated_points, detectfn = "HHN",
                                                detectpar = list(lambda0 = lambda0, sigma = sigma),
                                                noccasions = max(noccasions),
                                                nsessions = 1)
}
summary(capture_history_max_occasions, terse = TRUE)

i <- noccasions
capture_history <- subset(capture_history_max_occasions, occasions = 1:i)

summary(capture_history, terse = TRUE)

n <- dim(capture_history)[1]

# fit model specified by secr.fitformula
cfit2 <- secr.fit(capture_history, 
                  model = list(D~1, lambda0~1, sigma~1),  
                  mask = mlmesh, detectfn = "HHN", 
                  start = list(D = D, lambda0 = lambda0, sigma = sigma),
                  ncores = 1, method = "none", trace = TRUE)

# other stuff
tmp8 <- fx.total(cfit2)
plot(tmp8, covariate = 'D.sum', col = hcl.colors(13, "Blues"), plottype = 'shaded', legend = TRUE)
plot(detectors, add = TRUE)
plot(simulated_points, add = TRUE)

detectors6 <- make.grid(nx = 3, ny = 3, spacex = 8, spacey = 8, originxy = c(10, 16), detector = "count")
detectors7 <- make.grid(nx = 3, ny = 3, spacex = 8, spacey = 8, originxy = c(26, 16), detector = "count")
plot(tmp6, covariate = 'D.sum', col = hcl.colors(13, "Blues"), plottype = 'shaded', legend = TRUE)
plot(detectors6, add = TRUE)
segments(x0=18,x1=26,y0=8,y1=8,lty=2,col="black")
segments(x0=18,x1=26,y0=16,y1=16,lty=2,col="black")
segments(x0=18,x1=18,y0=8,y1=16,lty=2,col="black")
segments(x0=26,x1=26,y0=8,y1=16,lty=2,col="black")
plot(tmp7, covariate = 'D.sum', col = hcl.colors(13, "Blues"), plottype = 'shaded', legend = TRUE)
plot(detectors7, add = TRUE)
segments(x0=18,x1=26,y0=8,y1=8,lty=2,col="black")
segments(x0=18,x1=26,y0=16,y1=16,lty=2,col="black")
segments(x0=18,x1=18,y0=8,y1=16,lty=2,col="black")
segments(x0=26,x1=26,y0=8,y1=16,lty=2,col="black")

d0_many_3x3v1 <- tmp6
det_3x3v1 <- detectors6
d0_many_3x3v2 <- tmp7
det_3x3v2 <- detectors7
d0_many_7x7 <- tmp8
det_7x7 <- detectors
save(d0_many_3x3v1, d0_many_3x3v2, d0_many_7x7, det_3x3v1, det_3x3v2, det_7x7, file = "output/temp-newml-demo.RData")

d0_many <- rbind(data.frame(x = d0_many_3x3v1$x, y = d0_many_3x3v1$y, array_size = "3x3_1", D = c(covariates(d0_many_3x3v1)$D.sum)),
                 data.frame(x = d0_many_3x3v2$x, y = d0_many_3x3v2$y, array_size = "3x3_2", D = c(covariates(d0_many_3x3v2)$D.sum)),
                 data.frame(x = d0_many_7x7$x, y = d0_many_7x7$y, array_size = "7x7", D = c(covariates(d0_many_7x7)$D.sum)))

dets <- rbind(data.frame(x = det_3x3v1$x, y = det_3x3v1$y, array_size = "3x3_1"),
              data.frame(x = det_3x3v2$x, y = det_3x3v2$y, array_size = "3x3_2"),
              data.frame(x = det_7x7$x, y = det_7x7$y, array_size = "7x7"))

common_area <- data.frame(xmin = 18, xmax = 26, ymin = 8, ymax = 16)

p2a <- d0_many %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = D)) +
  scale_fill_gradientn(colours = hcl.colors(13, palette = "Blues")) +
  facet_wrap(. ~ array_size) +
  geom_rect(data = common_area, inherit.aes = F, colour = "white", fill = NA, size = 0.5, linetype = 2,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_point(data = dets, inherit.aes = T,
             colour = "red", pch = 4, size = 2) +
  coord_equal() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="right", legend.key.height = unit(0.7,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p2a
