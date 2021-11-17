# makes mona_input_dfs.RData, which contains 
# - original_densities: a data frame containing "densities: from the hi-res image of ML
# - mona_df: a data frame containing x and y co-ords for a low-res version of ML
# --- D: the "true" densities
# --- Dgood, Dblur: covariates in decreasing order of quality

library(dplyr)
library(ggplot2)
library(imager)

# load hi-res image of part of mona lisa created by make_hires_mona_jpg.R
mona <- load.image("output/hires_mona.jpg")
plot(mona)
# save as a data frame, need to flip on y-axis (from imager) and rescale to resized_x
original_densities <- as.data.frame(mona) %>% mutate(y = (max(y) + 1 - y) / 24,
                                                     x = x / 24)

# make a low res version, this is the true density surface
resized_x <- 50
small_mona <- resize(mona, resized_x, resized_x)
plot(small_mona)
table(small_mona) # 2 cells are zero, this will cause problems with logs later on; next smallest is 0.003
small_mona[small_mona == 0] <- 0.001

# make covariate surfaces

# good covariate has a small blur
good_small_mona <- isoblur(small_mona,1) 
plot(good_small_mona)

# moderate covariate has a bigger blur
blurry_small_mona <- isoblur(small_mona,4) 
plot(blurry_small_mona)

# turn into data frame, D is not yet true density just pixel intensity

good_small_mona <- as.data.frame(good_small_mona) %>% 
  mutate(y = max(y) - y + 1, Dgood_unstd = value) %>% dplyr::select(-value)

blurry_small_mona <- as.data.frame(blurry_small_mona) %>% 
  mutate(y = max(y) - y + 1, Dblur_unstd = value) %>% dplyr::select(-value)

mona_df <- as.data.frame(small_mona) %>% 
  mutate(y = max(y) - y + 1, 
         D_unstd = value) %>% dplyr::select(-value) %>%
  left_join(good_small_mona) %>%
  left_join(blurry_small_mona) 

# check correlations
cor(mona_df)

# scale D and covariates so that these can be interpreted as densities
# 1) units are animals / ha
# 2) covariates scaled the same number of expected activity centres as the true density surface

# each cell is 1m2 and D in secr is animals / 10000m2
# so expected number of points in each cell is D / 10000 and 
# expected total points generated is sum(D) / 10000. 
# If we want n_pts points then we need to multiply D by 10000 * (n_pts / sum(D))
# Works same way for covariate surfaces, except sum(D) becomes sum(Dgood) or sum(Dblur)
# We want 7500 in the many animals scenario, 80 in the few

mona_df <- mona_df %>%
  mutate(D_bigD = D_unstd * 7500 / sum(D_unstd) * 10000,
         Dgood_bigD = Dgood_unstd * 7500 / sum(Dgood_unstd) * 10000,
         Dblur_bigD = Dblur_unstd * 7500 / sum(Dblur_unstd) * 10000,
         D_smallD = D_unstd * 80 / sum(D_unstd) * 10000,
         Dgood_smallD = Dgood_unstd * 80 / sum(Dgood_unstd) * 10000,
         Dblur_smallD = Dblur_unstd * 80 / sum(Dblur_unstd) * 10000)

mona_df %>% summarize_all(mean)
mona_df %>% summarize_all(sum)
mona_df %>% summarize_all(min)
mona_df %>% summarize_all(max)

ggplot(mona_df, aes(x = x, y = y, fill = D_bigD)) + geom_tile() + scale_fill_viridis_c()
ggplot(mona_df, aes(x = x, y = y, fill = Dgood_bigD)) + geom_tile() + scale_fill_viridis_c()
ggplot(mona_df, aes(x = x, y = y, fill = Dblur_bigD)) + geom_tile() + scale_fill_viridis_c()

save(original_densities, mona_df, file = "output/mona_inputs.RData")

