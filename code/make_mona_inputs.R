# makes mona_input_dfs.RData, which contains 
# - original_densities: a data frame containing "densities: from the hi-res image of ML
# - mona_df: a data frame containing x and y co-ords for a low-res version of ML
# --- D: the "true" densities
# --- Dgood, Dblur, Dshift, Drept: covariates in decreasing order of quality

library(tidyverse)
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

# make covariate surfaces

# good covariate has a small blur
good_small_mona <- isoblur(small_mona,1) 
plot(good_small_mona)

# moderate covariate has a bigger blur
blurry_small_mona <- isoblur(small_mona,4) 
plot(blurry_small_mona)

# weak covariate shifts the moderate covariate around, with wrapping
blurry_shifted_small_mona <- imshift(blurry_small_mona,10,10,boundary=2) 
plot(blurry_shifted_small_mona)

# locally strong covariate is a mix of good covariate in top-right of image
# and weak covariate elsewehere
half_mona <- imappend(list(imsub(good_small_mona, y <= 2*height/3), 
                                           imsub(blurry_shifted_small_mona, y > 2*height/3)), 'y') 
plot(half_mona)

local_small_mona <- imappend(list(imsub(blurry_shifted_small_mona, x <= width/2), 
                                      imsub(half_mona, x > width/2)), 'x')  
plot(local_small_mona)

# Leonardo da Vinci as a covariate

# load hi-res image of part of mona lisa created by make_hires_mona_jpg.R
ldv <- load.image("output/hires_ldv.jpg")
plot(ldv)

# make a low res version, this is the true density surface
small_ldv <- resize(ldv, resized_x, resized_x, size_c = 1)
plot(small_ldv)

# turn into data frames

good_small_mona <- as.data.frame(good_small_mona) %>% 
  mutate(y = max(y) - y + 1, Dgood = value) %>% select(-value)

blurry_small_mona <- as.data.frame(blurry_small_mona) %>% 
  mutate(y = max(y) - y + 1, Dblur = value) %>% select(-value)

blurry_shifted_small_mona <- as.data.frame(blurry_shifted_small_mona) %>% 
  mutate(y = max(y) - y + 1, Dshift = value) %>% select(-value)

local_small_mona <- as.data.frame(local_small_mona) %>% 
  mutate(y = max(y) - y + 1, Drept = value) %>% select(-value)

small_ldv <- as.data.frame(small_ldv) %>% 
  mutate(y = max(y) - y + 1, Dldv = value) %>% select(-value)

mona_df <- as.data.frame(small_mona) %>% 
  mutate(y = max(y) - y + 1, 
         D = value) %>% select(-value) %>%
  left_join(good_small_mona) %>%
  left_join(blurry_small_mona) %>%
  left_join(blurry_shifted_small_mona) %>%
  left_join(local_small_mona) %>%
  left_join(small_ldv)

# check correlations
cor(mona_df)

save(original_densities, mona_df, file = "output/mona_inputs.RData")

