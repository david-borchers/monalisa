library(tidyverse)
library(secr)

load("output/simulated_capturehistories.RData")

ret_ch <- function(x, nocc){
  
  chs <- summary(x)
  
  data.frame(det_x = chs$trapsum$xrange[1],
             det_y = chs$trapsum$yrange[1],
             occasion = nocc,
             n_animals = chs$counts[4,nocc],
             n_detections = sum(chs$counts[6,1:nocc]))
}

# capture histories for Figure 5
ch_fig5 <- map_df(.x = capturehistory_lots, .f = ret_ch, nocc = 1)
ch_fig5 <- ch_fig5 %>% mutate(detector = paste(ifelse(det_y == 15, "Bottom", "Top"),
                                               ifelse(det_x == 15, "Left", "Right")))

# capture histories for Figure 6
ch_fig6 <- map_df(.x = capturehistory_lots, .f = ret_ch, nocc = 20)
ch_fig6 <- ch_fig6 %>% mutate(detector = paste(ifelse(det_y == 15, "Bottom", "Top"),
                                               ifelse(det_x == 15, "Left", "Right")))

# capture history for Figure 7 is as for Figure 5 (1 occasion, many animals)

# capture histories for Figure 8
ch_fig8a <- ret_ch(x = capturehistory_few, nocc = 1)
ch_fig8b <- ret_ch(x = capturehistory_few, nocc = 3)
ch_fig8c <- ret_ch(x = capturehistory_few, nocc = 10)
ch_fig8d <- ret_ch(x = capturehistory_few, nocc = 20)
ch_fig8 <- rbind(ch_fig8a, ch_fig8b, ch_fig8c, ch_fig8d)
rm(ch_fig8a, ch_fig8b, ch_fig8c, ch_fig8d)

# capture history for Figure 9 is as for Figure 8 

save(ch_fig5, ch_fig6, ch_fig8, file = "output/capthist_summaries.RData")
