#load pedigrees and SNV data
source("data-raw/VignetteTimingData.R")


#data frame of times to populate for each set of settings
time_dat <- data.frame(time_ped10 = rep(NA, 10),
                       time_ped50 = rep(NA, 10),
                       time_ped100 = rep(NA, 10),
                       time_ped150 = rep(NA, 10),
                       time_ped200 = rep(NA, 10))

simStudy_settings <- data.frame(affOnly = c(TRUE, TRUE, FALSE, FALSE),
                                remWild = c(TRUE, FALSE, TRUE, FALSE))

num_peds = c(10, 50, 100, 150, 200)

#one entry in time_res for each set of settings in simStudy_settings
#the jth column in time_dat is for the jth setting in num_peds
#the ith row in time_dat is for the ith run (of 10)
time_res = list()
time_res[[1]] = time_dat
time_res[[2]] = time_dat
time_res[[3]] = time_dat
time_res[[4]] = time_dat

set.seed(2304789)
for(i in 1:length(num_peds)){
  for(k in 1:10){

    # Must do this complicated sampling, because we must ensure that each
    # sampled family has a unique famID.
    tpeds = lapply(1:num_peds[i], function(x){
      time_peds[time_peds$FamID == sample(time_peds$FamID, size = 1), ]
    })

    #re-name famID to ensure that each family has unique FamID
    for (np in 1:length(tpeds)){
      tpeds[[np]]$FamID <- np
    }

    #combine into a single dataframe
    tpeds = do.call(rbind, tpeds)

    #time sim_RVstudy for each set of settings
    for(j in 1:nrow(simStudy_settings)){
      a = Sys.time()
      study_seq <- sim_RVstudy(ped_files = tpeds,
                               SNV_data = sout,
                               affected_only = simStudy_settings$affOnly[[j]],
                               remove_wild = simStudy_settings$remWild[[j]])
      b = Sys.time()

      time_res[[j]][k, i] = difftime(b, a, units = "mins")

    }
  }
}

time_res1 = time_res[[1]]
time_res2 = time_res[[2]]
time_res3 = time_res[[3]]
time_res4 = time_res[[4]]



# save(time_res1, file = "C:/Data/SimSeqTiming/time_res1.rdata", compress='xz')
# save(time_res2, file = "C:/Data/SimSeqTiming/time_res2.rdata", compress='xz')
# save(time_res3, file = "C:/Data/SimSeqTiming/time_res3.rdata", compress='xz')
# save(time_res4, file = "C:/Data/SimSeqTiming/time_res4.rdata", compress='xz')



load(file = "C:/Data/SimSeqTiming/time_res1.rdata")
load(file = "C:/Data/SimSeqTiming/time_res2.rdata")
load(file = "C:/Data/SimSeqTiming/time_res3.rdata")
load(file = "C:/Data/SimSeqTiming/time_res4.rdata")


apply(time_res1, 2, mean)
apply(time_res1, 2, sd)

apply(time_res2, 2, mean)
apply(time_res2, 2, sd)

apply(time_res3, 2, mean)
apply(time_res3, 2, sd)

apply(time_res4, 2, mean)
apply(time_res4, 2, sd)

simStudy_settings <- data.frame(affOnly = c(TRUE, TRUE, FALSE, FALSE),
                                remWild = c(TRUE, FALSE, TRUE, FALSE))

num_peds = c(10, 50, 100, 150, 200)



my.cols = c(SFUred = rgb(red = 166/225, green = 25/225, blue = 46/225), #SFUred
            SFUblack = rgb(red = 84/225, green = 88/225, blue = 90/225), #SFUblack
            SFUblue = rgb(red = 0/225, green = 112/225, blue = 150/225), #SFUblue
            SFUgold = rgb(red = 193/225, green = 160/225, blue = 30/225),#SFUgold
            SFUgreen = rgb(red = 112/225, green = 140/225, blue = 50/225),
            SFUteal = rgb(red = 0, green = 128/225, blue = 128/225),
            SFUmidnight = rgb(red = 25/225, green = 25/225, blue = 112/225)) #SFUgreen

my_pch = c(19, 17, 15, 18)
win.graph(h = 6, w = 6)
plot(x = num_peds, y = apply(time_res1, 2, mean),
     xlab = "Number of Pedigrees",
     xaxt = "n",
     ylab = "Computation Time (minutes)",
     main = "Compuation Time by Number of Pedigrees \n and Simulation Settings",
     col = my.cols[[1]],
     pch = my_pch[1],
     lwd = 1,
     type = "b",
     ylim = c(0, 10))
#create legend to detail population
legend("topleft", #title = "settings",
       legend = c("affected_only = TRUE, remove_wild = TRUE",
                  "affected_only = TRUE, remove_wild = FALSE",
                  "affected_only = FALSE, remove_wild = TRUE",
                  "affected_only = FALSE, remove_wild = FALSE"),
       cex = 1,
       pch = c(my_pch[1], my_pch[2], my_pch[3], my_pch[4]),
       col = c(my.cols[1], my.cols[3], my.cols[6], my.cols[7]))

axis(side = 1, at = num_peds,
     labels = as.character(num_peds),
     cex = 1.2, line = 0, lwd = 1)

points(x = num_peds, y = apply(time_res2, 2, mean),
       col = my.cols[[3]],
       pch = my_pch[2],
     type = "b")
points(x = num_peds, y = apply(time_res3, 2, mean),
       col = my.cols[[6]],
       pch = my_pch[3],
       type = "b")
points(x = num_peds, y = apply(time_res4, 2, mean),
       col = my.cols[[7]],
       pch = my_pch[4],
       type = "b")



#-------------#
# with ggplot #
#-------------#
time_res1$settings = "affected_only = T, remove_wild = T"
time_res2$settings = "affected_only = T, remove_wild = F"
time_res3$settings = "affected_only = F, remove_wild = T"
time_res4$settings = "affected_only = F, remove_wild = F"

tdf <- rbind(melt(time_res1, id.vars = 'settings'),
             melt(time_res2, id.vars = 'settings'),
             melt(time_res3, id.vars = 'settings'),
             melt(time_res4, id.vars = 'settings'))

tdf$num_peds <- ifelse(tdf$variable == "time_ped10", 10,
                       ifelse(tdf$variable == "time_ped50", 50,
                              ifelse(tdf$variable == "time_ped100", 100,
                                     ifelse(tdf$variable == "time_ped150", 150, 200))))

#obtain required statistics
library(dplyr)
tdf2 <- tdf %>%
  group_by(settings, num_peds) %>%
  summarize(mean = mean(value),
            sd = sd(value),
            n = n())

#plot with error bars
library(ggplot2)
ggplot(data = tdf2, aes(x = num_peds, y = mean,
                            group = settings, colour = factor(settings))) +
  geom_line(size = 0.9) + geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = mean - 1.96*sd/sqrt(n), ymax = mean + 1.96*sd/sqrt(n)), width = 1) +
  coord_cartesian(xlim = c(0, 210), ylim = c(0, 7)) +
  labs(title = "Computation Time by Number of Pedigrees and Settings",
       y = "Computation Time (minutes)",
       x = "Number of Pedigrees")

