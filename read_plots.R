# Reading from Excel files
library(xlsx)
# Data visualisation
library(ggplot2)
getwd()
setwd("C:/Users/11901250/Documents/SF/BIT09_HTA/")
reads_data <- read.csv(file="reads.csv", header = TRUE)

library(vioplot) # install.packages("vioplot")

sw1 <- reads_data$
sw2 <- reads_data$reads.trimmed.M.
sw3 <- reads_data$reads.mapped.M.
sw4 <- reads_data$reads.after.counting.M.
vioplot(sw1, sw2, sw3, sw4, names=c("start", "trimmed", "mapped","counted"),
        cex = 0.5,
        ylim = c(5,30),
        col = brewer.pal(4,"Dark2"))

vioplot(sw1, names="start",
        ylim = c(15,30),
        col = "#1B9E77")
vioplot(sw2, names="trimmed",
        ylim = c(12.5,27),
        col = "#D95F02")
vioplot(sw3, names="mapped",
        ylim = c(12.5,25),
        col = "#7570B3")
vioplot(sw4, names="counted",
        ylim = c(5,15),
        col = "#E7298A")

mr1 <- mean(reads_data$percentage.trimmed)
mr2 <- mean(reads_data$percentage.mapped)
mr3 <- mean(reads_data$percentage.counting)

barplot(c(mr1,mr2,mr3),
        names=c("trimmed","mapped","counted"),
        col = c("#D95F02","#7570B3","#E7298A"))

