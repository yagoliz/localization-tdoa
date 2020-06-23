library(R.matlab)
library(ggplot2)

# Files for data analysis
folder = "~/Imdea/git/localization/localization-tdoa/matlab/results_mlat/"

## IQ vs ABS
file_iq   = paste(folder, "iq1.mat" , sep="")
file_iq5  = paste(folder, "iq5.mat" , sep="")
file_iq10 = paste(folder, "iq10.mat", sep="")
file_iq20 = paste(folder, "iq20.mat", sep="")
file_abs  = paste(folder, "abs1.mat", sep="")
file_kiwi = paste(folder, "kiwi.mat", sep="")

# Opening file with data that was not corrected with LTESS
data_iq   = readMat(file_iq  )
data_iq5  = readMat(file_iq5 )
data_iq10 = readMat(file_iq10)
data_iq20 = readMat(file_iq20)
data_abs  = readMat(file_abs )
data_kiwi = readMat(file_kiwi)

# We introduce everything into a dataframe
errors_iq = data_iq$errors1/1000
errors_iq5 = data_iq5$errors5/1000
errors_iq10 = data_iq10$errors10/1000
errors_iq20 = data_iq20$errors20/1000
errors_abs = data_abs$errors1/1000
errors_kiwi = data_kiwi$errorsk;

# Let the plotting begin
labels = c(rep("IQ", length(errors_iq)), rep("ABS", length(errors_abs)))
errors = c(errors_iq, errors_abs)

number_ups = c(rep("1", length(errors_iq)), rep("2", length(errors_iq5)), rep("3", length(errors_iq10)), rep("4", length(errors_iq20)))
#labels_ups = c(rep("No Upsampling", length(errors_iq)), rep("5 Upsampling", length(errors_iq5)), rep("10 Upsampling", length(errors_iq10)), rep("20 Upsampling", length(errors_iq20)))
labels_ups = c("No Upsampling", "5", "10", "20")
errors_ups = c(errors_iq, errors_iq5, errors_iq10, errors_iq20)

number_kw = c(rep("1", length(errors_kiwi)), rep("2", length(errors_iq)))
labels_kw = c("KiwiSDR", "Fang's Algorithm")
errors_kw = c(errors_kiwi, errors_iq)

data_frame = data.frame(labels, errors)
data_frame_ups = data.frame(number_ups, labels_ups, errors_ups)
data_frame_kw = data.frame(number_kw, errors_kw)

library(tidyverse)
dataMedian = aggregate(data_frame$errors, list(data_frame$labels), median)
dataMedian = data.frame(dataMedian, c(rep("darkgrey", 2)))

dataMedianUps = aggregate(data_frame_ups$errors_ups, list(number_ups), median)
dataMedianUps = data.frame(dataMedianUps, c(rep("darkgrey", 4)))

dataMedianKw = aggregate(data_frame_kw$errors_k, list(number_kw), median)
dataMedianKw = data.frame(dataMedianKw, c(rep("darkgrey", 2)))

title = "IQ vs ABS correlation"
title_ups = "Different Upsampling factors"
title_kw = "Kiwi vs Fang's Algorithm"

# Ggplot magic
p <- ggplot(data_frame, aes(x = reorder(labels, errors, median), y = errors, color = reorder(labels, errors, median))) + geom_boxplot(outlier.shape=16, outlier.size=0, notch=FALSE) 
p <- p + xlab("Method") + ylab("Error (km)") + theme(legend.position="none") + ylim(0, 30) + ggtitle(title)
p <- p + geom_text(dataMedian, mapping = aes(x = reorder(Group.1, -x, median), y = x, color = reorder(Group.1, -x, median),label=sprintf("%0.2f", round(x, digits = 2))), position = position_dodge(width = 1), size = 15, vjust = -1)
p + theme(axis.text = element_text(size = 40)) + theme(axis.title.x = element_text(size = 40)) + theme(axis.title.y = element_text(size = 40))

#p + geom_text(dataMedian, mapping = aes(x = reorder(Group.1, -x, median), y = x, color = reorder(Group.1, -x, median),label=sprintf("%0.2f", round(x, digits = 2))), position = position_dodge(width = 1), vjust = -1)

pups <- ggplot(data_frame_ups, aes(x = number_ups, y = errors_ups, color = number_ups)) + geom_boxplot(outlier.shape=16, outlier.size=0, notch=FALSE) 
pups <- pups + xlab("Upsampling Factor") + ylab("Error (km)") + theme(legend.position="none") + ylim(0, 30) + ggtitle(title_ups)
pups <- pups + geom_text(dataMedianUps, mapping = aes(x = Group.1, y = x, color = Group.1,label=sprintf("%0.2f", round(x, digits = 2))), position = position_dodge(width = 1), size = 15, vjust = -1)
pups <- pups + scale_x_discrete(breaks=c("1", "2", "3", "4"), labels=labels_ups)
pups + theme(axis.text = element_text(size = 40)) + theme(axis.title.x = element_text(size = 40)) + theme(axis.title.y = element_text(size = 40))

pkw <- ggplot(data_frame_kw, aes(x = number_kw, y = errors_kw, color = number_kw)) + geom_boxplot(outlier.shape=16, outlier.size=0, notch=FALSE) 
pkw <- pkw + xlab("Multilateration Method") + ylab("Error (km)") + theme(legend.position="none") + ylim(0, 30) + ggtitle(title_kw)
pkw <- pkw + geom_text(dataMedianKw, mapping = aes(x = Group.1, y = x, color = Group.1,label=sprintf("%0.2f", round(x, digits = 2))), position = position_dodge(width = 1), size = 15, vjust = -1)
pkw <- pkw + scale_x_discrete(breaks=c("1", "2"), labels=labels_kw)
pkw + theme(axis.text = element_text(size = 40)) + theme(axis.title.x = element_text(size = 40)) + theme(axis.title.y = element_text(size = 40))

