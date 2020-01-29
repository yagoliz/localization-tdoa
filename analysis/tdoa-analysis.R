library(R.matlab)
library(ggplot2)

# Files for data analysis
folder = "~/github/tdoa-evaluation-rtlsdr/results/"
type = "dphase"

file_original = paste(folder, paste("original", paste(type, ".mat", sep=""), sep="_"), sep="")
file_fo = paste(folder, paste("fo_correction_no_interp", paste(type, ".mat", sep=""), sep="_"), sep="")
file_fo_interp5 = paste(folder, paste("lte_fo_correction_no_interp", paste(type, ".mat", sep=""), sep="_"), sep="")
file_fo_interp20 = paste(folder, paste("fo_correction_20_interp", paste(type, ".mat", sep=""), sep="_"), sep="")
file_fo_interp50 = paste(folder, paste("fo_correction_50_interp", paste(type, ".mat", sep=""), sep="_"), sep="")

# Resolution for the sampling rate
resolution = 3e8 / 1.92e6

# Opening file with data that was not corrected with LTESS
data_orig = readMat(file_original)
data_fo = readMat(file_fo)
data_fo_interp5 = readMat(file_fo_interp5)
data_fo_interp20 = readMat(file_fo_interp20)
data_fo_interp50 = readMat(file_fo_interp50)

# We introduce everything into a dataframe
doa_meters_orig = abs(data_orig$doa.meters) + resolution
doa_meters_fo = abs(data_fo$doa.meters.2) + resolution
doa_meters_fo_interp5 = abs(data_fo_interp5$doa.meters.2) + resolution
doa_meters_fo_interp20 = abs(data_fo_interp20$doa.meters.2) + resolution/20
doa_meters_fo_interp50 = abs(data_fo_interp50$doa.meters.2) + resolution/50
# Correct for extreme outliers
doa_meters_orig[doa_meters_orig > 1000] = 1000
doa_meters_fo[doa_meters_fo > 1000] = 1000
doa_meters_fo_interp5[doa_meters_fo_interp5 > 1000] = 1000
doa_meters_fo_interp20[doa_meters_fo_interp20 > 1000] = 1000
doa_meters_fo_interp50[doa_meters_fo_interp50 > 1000] = 1000

# Let the plotting begin
labels = c(rep("Panoradio", length(doa_meters_fo)), rep("interp1", length(doa_meters_fo)), rep("resample", length(doa_meters_fo_interp5)), rep("FOC + Int (20)", length(doa_meters_fo_interp20)), rep("FOC + Int (50)", length(doa_meters_fo_interp50)))
doa = c(doa_meters_orig, doa_meters_fo, doa_meters_fo_interp5, doa_meters_fo_interp20, doa_meters_fo_interp50)

data_frame = data.frame(labels, doa)

library(tidyverse)
dataMedian = aggregate(data_frame$doa, list(data_frame$labels), median)
dataMedian = data.frame(dataMedian, c(rep("darkgrey", 5)))

# Ggplot magic
p <- ggplot(data_frame, aes(x = reorder(labels, -doa, median), y = doa, color = reorder(labels, -doa, median))) + geom_boxplot(outlier.shape=16, outlier.size=0, notch=FALSE) 
p <- p + xlab("Method") + ylab("TDOA (m)") + theme(legend.position="none") + ylim(0, 1200)  + ggtitle("Results with interp1")
#p <- p + geom_text(dataMedian, mapping = aes(x = reorder(Group.1, -x, median), y = x, color = reorder(Group.1, -x, median),label=sprintf("%0.2f", round(x, digits = 2))), position = position_dodge(width = 1), size = 15, vjust = -1)
#p + theme(axis.text = element_text(size = 40)) + theme(axis.title.x = element_text(size = 40)) + theme(axis.title.y = element_text(size = 40))
p + geom_text(dataMedian, mapping = aes(x = reorder(Group.1, -x, median), y = x, color = reorder(Group.1, -x, median),label=sprintf("%0.2f", round(x, digits = 2))), position = position_dodge(width = 1), vjust = -1)
