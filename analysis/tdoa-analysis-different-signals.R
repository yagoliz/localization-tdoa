library(R.matlab)
library(ggplot2)

# Files for data analysis
folder = "~/Imdea/git/localization/localization-tdoa/matlab/results/"
signal = "196806"
num_samples = "2000000"
type = "dphase"

file_original = paste(folder, paste(signal, num_samples, "original", paste(type, ".mat", sep=""), sep="_"), sep="")
file_fo = paste(folder, paste(signal, num_samples,"fo_correction_no_interp", paste(type, ".mat", sep=""), sep="_"), sep="")
#file_fo_interp10 = paste(folder, paste(signal, num_samples,"fo_correction_no_interp", paste(type, ".mat", sep=""), sep="_"), sep="")
file_fo_interp10 = paste(folder, paste(signal, num_samples, "fo_correction_10_interp", paste(type, ".mat", sep=""), sep="_"), sep="")

# Resolution for the sampling rate
resolution = 3e8 / 2e6

# Opening file with data that was not corrected with LTESS
data_orig = readMat(file_original)
data_fo = readMat(file_fo)
data_fo_interp10 = readMat(file_fo_interp10)

# We introduce everything into a dataframe
doa_meters_orig = abs(data_orig$doa.meters.2) + resolution
doa_meters_orig_avg = abs(data_orig$doa.meters) + resolution
doa_meters_fo = abs(data_fo$doa.meters.2) + resolution
doa_meters_fo_interp10 = abs(data_fo_interp10$doa.meters.2) + resolution/10

# Correct for extreme outliers
doa_meters_orig_avg[doa_meters_orig_avg > 2000] = 2000
doa_meters_orig[doa_meters_orig > 2000] = 2000
doa_meters_fo[doa_meters_fo > 2000] = 2000
doa_meters_fo_interp10[doa_meters_fo_interp10 > 2000] = 2000

# Let the plotting begin
labels = c(rep("Panoradio (Original)", length(doa_meters_fo)), rep("Panoradio + SOC (mean)", length(doa_meters_orig_avg)), rep("CFO + SFO", length(doa_meters_fo)), rep("CFO + SFO + 10 interp", length(doa_meters_fo_interp10)))
doa = c(doa_meters_orig, doa_meters_orig_avg, doa_meters_fo, doa_meters_fo_interp10)

data_frame = data.frame(labels, doa)

library(tidyverse)
dataMedian = aggregate(data_frame$doa, list(data_frame$labels), median)
dataMedian = data.frame(dataMedian, c(rep("darkgrey", 4)))

title = ""
if (substr(signal, 1, 3) == "196") {
  title = "DAB"
} else if (substr(signal, 1, 3) == "610") {
  title = "DVB-T"
} else if (substr(signal, 1, 3) == "100") {
  title = "FM"
} else if (substr(signal, 1, 3) == "806") {
  title = "LTE"
}

# Ggplot magic
p <- ggplot(data_frame, aes(x = reorder(labels, doa, median), y = doa, color = reorder(labels, doa, median))) + geom_boxplot(outlier.shape=16, outlier.size=0, notch=FALSE) 
p <- p + xlab("Method") + ylab("TDOA (m)") + theme(legend.position="none") + ylim(0, 900) + ggtitle(paste("Results with", title, sep=" "))
#p <- p + geom_text(dataMedian, mapping = aes(x = reorder(Group.1, -x, median), y = x, color = reorder(Group.1, -x, median),label=sprintf("%0.2f", round(x, digits = 2))), position = position_dodge(width = 1), size = 15, vjust = -1)
#p + theme(axis.text = element_text(size = 40)) + theme(axis.title.x = element_text(size = 40)) + theme(axis.title.y = element_text(size = 40))

p + geom_text(dataMedian, mapping = aes(x = reorder(Group.1, -x, median), y = x, color = reorder(Group.1, -x, median),label=sprintf("%0.2f", round(x, digits = 2))), position = position_dodge(width = 1), vjust = -1)