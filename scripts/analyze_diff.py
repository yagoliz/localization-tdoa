#!/usr/bin/python3
import matplotlib.pyplot as plt

plt.style.use('ggplot')

data_rtl_0 = []
with open("device_0.txt") as f:
    for line in f:
        data_rtl_0.append(float(line[:-1]))

data_rtl_1 = []
with open("device_1.txt") as f:
    for line in f:
        data_rtl_1.append(float(line[:-1]))

data_rtl_2 = []
with open("device_2.txt") as f:
    for line in f:
        data_rtl_2.append(float(line[:-1]))

diff_01 = []
diff_02 = []
diff_12 = []
for i in range(len(data_rtl_0)):
    diff_01.append(1000*abs(data_rtl_0[i] - data_rtl_1[i]))
    diff_02.append(1000*abs(data_rtl_0[i] - data_rtl_2[i]))
    diff_12.append(1000*abs(data_rtl_1[i] - data_rtl_2[i]))

plt.plot(diff_01, 'r')
plt.plot(diff_02, 'g')
plt.plot(diff_12, 'b')
plt.legend(["RTL 0 vs 1", "RTL 0 vs 2", "RTL 1 vs 2"])
plt.ylabel("Difference (ms)")

plt.show()
