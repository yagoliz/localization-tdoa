%% Check drifts

%% Pre-alignment
load('data/dcf77_pre.mat');
t1 = input(1).t;
t2 = input(2).t;
t3 = input(3).t;

maxLength = min([length(t1), length(t2), length(t3)]);
% Drifts
td12 = t1(1:maxLength) - t2(1:maxLength);
td12 = (td12 - td12(1))*1e6;
td13 = t1(1:maxLength) - t3(1:maxLength);
td13 = (td13 - td13(1))*1e6;
td23 = t2(1:maxLength) - t3(1:maxLength);
td23 = (td23 - td23(1))*1e6;

time = (0:length(td12)-1)/12000;
figure(); hold on; grid on;
plot(time, td12); plot(time, td13); plot(time, td23);
legend('1 & 2', '1 & 3', '2 & 3');
title('Drift in time between different pairs (f_s = 12kHz)');
xlabel('Time (s)'); ylabel('Time difference (\mu{}s)');
xlim([time(1) time(end)]);

%% Drift per second
mean12 = 1e6*abs(td12(1) - td12(end))/30;
mean13 = 1e6*abs(td13(1) - td13(end))/30;
mean23 = 1e6*abs(td23(1) - td23(end))/30;