%% Here we will analyze the correlation offsets
clear; clc; close all;

% Load data from DCF77
load('data/dcf77.mat');

fs = input(1).fs;
% Load data by sensor
iq1 = input(1).z;
td1 = input(1).t;

iq2 = input(2).z;
td2 = input(2).t;

iq3 = input(3).z;
td3 = input(3).t;

% Largest index to analyze
jump = 10000;
maxValue = 300000;
numRows = maxValue/jump;

%% Correlations

% IQ
iqcorrelate12 = zeros(numRows,1);
data = 1;
for i = 1:jump:maxValue-jump+1
    [r, lags] = xcorr(iq1(i:i + (jump-1)), iq2(i:i + (jump-1)), round(6371e3*pi/fs), 'coeff');
    [~,idx] = max(abs(r));
    iqcorrelate12(data) = lags(idx);
    data = data + 1;
end

iqcorrelate12 = iqcorrelate12 / 12000 * 1e3;
plot(iqcorrelate12); xlabel('Chunk number'); ylabel('Delay (ms)'); grid on;
xlim([1,30]);
ylim([-1,0]); title('With time alignment');
% title('No time alignment');