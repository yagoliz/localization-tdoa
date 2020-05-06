%% Analyze DCF77 and GYN2BR

clear all; close all;
% DCF77
load('data/dcf77_pre.mat');

iqdata = input(3).z;
tdata = input(3).t - tmin(3);

fc = input(3).freq;
freqs = linspace(-6, 6, 512) + fc;
y = zeros(length(iqdata)/512, 512);
pos = 1;
for i = 1:512:length(iqdata)-511
    y(pos,:) = 10 * log10(abs(fftshift(fft(iqdata(i:i+511)))));
    pos = pos + 1;
end

figure();
imagesc(freqs, tdata, y);
xlabel('Frequency (kHz)');
ylabel('Time (s)');

% Both
load('data/both_pre.mat');

iqdata = input(2).z;
tdata = input(2).t - tmin(2);

fc = input(3).freq;
freqs = linspace(-6, 6, 512) + fc;
y = zeros(length(iqdata)/512, 512);
pos = 1;
for i = 1:512:length(iqdata)-511
    y(pos,:) = 10 * log10(abs(fftshift(fft(iqdata(i:i+511)))));
    pos = pos + 1;
end

figure();
imagesc(freqs, tdata, y);
xlabel('Frequency (kHz)');
ylabel('Time (s)');