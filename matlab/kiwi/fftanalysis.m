%% Analyze DCF77 and GYN2BR

clear all; close all;
% DCF77
load('data/dcf77_1_pre.mat');

iqdata = input(3).z;
tdata = input(3).t - input(3).t(1);

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
load('data/dcf77_2_pre.mat');

iqdata = input(3).z;
tdata = input(3).t - input(3).t(1);

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