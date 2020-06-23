clear all; clc;
fRS_MHz = 938.8;
fUS_MHz = 806;
sampling_rate = 2e6;
num_samples_per_freq = 2000000;
folder_identifier_tdoa = ['/home/ygglc/Imdea/git/localization/localization-tdoa/data/localization/', num2str(fRS_MHz), '_', num2str(fUS_MHz), '_', num2str(num_samples_per_freq), '/', 'D0-E2-localization.dat'];

addpath('ltess');
signal = spec_load(folder_identifier_tdoa);

signalc = signal(1:1500000);

ppm = 0.0   ;
Ts = 1/sampling_rate;
phi = ppm * 1e-6;
t = 0:Ts/(1+phi):Ts/(1+phi)*(length(signalc) - 1);

signalccent = signalc .* (exp(-1i * 2 * pi .* t * phi * fRS_MHz * 1e6)');
signalmod = resample(signalccent, t, 2e6);

[r,lags] = xcorr(signalc, signalmod,1000);
[~,idx] = max(abs(r));
lags(idx)
