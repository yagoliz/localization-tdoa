%% PPM test for a single RTL-SDR
% Author: Yago Lizarribar
%% Clean up the workspace
clearvars; clc;

%% Main variables
addpath('ltess');
samplingRate = 1.92e6;
fRS = 806e6;

%% Load the signals
folder = '/home/yago/Imdea/git/localization/localization-scripts/data/ppm_tests/';
filename = 'D0-E1.dat';

signal = spec_load([folder, filename]);
numChunks = length(signal)/samplingRate;
samplingVector = 0:samplingRate-1;
Ts = 1/samplingRate;
t = 0:Ts:Ts*(samplingRate - 1);

ppmVector = zeros(1, numChunks);
ppmVectorCorrected = zeros(1, numChunks);

%% Main loop
firstChunk = signal(1:samplingRate);
[ppmOriginal, ~] = ltess(firstChunk, samplingRate);
delta_f = ppmOriginal * 1e-6 * fRS;
samplingVectorOffset = samplingVector - samplingVector.*(ppmOriginal*1e-6);

for i = 1:numChunks
    disp(['Iteration: ', num2str(i)]);
    % Each iteration we analyze one chunk of 1 second duration
    signalChunk = signal((i-1)*samplingRate+1:i*samplingRate);
    [ppmVector(i), ~] = ltess(signalChunk, samplingRate);
    
    % Now we correct the signal and calculate the PPM for the corrected
    % chunk
%     signalChunk = signalChunk .* exp(-1i * 2 * pi .* t * ppmOriginal * (1e-6) * fRS)';
    signalCorrected = interp1(samplingVectorOffset, signalChunk, samplingVector)';
    [ppmVectorCorrected(i), ~] = ltess(signalCorrected, samplingRate);
end

%% Plotting and saving
pause(1);