% =========================================================================
%  Batch test realization
%  Author: Yago Lizarribar
% =========================================================================

clear;
clc;
close all;

% adds subfolder with functions to PATH
[p,n,e] = fileparts(mfilename('fullpath'));
addpath([p '/tdoa']);
addpath([p '/ltess']);

%% Read Parameters from config qfile, that specifies all parameters
%---------------------------------------------
ppm_analysis_config;

% Get all files
filedev1_ltess = [folder_identifier, 'D0_ltess.dat'];
filedev2_ltess = [folder_identifier, 'D1_ltess.dat'];
filedev3_ltess = [folder_identifier, 'D2_ltess.dat'];

filedev1_tdoa = [folder_identifier, 'D0_ppm.dat'];
filedev2_tdoa = [folder_identifier, 'D1_ppm.dat'];
filedev3_tdoa = [folder_identifier, 'D2_ppm.dat'];

%% Main loop

% Read Signals from File
disp('______________________________________________________________________________________________');
disp('READ LTESS DATA FROM FILES');
signal1_ltess = spec_load(filedev1_ltess);
disp('File 1 Successfully loaded');
signal2_ltess = spec_load(filedev2_ltess);
disp('File 2 Successfully loaded');
signal3_ltess = spec_load(filedev3_ltess);
disp('File 3 Successfully loaded');

disp(' ');
disp('______________________________________________________________________________________________');
disp('LO CORRECTION 1');
% Calculation of LO for device 1
[PPM_1, PPM2_1] = ltess(signal1_ltess, sampling_rate_ltess);
fprintf('- Device 1 --> PPM: %f - PPM2: %f\n', PPM_1, PPM2_1);

disp('______________________________________________________________________________________________');
disp('LO CORRECTION 2');

% Calculation of LO for device 2
[PPM_2, PPM2_2] = ltess(signal2_ltess, sampling_rate_ltess);
fprintf('- Device 2 --> PPM: %f - PPM2: %f\n', PPM_2, PPM2_2);

disp('______________________________________________________________________________________________');
disp('LO CORRECTION 3');

% Calculation of LO for device 1
[PPM_3, PPM2_3] = ltess(signal3_ltess, sampling_rate_ltess);
fprintf('- Device 3 --> PPM: %f - PPM2: %f\n', PPM_1, PPM2_1);

%% Calculate TDOA
fd1 = fopen(filedev1_tdoa,'r');
fd2 = fopen(filedev2_tdoa,'r');
fd3 = fopen(filedev3_tdoa,'r');

doa_meters_12 = zeros(100,1);
doa_meters_13 = zeros(100,1);
doa_meters_23 = zeros(100,1);
doa_meters_corr_12 = zeros(100,1);
doa_meters_corr_13 = zeros(100,1);
doa_meters_corr_23 = zeros(100,1);

doa_samples_12 = zeros(100,1);
doa_samples_13 = zeros(100,1);
doa_samples_23 = zeros(100,1);
doa_samples_corr_12 = zeros(100,1);
doa_samples_corr_13 = zeros(100,1);
doa_samples_corr_23 = zeros(100,1);

for ii = 1:100
    % We need to read chunk by chunk, otherwise machine will crash
    disp('______________________________________________________________________________________________');
    disp('READ TDOA DATA FROM FILES');
    
    % Signal 1
    signal1 = chunk_load(fd1, num_samples_per_freq);
    disp('File 1 Successfully loaded');
    
    % Signal 2
    signal2 = chunk_load(fd2, num_samples_per_freq);
    disp('File 2 Successfully loaded');
    
    % Signal 3
    signal3 = chunk_load(fd3, num_samples_per_freq);
    disp('File 3 Successfully loaded');
    
    % Correct signals
    signal1_corr = correct_chunk(signal1, PPM_1, sampling_rate_tdoa, fRS);
    signal2_corr = correct_chunk(signal2, PPM_2, sampling_rate_tdoa, fRS);
    signal3_corr = correct_chunk(signal3, PPM_3, sampling_rate_tdoa, fRS);
    
    disp(' ');
    disp('______________________________________________________________________________________________');
    disp('CORRELATION 1 & 2');
    [doa_meters_12(ii), doa_samples_12(ii)] = tdoa(signal1, signal2, sampling_rate_tdoa, num_samples_per_slice, smoothing_factor, corr_type);
    [doa_meters_corr_12(ii), doa_samples_corr_12(ii)] = tdoa(signal1_corr, signal2_corr, sampling_rate_tdoa, num_samples_per_slice, smoothing_factor, corr_type);

    disp(' ');
    disp('______________________________________________________________________________________________');
    disp('CORRELATION 1 & 3');
    [doa_meters_13(ii), doa_samples_13(ii)] = tdoa(signal1, signal3, sampling_rate_tdoa, num_samples_per_slice, smoothing_factor, corr_type);
    [doa_meters_corr_13(ii), doa_samples_corr_13(ii)] = tdoa(signal1_corr, signal3_corr, sampling_rate_tdoa, num_samples_per_slice, smoothing_factor, corr_type);

    disp(' ');
    disp('______________________________________________________________________________________________');
    disp('CORRELATION 2 & 3');
    [doa_meters_23(ii), doa_samples_23(ii)] = tdoa(signal2, signal3, sampling_rate_tdoa, num_samples_per_slice, smoothing_factor, corr_type);
    [doa_meters_corr_23(ii), doa_samples_corr_23(ii)] = tdoa(signal2_corr, signal3_corr, sampling_rate_tdoa, num_samples_per_slice, smoothing_factor, corr_type);
end

fclose(fd1); fclose(fd2); fclose(fd3);

save(['results/', results_filename], 'doa_samples', 'doa_meters', 'doa_samples_2', 'doa_meters_2', 'correlation_value', 'correlation_value_interp');