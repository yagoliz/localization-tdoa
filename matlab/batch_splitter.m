% =========================================================================
%  Batch test realization
%  Author: Yago Lizarribar
% =========================================================================

clear;
clc;
close all;
warning ('off','all');

% adds subfolder with functions to PATH
[p,n,e] = fileparts(mfilename('fullpath'));
addpath([p '/tdoa']);
addpath([p '/ltess']);

global c;
c = 299792458; % m/s

%% Read Parameters from config file
%---------------------------------------------
batch_config_splitter;

% Check interpol parameter (needs to be at least 1)
if interpol_factor <= 0
    interpol_factor = 1;
end

% Get all files
filesdev1_ltess = dir([folder_identifier, 's1/', '*ltess*']);
filesdev2_ltess = dir([folder_identifier, 's2/', '*ltess*']);

filesdev1_tdoa = dir([folder_identifier, 's1/', '*localization*']);
filesdev2_tdoa = dir([folder_identifier, 's2/', '*localization*']);

% Create state variables
doa_samples = zeros(length(filesdev1_tdoa),1);
doa_meters  = zeros(length(filesdev1_tdoa),1);
doa_samples_2 = zeros(length(filesdev1_tdoa),1);
doa_meters_2 = zeros(length(filesdev1_tdoa),1);
correlation_value = zeros(length(filesdev1_tdoa),3);
correlation_value_interp = zeros(length(filesdev1_tdoa),3);

%% Main loop
if compute_ppm
    ppm1 = zeros(length(filesdev1_tdoa),1);
    ppm2 = ppm1;
    ppm3 = ppm2;
else
    try
        load(fileppm);
    catch
        disp('WARNING: No PPM file');
        compute_ppm = true; % We'll have to compute the drift
    end
end
PPM_prev_1 = 0;
PPM_prev_2 = 0;
for kk = 38:length(num_samples_per_slice)
samples_slice = num_samples_per_slice(kk);
    for filenum=1:length(filesdev1_tdoa)
        % Read TDOA Signals from File
        file1_tdoa = filesdev1_tdoa(filenum).name;
        file2_tdoa = filesdev2_tdoa(filenum).name;
    %     file3_tdoa = filesdev3_tdoa(filenum).name;
        
        disp('______________________________________________________________________________________________');
        disp('READ TDOA DATA FROM FILES');
        signal1 = spec_load([folder_identifier, 's1/', file1_tdoa]);
        disp('File 1 Successfully loaded');
        signal2 = spec_load([folder_identifier, 's2/', file2_tdoa]);
        disp('File 2 Successfully loaded');
    %     signal3 = spec_load([folder_identifier, 's3/', file3_tdoa]);
    %     disp('File 3 Successfully loaded');
    
        %% Local Oscillator Offset correction
        if lo_correction && compute_ppm
            file1 = filesdev1_ltess(filenum).name;
            file2 = filesdev2_ltess(filenum).name;
    %         file3 = filesdev3_ltess(filenum).name;
    
            % Read Signals from File
            disp('______________________________________________________________________________________________');
            disp('READ LTESS DATA FROM FILES');
            signal1_ltess = spec_load([folder_identifier, 's1/', file1]);
            disp('File 1 Successfully loaded');
            signal2_ltess = spec_load([folder_identifier, 's2/', file2]);
            disp('File 2 Successfully loaded');
    
            disp(' ');
            disp('______________________________________________________________________________________________');
            disp('LO CORRECTION 1');
            % Calculation of LO for device 1
            [PPM, PPM2] = ltess(signal1_ltess, sampling_rate_ltess);
            fprintf('- Device 1 --> PPM: %f - PPM2: %f\n', PPM, PPM2);
            if isnan(PPM)
                PPM = PPM_prev_1;
                ppm1(filenum) = PPM_prev_1;
            else
                PPM_prev_1 = PPM;
                ppm1(filenum) = PPM;
            end
            % Correction of FO and sampling rate
            signal1 = correct_fo(signal1, PPM, sampling_rate_tdoa, fRS, fUS);
    
            disp('______________________________________________________________________________________________');
            disp('LO CORRECTION 2');
            
            % Calculation of LO for device 2
            [PPM, PPM2] = ltess(signal2_ltess, sampling_rate_ltess);
            fprintf('- Device 2 --> PPM: %f - PPM2: %f\n', PPM, PPM2);
            if isnan(PPM)
                PPM = PPM_prev_2;
                ppm2(filenum) = PPM_prev_2;
            else
                PPM_prev_2 = PPM;
                ppm2(filenum) = PPM;
            end
            % Correction of FO and sampling rate
            signal2 = correct_fo(signal2, PPM, sampling_rate_tdoa, fRS, fUS);
    
        elseif lo_correction && ~compute_ppm
            signal1 = correct_fo(signal1, ppm1(filenum), sampling_rate_tdoa, fRS, fUS);
            signal2 = correct_fo(signal2, ppm2(filenum), sampling_rate_tdoa, fRS, fUS);        
        end
    
        %% Calculate TDOA
            disp(' ');
            disp('______________________________________________________________________________________________');
            disp('CORRELATION 1 & 2');
            [doa_meters(filenum), doa_samples(filenum), doa_meters_2(filenum), doa_samples_2(filenum), correlation_value(filenum,:), correlation_value_interp(filenum,:)] = ...
                tdoa2(signal1, signal2, num_samples_per_freq, samples_slice, sampling_rate_tdoa, 0, smoothing_factor, corr_type, report_level, signal_bandwidth_khz, interpol_factor, false, 0);
        
    end
    results_file = [results_filename,'_',num2str(samples_slice)];
    save(['convergence/', results_file], 'doa_samples', 'doa_meters', 'doa_samples_2', 'doa_meters_2', 'correlation_value', 'correlation_value_interp');
end

if lo_correction && compute_ppm
    save(['ppm/',fileppm], 'ppm1', 'ppm2');
end