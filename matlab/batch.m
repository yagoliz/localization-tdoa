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
batch_config;

% Check interpol parameter (needs to be at least 1)
if interpol_factor <= 0
    interpol_factor = 1;
end

% Get all files
filesdev1_ltess = dir([folder_identifier_ltess, 'D0*.dat']);
filesdev2_ltess = dir([folder_identifier_ltess, 'D1*.dat']);
filesdev3_ltess = dir([folder_identifier_ltess, 'D2*.dat']);

filesdev1_tdoa = dir([folder_identifier_tdoa, 'D0*.dat']);
filesdev2_tdoa = dir([folder_identifier_tdoa, 'D1*.dat']);
filesdev3_tdoa = dir([folder_identifier_tdoa, 'D2*.dat']);

% calculate geodetic reference point as mean center of all RX positions
geo_ref_lat  = mean([rx1_lat, rx2_lat, rx3_lat]);
geo_ref_long = mean([rx1_long, rx2_long, rx3_long]);
disp(['geodetic reference point (mean of RX positions): lat=' num2str(geo_ref_lat, 8) ', long=' num2str(geo_ref_long, 8) ])

% known signal path differences between two RXes to Ref (sign of result is important!)
rx_distance_diff12 = dist_latlong(tx_ref_lat, tx_ref_long, rx1_lat, rx1_long, geo_ref_lat, geo_ref_long) - dist_latlong(tx_ref_lat, tx_ref_long, rx2_lat, rx2_long, geo_ref_lat, geo_ref_long); % (Ref to RX1 - Ref to RX2) in meters
rx_distance_diff13 = dist_latlong(tx_ref_lat, tx_ref_long, rx1_lat, rx1_long, geo_ref_lat, geo_ref_long) - dist_latlong(tx_ref_lat, tx_ref_long, rx3_lat, rx3_long, geo_ref_lat, geo_ref_long); % (Ref to RX1 - Ref to RX3) in meters
rx_distance_diff23 = dist_latlong(tx_ref_lat, tx_ref_long, rx2_lat, rx2_long, geo_ref_lat, geo_ref_long) - dist_latlong(tx_ref_lat, tx_ref_long, rx3_lat, rx3_long, geo_ref_lat, geo_ref_long); % (Ref to RX2 - Ref to RX3) in meters

% distance between two RXes in meters
rx_distance12 = dist_latlong(rx1_lat, rx1_long, rx2_lat, rx2_long, geo_ref_lat, geo_ref_long);
rx_distance13 = dist_latlong(rx1_lat, rx1_long, rx3_lat, rx3_long, geo_ref_lat, geo_ref_long);
rx_distance23 = dist_latlong(rx2_lat, rx2_long, rx3_lat, rx3_long, geo_ref_lat, geo_ref_long);

% Create state variables
doa_samples = zeros(length(filesdev1_tdoa)*3,1);
doa_meters  = zeros(length(filesdev1_tdoa)*3,1);
doa_samples_2 = zeros(length(filesdev1_tdoa)*3,1);
doa_meters_2 = zeros(length(filesdev1_tdoa)*3,1);
correlation_value = zeros(3,3,length(filesdev1_tdoa)*3);
correlation_value_interp = zeros(3,3,length(filesdev1_tdoa)*3);

%% Main loop
PPM_prev_1 = 0;
PPM_prev_2 = 0;
PPM_prev_3 = 0;
for filenum=1:length(filesdev1_tdoa)
    % Read TDOA Signals from File
    file1_tdoa = filesdev1_tdoa(filenum).name;
    file2_tdoa = filesdev2_tdoa(filenum).name;
    file3_tdoa = filesdev3_tdoa(filenum).name;
    
    disp('______________________________________________________________________________________________');
    disp('READ TDOA DATA FROM FILES');
    signal1 = spec_load([folder_identifier_tdoa, file1_tdoa]);
    disp('File 1 Successfully loaded');
    signal2 = spec_load([folder_identifier_tdoa, file2_tdoa]);
    disp('File 2 Successfully loaded');
    signal3 = spec_load([folder_identifier_tdoa, file3_tdoa]);
    disp('File 3 Successfully loaded');

    %% Local Oscillator Offset correction
    if lo_correction
        file1 = filesdev1_ltess(filenum).name;
        file2 = filesdev2_ltess(filenum).name;
        file3 = filesdev3_ltess(filenum).name;

        % Read Signals from File
        disp('______________________________________________________________________________________________');
        disp('READ LTESS DATA FROM FILES');
        signal1_ltess = spec_load([folder_identifier_ltess, file1]);
        disp('File 1 Successfully loaded');
        signal2_ltess = spec_load([folder_identifier_ltess, file2]);
        disp('File 2 Successfully loaded');
        signal3_ltess = spec_load([folder_identifier_ltess, file3]);
        disp('File 3 Successfully loaded');

        disp(' ');
        disp('______________________________________________________________________________________________');
        disp('LO CORRECTION 1');
        % Calculation of LO for device 1
        [PPM, PPM2] = ltess(signal1_ltess, sampling_rate_ltess);
        fprintf('- Device 1 --> PPM: %f - PPM2: %f\n', PPM, PPM2);
        if isnan(PPM)
            PPM = PPM_prev_1;
        else
            PPM_prev_1 = PPM;
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
        else
            PPM_prev_2 = PPM;
        end
        % Correction of FO and sampling rate
        signal2 = correct_fo(signal2, PPM, sampling_rate_tdoa, fRS, fUS);

        disp('______________________________________________________________________________________________');
        disp('LO CORRECTION 3');

        % Calculation of LO for device 1
        [PPM, PPM2] = ltess(signal3_ltess, sampling_rate_ltess);
        if isnan(PPM)
            PPM = PPM_prev_3;
        else
            PPM_prev_3 = PPM;
        end
        fprintf('- Device 3 --> PPM: %f - PPM2: %f\n', PPM, PPM2);
        % Correction of FO and sampling rate
        signal3 = correct_fo(signal3, PPM, sampling_rate_tdoa, fRS, fUS);
    end

    %% Calculate TDOA

    disp(' ');
    disp('______________________________________________________________________________________________');
    disp('CORRELATION 1 & 2');
    [doa_meters(filenum*3 - 2), doa_samples(filenum*3 - 2), doa_meters_2(filenum*3 - 2), doa_samples_2(filenum*3 - 2), correlation_value(1,:,filenum), correlation_value_interp(1,:,filenum)] = tdoa2(signal1, signal2, rx_distance_diff12, rx_distance12, smoothing_factor, corr_type, report_level, signal_bandwidth_khz, interpol_factor);

    disp(' ');
    disp('______________________________________________________________________________________________');
    disp('CORRELATION 1 & 3');
    [doa_meters(filenum*3 - 1), doa_samples(filenum*3 - 1), doa_meters_2(filenum*3 - 1), doa_samples_2(filenum*3 - 1), correlation_value(2,:,filenum), correlation_value_interp(2,:,filenum)] = tdoa2(signal1, signal3, rx_distance_diff13, rx_distance13, smoothing_factor, corr_type, report_level, signal_bandwidth_khz, interpol_factor);

    disp(' ');
    disp('______________________________________________________________________________________________');
    disp('CORRELATION 2 & 3');
    [doa_meters(filenum*3), doa_samples(filenum*3), doa_meters_2(filenum*3), doa_samples_2(filenum*3), correlation_value(3,:,filenum), correlation_value_interp(3,:,filenum)] = tdoa2(signal2, signal3, rx_distance_diff23, rx_distance23, smoothing_factor, corr_type, report_level, signal_bandwidth_khz, interpol_factor);
end

save(['results/', results_filename], 'doa_samples', 'doa_meters', 'doa_samples_2', 'doa_meters_2', 'correlation_value', 'correlation_value_interp');