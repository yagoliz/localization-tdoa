% =========================================================================
%  Experimental Evaluation Script for RTL-SDR based TDOA
%  DC9ST, 2017-2019
% =========================================================================

clear;
clc;
close all;

% adds subfolder with functions to PATH
[p,n,e] = fileparts(mfilename('fullpath'));
addpath([p '/tdoa']);
addpath([p '/ltess']);
addpath([p '/test']); % only required for the test setups


%% Read Parameters from config qfile, that specifies all parameters
%---------------------------------------------
test_config;

% create filenames
dateiname1 = [folder_identifier, 'D0-E4-ltess.dat'];
dateiname2 = [folder_identifier, 'D1-E4-ltess.dat'];
dateiname3 = [folder_identifier, 'D2-E4-ltess.dat'];

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

%% Read Signals from File
disp('______________________________________________________________________________________________');
disp('READ DATA FROM FILES');
signal1 = read_file_iq(dateiname1);
disp('File 1 Successfully loaded');
signal2 = read_file_iq(dateiname2);
disp('File 2 Successfully loaded');
signal3 = read_file_iq(dateiname3);
disp('File 3 Successfully loaded');

if (report_level > 1)
    % display raw signals
    num_samples_total = length(signal1);
    inphase1 = real(signal1);
    quadrature1 = imag(signal1);
    inphase2 = real(signal2);
    quadrature2 = imag(signal2);
    inphase3 = real(signal3);
    quadrature3 = imag(signal3);
    
    figure;
    subplot(3,1,1);
%     plot(1:num_samples_total, inphase1(1:num_samples_total), 1:num_samples_total, quadrature1(1:num_samples_total));
    plot(1:num_samples_total, abs(signal1));
    title('raw RX 1: I and Q');
    subplot(3,1,2);
%     plot(1:num_samples_total, inphase2(1:num_samples_total), 1:num_samples_total, quadrature2(1:num_samples_total));
    plot(1:num_samples_total, abs(signal2));
    title('raw RX 2: I and Q');    
    subplot(3,1,3);
%     plot(1:num_samples_total, inphase3(1:num_samples_total), 1:num_samples_total, quadrature3(1:num_samples_total));
    plot(1:num_samples_total, abs(signal3));
    title('raw RX 3: I and Q');
end

if (report_level > 1)
    % calculate and show spectrogram
    nfft = 256;
    overlap = -1000;
    
    figure;
    subplot(4,2,1);
    complex_signal = detrend(signal1);
    [~,F,T,P] = spectrogram(complex_signal, nfft, overlap, nfft, 2e6 );
    spectrum = fftshift(fliplr(10*log10(abs(P))'), 2);
    for i=1:nfft
        spectrum(:,i) = smooth(spectrum(:,i),9);
    end
    surf(T,F, spectrum', 'edgecolor', 'none');
    axis tight;
    view(0,90);
    title('RX 1');
    xlabel('time');
    ylabel('frequency');

    subplot(4,2,3);
    complex_signal = detrend(signal2);
    [~,F,T,P] = spectrogram(complex_signal, nfft, overlap, nfft, 2e6 );
    spectrum = fftshift(fliplr(10*log10(abs(P))'), 2);
    for i=1:nfft
        spectrum(:,i) = smooth(spectrum(:,i),9);
    end
    surf(T,F, spectrum', 'edgecolor', 'none');
    axis tight;
    view(0,90);
    title('RX 2');
    xlabel('time');
    ylabel('frequency');

    
    subplot(4,2,5);
    complex_signal = detrend(signal3);
    [S,F,T,P] = spectrogram(complex_signal, nfft, overlap, nfft, 2e6 );
    spectrum = fftshift(fliplr(10*log10(abs(P))'), 2);
    for i=1:nfft
        spectrum(:,i) = smooth(spectrum(:,i),9);
    end
    surf(T,F, spectrum', 'edgecolor', 'none');
    axis tight;
    view(0,90);
    title('RX 3');
    xlabel('time');
    ylabel('frequency');
    
    % display spectrum
    spectrum_smooth_factor  = 201; 
    subplot(4,2,2);
    spectrum_single1 = 10*log10(abs(fftshift(fft(signal1(1.7e6 : 1.7e6 + 2^18)))));
    spectrum_single1 = smooth(spectrum_single1, spectrum_smooth_factor);
    plot(spectrum_single1);
    title('Measurement RX 1');
    grid;

    subplot(4,2,4);
    spectrum_single2 = 10*log10(abs(fftshift(fft(signal2(1.7e6 : 1.7e6 + 2^18)))));
    spectrum_single2 = smooth(spectrum_single2, spectrum_smooth_factor);
    plot(spectrum_single2);
    title('Measurement RX 2');
    grid;

    subplot(4,2,6);
    spectrum_single3 = 10*log10(abs(fftshift(fft(signal3(1.7e6 : 1.7e6 + 2^18)))));
    spectrum_single3 = smooth(spectrum_single3, spectrum_smooth_factor);
    plot(spectrum_single3);
    title('Measurement RX 3');
    grid;

    subplot(4,2,7:8);
    freq_axis = -(length(spectrum_single1)/2) : 1 : ((length(spectrum_single1)/2)-1);
    plot(freq_axis, spectrum_single1, freq_axis, spectrum_single2, freq_axis, spectrum_single3);
    title('Measurement Signal RX 1,2 & 3');
    grid;
end

%% Chunk division
signal1_ltess = signal1(1:ltess_interval);
signal2_ltess = signal2(1:ltess_interval);
signal3_ltess = signal3(1:ltess_interval);

signal1 = signal1(ltess_interval + 1 : end);
signal2 = signal2(ltess_interval + 1 : end);
signal3 = signal3(ltess_interval + 1 : end);

%% Local Oscillator Offset correction
if lo_correction
    disp(' ');
    disp('______________________________________________________________________________________________');
    disp('LO CORRECTION 1');
    % Calculation of LO for device 1
    [PPM, PPM2] = ltess(signal1_ltess, sampling_rate_ltess);
    fprintf('- Device 1 --> PPM: %f - PPM2: %f\n', PPM, PPM2);
    % Correction of FO and sampling rate
    signal1 = correct_fo(signal1, PPM, sampling_rate_tdoa, fRS, fUS);

    disp('______________________________________________________________________________________________');
    disp('LO CORRECTION 2');
    % Calculation of LO for device 2
    [PPM, PPM2] = ltess(signal2_ltess, sampling_rate_ltess);
    fprintf('- Device 2 --> PPM: %f - PPM2: %f\n', PPM, PPM2);
    % Correction of FO and sampling rate
    signal2 = correct_fo(signal2, PPM, sampling_rate_tdoa, fRS, fUS);

    disp('______________________________________________________________________________________________');
    disp('LO CORRECTION 3');
    % Calculation of LO for device 1
    [PPM, PPM2] = ltess(signal3_ltess, sampling_rate_ltess);
    fprintf('- Device 3 --> PPM: %f - PPM2: %f\n', PPM, PPM2);
    % Correction of FO and sampling rate
    signal3 = correct_fo(signal3, PPM, sampling_rate_tdoa, fRS, fUS);
end

%% Calculate TDOA

disp(' ');
disp('______________________________________________________________________________________________');
disp('CORRELATION 1 & 2');
[doa_meters12, doa_samples12, reliability12 ] = tdoa2(signal1, signal2, rx_distance_diff12, rx_distance12, smoothing_factor, corr_type, report_level, signal_bandwidth_khz, interpol_factor);

disp(' ');
disp('______________________________________________________________________________________________');
disp('CORRELATION 1 & 3');
[doa_meters13, doa_samples13, reliability13 ] = tdoa2(signal1, signal3, rx_distance_diff13, rx_distance13, smoothing_factor, corr_type, report_level, signal_bandwidth_khz, interpol_factor);

disp(' ');
disp('______________________________________________________________________________________________');
disp('CORRELATION 2 & 3');
[doa_meters23, doa_samples23, reliability23 ] = tdoa2(signal2, signal3, rx_distance_diff23, rx_distance23, smoothing_factor, corr_type, report_level, signal_bandwidth_khz, interpol_factor);


%% Generate html map
if generate_map
    disp(' ');
    disp('______________________________________________________________________________________________');
    disp('GENERATE HYPERBOLAS');

    [points_lat1, points_long1] = gen_hyperbola(doa_meters12, rx1_lat, rx1_long, rx2_lat, rx2_long, geo_ref_lat, geo_ref_long);
    [points_lat2, points_long2] = gen_hyperbola(doa_meters13, rx1_lat, rx1_long, rx3_lat, rx3_long, geo_ref_lat, geo_ref_long);
    [points_lat3, points_long3] = gen_hyperbola(doa_meters23, rx2_lat, rx2_long, rx3_lat, rx3_long, geo_ref_lat, geo_ref_long);

    disp(' ');
    disp('______________________________________________________________________________________________');
    disp('GENERATE HTML');
    rx_lat_positions  = [rx1_lat   rx2_lat   rx3_lat ];
    rx_long_positions = [rx1_long  rx2_long  rx3_long];

    hyperbola_lat_cell  = {points_lat1,  points_lat2, points_lat3};
    hyperbola_long_cell = {points_long1, points_long2, points_long3};

    [heatmap_long, heatmap_lat, heatmap_mag] = create_heatmap(doa_meters12, doa_meters13, doa_meters23, rx1_lat, rx1_long, rx2_lat, rx2_long, rx3_lat, rx3_long, heatmap_resolution, geo_ref_lat, geo_ref_long); % generate heatmap
    heatmap_cell = {heatmap_long, heatmap_lat, heatmap_mag};

    if strcmp(map_mode, 'google_maps')
        % for google maps
        create_html_file_gm( ['ergebnisse/map_' file_identifier '_' corr_type '_interp' num2str(interpol_factor) '_bw' int2str(signal_bandwidth_khz) '_smooth' int2str(smoothing_factor) '_gm.html'], rx_lat_positions, rx_long_positions, hyperbola_lat_cell, hyperbola_long_cell, heatmap_cell, heatmap_threshold);
    else
        % for open street map
        create_html_file_osm( ['ergebnisse/map_' file_identifier '_' corr_type '_interp' num2str(interpol_factor) '_bw' int2str(signal_bandwidth_khz) '_smooth' int2str(smoothing_factor) '_osm.html'], rx_lat_positions, rx_long_positions, hyperbola_lat_cell, hyperbola_long_cell, heatmap_cell, heatmap_threshold);
    end
    disp('______________________________________________________________________________________________');
end