%% Config File for TDOA setup

% RX and Ref TX Position
rx1_lat = 0; % RX 1
rx1_long = 0;

rx2_lat = 0; % RX 2
rx2_long = 0;

rx3_lat = 0; % RX 3
rx3_long = 0;

tx_ref_lat = -10; % Referenz: Rotenberg DAB
tx_ref_long = 7.77116;

% Variables for composing
fRS_MHz = 100;
fUS_MHz = 806;

lo_correction = false;
interpol_factor = 10;
corr_type = 'dphase';

% Intervals for data
sampling_rate_ltess = 1.92e6;
sampling_rate_tdoa = 2e6;

% signal processing parameters
signal_bandwidth_khz = 0;  % 400, 200, 40, 12, 0(no)
smoothing_factor = 0;

%% Variable formation
% IQ Data Files
folder_identifier_ltess = ['/home/yago/Imdea/git/localization/localization-scripts/data/ltess/', num2str(fRS_MHz), '_', num2str(fUS_MHz), '/'];
folder_identifier_tdoa = ['/home/yago/Imdea/git/localization/localization-scripts/data/localization/', num2str(fRS_MHz), '_', num2str(fUS_MHz), '/'];

% Frequencies for Reference and Unknown signals
fRS = fRS_MHz * 1e6;
fUS = fUS_MHz * 1e6;

% filename for results
if lo_correction
    if interpol_factor < 2
        results_filename = [num2str(round(fRS_MHz)), num2str(fUS_MHz), '_fo_correction_no_interp_', corr_type];
    elseif interpol_factor < 10
        results_filename = [num2str(round(fRS_MHz)), num2str(fUS_MHz), '_fo_correction_0', num2str(interpol_factor), '_interp_', corr_type];
    else
        results_filename = [num2str(round(fRS_MHz)), num2str(fUS_MHz), '_fo_correction_', num2str(interpol_factor), '_interp_', corr_type];
    end
else
    results_filename = [num2str(round(fRS_MHz)), num2str(fUS_MHz), '_original_', corr_type];
end

% 0: no plots
% 1: show correlation plots
% 2: show also input spcetrograms and spectra of input meas
% 3: show also before and after filtering
report_level = 1;

% For parallel processing
M = 4;