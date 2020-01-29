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
fRS_MHz = 806;
fUS_MHz = fRS_MHz;

lo_correction = true;
interpol_factor = 1;
corr_type = 'dphase';

% Intervals for data
num_samples_per_freq = 1.92e6;
num_samples_per_slice = 1.2e6;
sampling_rate = 1.92e6;

% signal processing parameters
signal_bandwidth_khz = 0;  % 400, 200, 40, 12, 0(no)
smoothing_factor = 0;

%% Variable formation
% IQ Data Files
folder_identifier = '/home/yago/Imdea/git/localization/localization-tdoa/data/ltess/20191122/';

% Frequencies for Reference and Unknown signals
fRS = fRS_MHz * 1e6;
fUS = fUS_MHz * 1e6;

% filename for results
if lo_correction
    if interpol_factor < 2
        results_filename = ['lte_fo_correction_no_interp_', corr_type];
    elseif interpol_factor < 10
        results_filename = ['lte_fo_correction_0', num2str(interpol_factor), '_interp_', corr_type];
    else
        results_filename = ['lte_fo_correction_', num2str(interpol_factor), '_interp_', corr_type];
    end
else
    results_filename = ['lte_original_', corr_type];
end

% 0: no plots
% 1: show correlation plots
% 2: show also input spcetrograms and spectra of input meas
% 3: show also before and after filtering
report_level = 1;

% For parallel processing
M = 4;
