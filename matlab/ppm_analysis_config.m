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
f_MHz = 806;

corr_type = 'dphase';

% Intervals for data
num_samples_per_freq = 2e6; 
num_samples_per_slice = 0.8 * num_samples_per_freq;
sampling_rate_ltess = 1.92e6;
sampling_rate_tdoa = 2e6;

% signal processing parameters
signal_bandwidth_khz = 0;  % 400, 200, 40, 12, 0(no)
smoothing_factor = 0;

%% Variable formation
% IQ Data Files
folder_identifier = '/home/yago/Imdea/git/localization/localization-tdoa/data/ppm_stability/';


% Frequencies for Reference and Unknown signals
fRS = f_MHz * 1e6;
