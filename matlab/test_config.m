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


% IQ Data Files
file_identifier = 'test.dat';
folder_identifier = '/home/yago/Imdea/git/localization/localization-scripts/data/ltess/';

% Intervals for data
sampling_rate_ltess = 1.92e6;
sampling_rate_tdoa = 1.92e6;
ltess_interval = 2 * sampling_rate_ltess;
tdoa_interval = 3 * sampling_rate_ltess;

% LO correction
lo_correction = true;

% Frequencies for Reference and Unknown signals
fRS = 806e6;
fUS = 806e6;

% signal processing parameters
signal_bandwidth_khz = 0;  % 400, 200, 40, 12, 0(no)
smoothing_factor = 0;
corr_type = 'dphase';  %'abs' or 'dphase'
interpol_factor = 10;

% 0: no plots
% 1: show correlation plots
% 2: show also input spcetrograms and spectra of input meas
% 3: show also before and after filtering
report_level = 1;

% map output
generate_map = false;
% 'open_street_map' (default) or 'google_maps'
map_mode = 'open_street_map';

% heatmap (only with google maps)
heatmap_resolution = 400; % resolution for heatmap points
heatmap_threshold = 0.1;  % heatmap point with lower mag are suppressed for html output