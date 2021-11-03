%% Config File for TDOA setup

% Variables for composing
tx_type = 'dvbt';
if strcmp(tx_type, 'dvbt')
    fRS_MHz = 627;
elseif strcmp(tx_type, 'dab')
    fRS_MHz = 196;
elseif strcmp(tx_type, 'gsm')
    fRS_MHz = 932;
end

us_type = 'gsm';
if strcmp(us_type, 'dvbt')
    fUS_MHz = 627;
elseif strcmp(us_type, 'dab')
    fUS_MHz = 196;
elseif strcmp(us_type, 'gsm')
    fUS_MHz = 932;
elseif strcmp(us_type, 'lte')
    fUS_MHz = 806;
end

lo_correction = true;
compute_ppm = false;
fileppm = ['ppm/', 'ppm_rs_', tx_type, '_us_', us_type];
interpol_factor = 10;
corr_type = 'dphase';

% Intervals for data
num_samples_per_freq = 1e6; 
num_samples_per_slice = 0.7 * num_samples_per_freq;
sampling_rate_ltess = 1.92e6;
sampling_rate_tdoa = 2e6;

% signal processing parameters
signal_bandwidth_khz = 0;  % 400, 200, 40, 12, 0(no)
smoothing_factor = 0;

%% Variable formation
% IQ Data Files
folder_identifier = ['/home/ygglc/desktop/phd/research/localization/core/localization-tdoa/data/splitter/', tx_type, '/', us_type, '/'];

% Frequencies for Reference and Unknown signals
fRS = fRS_MHz * 1e6;
fUS = fUS_MHz * 1e6;

% filename for results
if lo_correction
    if interpol_factor < 2
        results_filename = ['splitter_', tx_type, '_', us_type, '_fo_correction_no_interp_', corr_type];
    elseif interpol_factor < 10
        results_filename = ['splitter_', tx_type, '_', us_type, '_fo_correction_0', num2str(interpol_factor), '_interp_', corr_type];
    else
        results_filename = ['splitter_', tx_type, '_', us_type, '_fo_correction_', num2str(interpol_factor), '_interp_', corr_type];
    end
else
    results_filename = ['splitter_', tx_type, '_', us_type, '_original_', corr_type];
end

% 0: no plots
% 1: show correlation plots
% 2: show also input spcetrograms and spectra of input meas
% 3: show also before and after filtering
report_level = 1;
