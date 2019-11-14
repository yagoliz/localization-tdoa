%% Configuration file for the analysis of data for TDOA localization

% General parameters
lo_correction = 1;
tdoa_calculation = 1;
second_correlation = 0; % To do the secondary check

% Receivers and detected frequencies
RECEIVERS = 3;
fRS = 806e6;
fUS = 806e6;

% Sampling rates
SAMPLING_RATE_LTESS = 1.92e6;
SAMPLING_RATE_TDOA = 1.92e6;
GUARD_INTERVAL = 0.1e6;

% Methods to use on the correlation
methods = {'abs', 'dphase'};

% Paths for functions
addpath('common');
addpath('tdoa');

% Path for data
data_path = '/home/yago/Imdea/git/localization/localization-scripts/data/ltess/';
index_to_analyze = 8;

% Where to analyze the signals
index_samples = zeros(4,1);

% Divide unknown signal in smaller chunks
chunks = 5;

% Resampling factor
resampling_factor = 5;
          
          
