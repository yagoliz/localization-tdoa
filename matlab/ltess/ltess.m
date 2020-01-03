%% Calculate the Local Oscillator (LO) offset for a given chunk of data
function [PPM, PPM2] = ltess(chunk, sampling_rate)
% Configuration parameters
RESAMPLE_FACTOR = 60;
PSS_STEP = 9600;
SEARCH_WINDOW = 150;
CORRELATION_FACTOR = 0.1;
CORRELATION_FACTOR_PREAMBLE = 0.1;
PREAMBLE=20; % Number of analyzed PSS before start to jump
FLIP=0; % 0 disabled 1: enabled  (testing purposes)
POLYNOMIAL_DEGREE=1;

[Z, Z_t] = get_Zadoof();

% First run
[PPM,PSS_percent,~,~,~,PPM2] = getDrift(chunk,sampling_rate, Z, Z_t, PREAMBLE, ...
            PSS_STEP, SEARCH_WINDOW, CORRELATION_FACTOR, RESAMPLE_FACTOR, FLIP, POLYNOMIAL_DEGREE);
fprintf('Zadoof adaptation --> PPM: %f [%f] - PSS detected: %f\n', PPM, PPM2, PSS_percent)

% Second run
T_s = 1/sampling_rate;
delta_f=(PPM*1e-6)*806e6;
Z_t_rotated = {};
Z_t_rotated{1} = Z_t{1}.*exp(-1i*2*pi*T_s*delta_f*(1:length(Z_t{1}))');
Z_t_rotated{2} = Z_t{2}.*exp(-1i*2*pi*T_s*delta_f*(1:length(Z_t{2}))');
Z_t_rotated{3} = Z_t{3}.*exp(-1i*2*pi*T_s*delta_f*(1:length(Z_t{3}))');

[PPM,~,~,~,~,PPM2, ~, ~, ~, ~ ] = getDrift(chunk,sampling_rate, Z, Z_t_rotated, ...
    PREAMBLE, PSS_STEP, SEARCH_WINDOW, CORRELATION_FACTOR, RESAMPLE_FACTOR, FLIP, POLYNOMIAL_DEGREE);

fprintf('Final result --> PPM: %f [%f] - PSS detected: %f\n', PPM, PPM2, PSS_percent)
end