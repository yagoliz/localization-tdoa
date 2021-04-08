% =========================================================================
%  Real test realization
%  Author: Yago Lizarribar
% =========================================================================
clear; clc; close all; warning ('off','all');

% adds subfolder with functions to PATH
[p,n,e] = fileparts(mfilename('fullpath'));
addpath([p '/tdoa']);
addpath([p '/ltess']);

%% Load configuration
filename = 'real_tests/locations.json';
problem_data = json_load(filename);

% Unknown transmitter
fUS = problem_data.transmitters.unknown.freq * 1e6;
sr_local = 2e6;
us_lat = problem_data.transmitters.unknown.coord(1);
us_lon = problem_data.transmitters.unknown.coord(2);

% Reference transmitter
fRS = problem_data.transmitters.reference.freq * 1e6;
sr_ltess = 1.92e6;
rs_lat = problem_data.transmitters.reference.coord(1);
rs_lon = problem_data.transmitters.reference.coord(2);

% Sensors
NUM_SENSORS = length(problem_data.sensors);
NUM_FILES_SENSOR = 5;
EXPERIMENT = 5;

folder_location = '/home/ygglc/Imdea/git/localization/localization-tdoa/data/real_tests';
localization_files = cell(NUM_SENSORS, NUM_FILES_SENSOR);
ltess_files = cell(NUM_SENSORS, NUM_FILES_SENSOR);

sensors = zeros(NUM_SENSORS, 2);
distances_rs = zeros(NUM_SENSORS,1);
distances_us = distances_rs;

for ii = 1:NUM_SENSORS
    sensors(ii,:) = problem_data.sensors(ii).coordinates;
    distances_rs(ii) = 1000*deg2km(distance(rs_lat, rs_lon, sensors(ii,1), sensors(ii,2)));
    distances_us(ii) = 1000*deg2km(distance(us_lat, us_lon, sensors(ii,1), sensors(ii,2)));
    
    localization_struct = dir([folder_location, filesep, problem_data.sensors(ii).name, filesep, '*localization*']);
    ltess_struct = dir([folder_location, filesep, problem_data.sensors(ii).name, filesep, '*ltess*']);
    for jj = 1:NUM_FILES_SENSOR
        localization_files{ii,jj} = [localization_struct(jj).folder, filesep, localization_struct(jj).name];
        ltess_files{ii,jj} = [ltess_struct(jj).folder, filesep, ltess_struct(jj).name];
    end
end

%% PPM correction
signal_local = cell(NUM_SENSORS,1);
signal_ltess = cell(NUM_SENSORS,1);
tol = 1e-2;
for ii = 1:NUM_SENSORS
    signal_sensor = spec_load(localization_files{ii,EXPERIMENT});
    signal_ltess{ii} = spec_load(ltess_files{ii,EXPERIMENT});

    PPM = inf; count = 0; MAX = 10;
    while abs(PPM) > tol && count < MAX
        [PPM, PPM2] = ltess(signal_ltess{ii}, sr_ltess);
        signal_ltess{ii} = correct_fo_ltess(signal_ltess{ii},PPM,sr_ltess, 806e6);
        signal_sensor = correct_fo(signal_sensor, PPM, sr_local, fRS, fUS);
        
        count = count + 1;
    end
    signal_local{ii} = signal_sensor;
end

%% TDOA estimation
combinations = nchoosek(1:NUM_SENSORS,2);
N = size(combinations,1);

% Signal parameters
max_lag = round(100*1000*sr_local/3e8);
corr_type = 'dphase';
report_level = 1;
signal_bandwidth_khz = 0;
interpol_factor = 1;
ns_freq = 2e6;
ns_slice = 0.8 * ns_freq;

doa_meters = zeros(N,1);
doa_samples = zeros(N,1);
us_doa = zeros(N,1);
keep = zeros(N,1);
for ii = 1:N
    % Select the sensors and compute distance differences
    si = combinations(ii,1);
    sj = combinations(ii,2);
    rs_diff_ij = distances_rs(si) - distances_rs(sj);
    us_doa(ii) = distances_us(si) - distances_us(sj);
    
    % Compute TDOA and store it
    [doa_meters(ii), doa_samples(ii), ~, ~, ~, ~, keep(ii)] = ...
        tdoa2(signal_local{si}, signal_local{sj}, ns_freq, ns_slice, sr_local, rs_diff_ij, ...
              max_lag, corr_type, report_level, signal_bandwidth_khz, interpol_factor);
end

%% Multilateration
keep = keep > 0.5;
X0 = mean(sensors,1);
options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt');

optimfun = @(X) optimnlls(X,doa_meters,sensors,combinations,keep);
X = lsqnonlin(optimfun,X0,[],[],options);
X

%% Non linear LeastSquares function
function F = optimnlls(X,doa,sensors,combinations,keep)
    lat = X(1);
    lon = X(2);
    
    d = zeros(size(sensors,1),1);
    for ii = 1:size(sensors,1)
        d(ii) = deg2km(distance(lat, lon, sensors(ii,1), sensors(ii,2)));
    end
    
    t = zeros(size(combinations,1),1);
    for ii = 1:size(t,1)
        si = combinations(ii,1);
        sj = combinations(ii,2);
        t(ii) = (d(si) - d(sj))*1000;
    end
    
    F = (doa(keep) - t(keep));
end
