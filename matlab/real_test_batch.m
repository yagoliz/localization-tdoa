% =========================================================================
%  Real test realization
%  Author: Yago Lizarribar
% =========================================================================
clear; clc; close all; warning ('off','all');

% Speed of light
global c err_prev; %m/s
c = 299792458;

% adds subfolder with functions to PATH
[p,n,e] = fileparts(mfilename('fullpath'));
addpath([p '/tdoa']);
addpath([p '/ltess']);
addpath([p, '/geodesy']);
addpath([p, '/optimization']);

%% Load configuration
filename = 'real_tests/locations_dab_dab.json';
problem_data = json_load(filename);

% Unknown transmitter
fUS = problem_data.transmitters.unknown.freq * 1e6;
sr_local = problem_data.config.sample_rate;
us_lat = problem_data.transmitters.unknown.coord(1);
us_lon = problem_data.transmitters.unknown.coord(2);

% Reference transmitter
fRS = problem_data.transmitters.reference.freq * 1e6;
sr_ltess = problem_data.config.sample_rate_ltess;
rs_lat = problem_data.transmitters.reference.coord(1);
rs_lon = problem_data.transmitters.reference.coord(2);

% Sensors
NUM_SENSORS = length(problem_data.sensors);
NUM_FILES_SENSOR = problem_data.config.files_per_sensor;

folder_location = '/home/ygglc/desktop/phd/research/localization/core/localization-tdoa/data/real_tests';
localization_files = cell(NUM_SENSORS, NUM_FILES_SENSOR);
ltess_files = cell(NUM_SENSORS, NUM_FILES_SENSOR);

sensors = zeros(NUM_SENSORS, 3);
distances_rs = zeros(NUM_SENSORS,1);
distances_us = distances_rs;

X = zeros(NUM_FILES_SENSOR,3);
lls = zeros(NUM_FILES_SENSOR,2);
nlls = zeros(NUM_FILES_SENSOR,3);
nlop = zeros(NUM_FILES_SENSOR,3);
nlop_l1 = zeros(NUM_FILES_SENSOR,3);
nlop_l1l2 = zeros(NUM_FILES_SENSOR,3);

combinations = nchoosek(1:NUM_SENSORS,2);
% combinations = combinations(1:NUM_SENSORS-1,:);
N = size(combinations,1);
doa_meters = zeros(N,NUM_FILES_SENSOR);
doa_samples = zeros(N,NUM_FILES_SENSOR);

ell = referenceEllipsoid('WGS84');

correct_offset = true;
for EXPERIMENT = 1:NUM_FILES_SENSOR

    for ii = 1:NUM_SENSORS
        sensors(ii,1:2) = problem_data.sensors(ii).coordinates;
        sensors(ii,3) = problem_data.sensors(ii).height;
    %     distances_rs(ii) = 1000*deg2km(distance(rs_lat, rs_lon, sensors(ii,1), sensors(ii,2)));
    %     distances_us(ii) = 1000*deg2km(distance(us_lat, us_lon, sensors(ii,1), sensors(ii,2)));
    %     
        distances_rs(ii) = havdist([rs_lat, rs_lon], sensors(ii,1:2));
%         [distances_rs(ii),~] = distance(rs_lat,rs_lon,sensors(ii,1),sensors(ii,2),ell);
        distances_us(ii) = havdist([rs_lat, rs_lon], sensors(ii,1:2));

        localization_struct = dir([folder_location, filesep,  problem_data.config.folder, filesep, problem_data.sensors(ii).name, filesep, '*localization*']);
        ltess_struct = dir([folder_location, filesep, problem_data.config.folder, filesep, problem_data.sensors(ii).name, filesep, '*ltess*']);
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

        if correct_offset
            PPM = inf; count = 0; MAX = 10;
            while abs(PPM) > tol && count < MAX
                [PPM, PPM2] = ltess(signal_ltess{ii}, sr_ltess);
                if ~isnan(PPM)
                    signal_ltess{ii} = correct_fo_ltess(signal_ltess{ii},PPM,sr_ltess, 806e6);
                    signal_sensor = correct_fo(signal_sensor, PPM, sr_local, fRS, fUS);
                else
                    disp('[WARNING] Invalid file');
                end

                count = count + 1;
            end
        end
        signal_local{ii} = signal_sensor;
    end

    %% TDOA estimation
    % Signal parameters
    max_lag = round(100*1000*sr_local/3e8);
    corr_type = 'dphase';
    report_level = 1;
    signal_bandwidth_khz = 0;
    interpol_factor = 5;
    ns_freq = problem_data.config.samples_per_slice;
    ns_slice = 0.7 * ns_freq;

    us_doa = zeros(N,1);
    keep = zeros(N,1);
    for ii = 1:N
        % Select the sensors and compute distance differences
        si = combinations(ii,1);
        sj = combinations(ii,2);
        rs_diff_ij = distances_rs(si) - distances_rs(sj);
        us_doa(ii) = distances_us(si) - distances_us(sj);

        % Compute TDOA and store it
        [doa_meters(ii,EXPERIMENT), doa_samples(ii,EXPERIMENT), ~, ~, ~, ~, keep(ii)] = ...
            tdoa2(signal_local{si}, signal_local{sj}, ns_freq, ns_slice, sr_local, rs_diff_ij, ...
                  max_lag, corr_type, report_level, signal_bandwidth_khz, interpol_factor);
    end

    %% Multilateration
    % Data conversion
    sensors_ecef = llh2ecef(sensors);
    center_ecef = mean(sensors_ecef,1);
    sensors_ecef = sensors_ecef - center_ecef;

    % Optimization routines
    % With latitude/longitude
    X0 = mean(sensors,1);
    [xapprox, yapprox] = latlong2xy(sensors(:,1),sensors(:,2),X0(1),X0(2));
    % Fixing higher delay than distance
    xapprox = xapprox/0.995;
    yapprox = yapprox/0.995;

    % With ECEF coordinates
    [x0_init,y0_init] = solution2d(doa_meters(1:NUM_SENSORS-1,EXPERIMENT),[xapprox,yapprox]);
    [lat,lon] = xy2latlong(x0_init(1),y0_init(1),X0(1),X0(2));
    
    lls(EXPERIMENT,:) = [lat,lon]

    X0_init = llh2ecef([lat,lon,0])-center_ecef;
    % X0 = mean(sensors_ecef,1);
    % X0 = llh2ecef([rs_lat, rs_lon, 600]) - center_ecef;
    options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','Display','none','OptimalityTolerance',1e-10,'StepTolerance',1e-10);

    optimfun = @(X) optimnlls_ecef(X, doa_meters(:,EXPERIMENT), sensors_ecef, combinations);
    X_ecef = lsqnonlin(optimfun,X0_init,[],[],options);
    nlls(EXPERIMENT,:) = ecef2llh(X_ecef + center_ecef)

%     % Unconstrained minimization
%     options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'Display','none','OptimalityTolerance',1e-10,'StepTolerance',1e-14);
%     minimfun = @(X) optimmin_ecef(X, doa_meters(:,EXPERIMENT), sensors_ecef, combinations, 1e-6);
%     X_min_ecef = fminunc(minimfun,X0_init,options);
%     nlop(EXPERIMENT,:) = ecef2llh(X_min_ecef + center_ecef)
%     
%     % L1 norm minimization
%     options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'Display','iter','OptimalityTolerance',1e-14,'StepTolerance',1e-14);
%     minimfunl1 = @(X) optimminl1_ecef(X, doa_meters(:,EXPERIMENT), sensors_ecef, combinations);
%     X_min_l1_ecef = fminunc(minimfunl1,X0_init,options);
%     nlop_l1(EXPERIMENT,:) = ecef2llh(X_min_l1_ecef + center_ecef)
%     
%     % L1-L2 minimization
%     options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'Display','iter','OptimalityTolerance',1e-14,'StepTolerance',1e-14);
% 
%     % Initial error
%     d = ecef_distance(sensors_ecef, X0_init);
%     t = zeros(size(combinations,1),1);
%     si = combinations(:,1);
%     sj = combinations(:,2);
% 
%     t = d(si) - d(sj);
% 
%     err_prev = doa_meters - t;
% 
%     minimfunl1l2 = @(X) optimminl1l2_ecef(X, doa_meters(:,EXPERIMENT), sensors_ecef, combinations);
%     X_min_l1l2_ecef = fminunc(minimfunl1l2,X0_init,options);
%     nlop_l1l2(EXPERIMENT,:) = ecef2llh(X_min_l1l2_ecef + center_ecef)

end