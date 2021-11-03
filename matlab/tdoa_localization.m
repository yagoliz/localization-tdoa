% =========================================================================
%  Real test realization
%  Author: Yago Lizarribar
% =========================================================================
function tdoa_res = tdoa_localization(config_file)

warning ('off','all');

% Speed of light
global c err_prev; 
c = 299792458; % m/s

% adds subfolder with functions to PATH
[p,n,e] = fileparts(mfilename('fullpath'));
addpath([p '/tdoa']);
addpath([p '/ltess']);
addpath([p, '/geodesy']);
addpath([p, '/optimization']);

%% Load configuration
problem_data = json_load(config_file);

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

EXP = problem_data.config.files_per_sensor;

folder_location = problem_data.config.folder;
localization_files = cell(NUM_SENSORS,1);
ltess_files = cell(NUM_SENSORS,1);

sensors = zeros(NUM_SENSORS, 3);
distances_rs = zeros(NUM_SENSORS,1);

lls = zeros(2,1);
nlls = zeros(3,1);

combinations = nchoosek(1:NUM_SENSORS,2);
N = size(combinations,1);
doa_meters = zeros(N,1);
doa_samples = zeros(N,1);

ell = referenceEllipsoid('WGS84');
correct_offset = problem_data.config.correct;
for ii = 1:NUM_SENSORS
    sensors(ii,1:2) = problem_data.sensors(ii).coordinates;
    sensors(ii,3) = problem_data.sensors(ii).height;

%     distances_rs(ii) = havdist([rs_lat, rs_lon], sensors(ii,1:2));
%     distances_us(ii) = havdist([rs_lat, rs_lon], sensors(ii,1:2));
    distances_rs(ii) = distance(rs_lat, rs_lon, sensors(ii,1), sensors(ii,2), ell);

    localization_struct = dir([folder_location, filesep,  problem_data.config.folder_date, filesep, problem_data.sensors(ii).name, filesep, 'E', num2str(EXP), '-', '*localization*']);
    ltess_struct = dir([folder_location, filesep, problem_data.config.folder_date, filesep, problem_data.sensors(ii).name, filesep, 'E', num2str(EXP), '-', '*ltess*']);
    
    localization_files{ii} = [localization_struct.folder, filesep, localization_struct.name];
    ltess_files{ii} = [ltess_struct.folder, filesep, ltess_struct.name];
end

%% PPM correction
signal_local = cell(NUM_SENSORS,1);
signal_ltess = cell(NUM_SENSORS,1);
for ii = 1:NUM_SENSORS
    signal_sensor = spec_load(localization_files{ii});
    signal_ltess{ii} = spec_load(ltess_files{ii});

    if correct_offset

        [PPM, PPM2] = ltess(signal_ltess{ii}, sr_ltess);
        if ~isnan(PPM)
            signal_sensor = correct_fo(signal_sensor, PPM, sr_local, fRS, fUS);
        else
            disp('[WARNING] Invalid file');
        end
    end
    signal_local{ii} = signal_sensor;
end

%% TDOA estimation
% Signal parameters
max_lag = round(100*1000*sr_local/3e8);
corr_type = problem_data.config.corr_method;
report_level = 1;
signal_bandwidth_khz = 0;
interpol_factor = problem_data.config.interp;
ns_freq = problem_data.config.samples_per_slice;
ns_slice = 0.5 * ns_freq;

keep = zeros(N,1);
for ii = 1:N
    % Select the sensors and compute distance differences
    si = combinations(ii,1);
    sj = combinations(ii,2);
    rs_diff_ij = distances_rs(si) - distances_rs(sj);

    % Compute TDOA and store it
    [doa_meters(ii), doa_samples(ii), ~, ~, ~, ~, keep(ii)] = ...
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
[x0_init,y0_init] = solution2d(doa_meters(1:NUM_SENSORS-1),[xapprox,yapprox]);
[lat,lon] = xy2latlong(x0_init(1),y0_init(1),X0(1),X0(2));

lls(:) = [lat,lon];

X0_init = llh2ecef([lat,lon,600])-center_ecef;
% X0 = mean(sensors_ecef,1);
% X0 = llh2ecef([rs_lat, rs_lon, 600]) - center_ecef;
options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','Display','none','OptimalityTolerance',1e-10,'StepTolerance',1e-10);

optimfun = @(X) optimnlls_ecef(X, doa_meters, sensors_ecef, combinations);
X_ecef = lsqnonlin(optimfun,X0_init,[],[],options);
nlls(:) = ecef2llh(X_ecef + center_ecef);

% Calculate heatmap
lat_min = nlls(1)-0.02; lat_max = nlls(1)+0.02;
lon_min = nlls(2)-0.02; lon_max = nlls(2)+0.02;
[xrange,yrange] = latlong2xy([lat_min;lat_max],[lon_min;lon_max],X0(1),X0(2));
[p_x, p_y, mse_doa ] = heatmap(doa_meters, [xapprox, yapprox], xrange, yrange, combinations, 200);
[lats, lons] = xy2latlong(p_x, p_y, X0(1), X0(2));

idx = mse_doa > 0.1;
lats = lats(idx);
lons = lons(idx);
mse_doa = mse_doa(idx);
hm = [lats,lons,mse_doa];

% Prepare output
tdoa_res = struct();
tdoa_res.hm = hm;
tdoa_res.res_linear = lls;
tdoa_res.res_accurate = nlls;
end