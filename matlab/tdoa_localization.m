% =========================================================================
%  Real test realization
%  Author: Yago Lizarribar
% =========================================================================
function tdoa_res = tdoa_localization(config_file) %#codegen

warning ('off','all');

% Speed of light
global c; 
c = 299792458; % m/s

% adds subfolder with functions to PATH
[p,~,~] = fileparts(mfilename('fullpath'));
addpath([p '/tdoa']);
addpath([p '/ltess']);
addpath([p, '/geodesy']);
addpath([p, '/optimization']);

%% Load configuration
problem_data = json_load(config_file);

% Unknown transmitter
fUS = problem_data.transmitters.unknown.freq * 1e6;
sr_local = problem_data.config.sample_rate;

% Reference transmitter
fRS = problem_data.transmitters.reference.freq * 1e6;
sr_ltess = problem_data.config.sample_rate_ltess;
rs_lat = problem_data.transmitters.reference.coord(1);
rs_lon = problem_data.transmitters.reference.coord(2);

% Sensors
NUM_SENSORS = length(problem_data.sensors);

EXP = problem_data.config.filenum;

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
doa_meters2 = zeros(N,1);
doa_samples = zeros(N,1);
doa_samples2 = zeros(N,1);

% Add correction
ell = referenceEllipsoid('WGS84');
correct_offset = problem_data.config.correct;

% Plot heatmap
if ~isfield(problem_data.config, "plot_heatmap")
    plot_heatmap = false;
else
    plot_heatmap = problem_data.config.plot_heatmap;
end

% Specify hyperbola
if ~isfield(problem_data.config, "generate_hyperbola")
    generate_hyperbola = false;
else
    generate_hyperbola = problem_data.config.generate_hyperbola;
end

% Opening loop
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
max_lag = round(50*1000*sr_local/c);
if ~isfield(problem_data.config,"corr_method")
    corr_type = "dphase";
else
    corr_type = problem_data.config.corr_method;
end

report_level = 1;
if ~isfield(problem_data.config, "bw_rs")
    bw_rs = 0;
else
    bw_rs = problem_data.config.bw_rs;
end
if ~isfield(problem_data.config, "bw_us")
    bw_us = 0;
else
    bw_us = problem_data.config.bw_us;
end
interpol_factor = problem_data.config.interp;
ns_freq = problem_data.config.samples_per_slice;
ns_slice = 0.7 * ns_freq;

keep = zeros(N,1);
for ii = 1:N
    % Select the sensors and compute distance differences
    si = combinations(ii,1);
    sj = combinations(ii,2);
    rs_diff_ij = distances_rs(si) - distances_rs(sj);

    % Compute TDOA and store it
    [doa_meters(ii), doa_samples(ii), doa_meters2(ii), doa_samples2(ii), ~, ~, keep(ii)] = ...
        tdoa2(signal_local{si}, signal_local{sj}, ns_freq, ns_slice, sr_local, rs_diff_ij, ...
              max_lag, corr_type, report_level, bw_rs, bw_us, interpol_factor);
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

% With ECEF coordinates
[x0_init,y0_init] = solution2d(doa_meters(1:NUM_SENSORS-1),[xapprox,yapprox]); % Fixing higher delay than distance
[lat,lon] = xy2latlong(x0_init(1),y0_init(1),X0(1),X0(2));

lls(:) = [lat,lon];

X0_init = llh2ecef([lat,lon,X0(3)])-center_ecef;
options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','Display','none','OptimalityTolerance',1e-12,'StepTolerance',1e-12);

optimfun = @(X) optimnlls_ecef(X, doa_meters, sensors_ecef, combinations);
X_ecef = lsqnonlin(optimfun,X0_init,[],[],options);
nlls(:) = ecef2llh(X_ecef+center_ecef);

%% Pretty results
% Calculate heatmap
lat_min = lls(1)-0.1; lat_max = lls(1)+0.1;
lon_min = lls(2)-0.1; lon_max = lls(2)+0.1;
[xrange,yrange] = latlong2xy([lat_min;lat_max],[lon_min;lon_max],X0(1),X0(2));
[p_x, p_y, mse_doa ] = heatmap(doa_meters(1:NUM_SENSORS-1), [xapprox, yapprox], xrange, yrange, combinations(1:NUM_SENSORS-1,:), 200);
[lats, lons] = xy2latlong(p_x, p_y, X0(1), X0(2));

if plot_heatmap
    figure();
    lati = unique(lats); loni = unique(lons);
    [LT, LN] = meshgrid(lati,loni);
    Z = reshape(mse_doa,size(LT));
    surf(LT,LN,Z);
    xlabel('Latitude'); ylabel('Longitude');
end

idx = mse_doa > 0.1;
lats = lats(idx);
lons = lons(idx);
mse_doa = mse_doa(idx);
hm = [lats,lons,mse_doa];

% Generate hyperbola
if generate_hyperbola
    t = 0:0.001:2;
    hyperbolas = zeros(2*(NUM_SENSORS-1),length(t)*4);
    for ii = 1:NUM_SENSORS-1
        perturbation = 0.5 / sr_local / interpol_factor * c;
        si = combinations(ii,1);
        sj = combinations(ii,2);
        h1 = hyperbola2d(doa_meters(ii)+perturbation,[xapprox(si),yapprox(si)],[xapprox(sj),yapprox(sj)],t);
        h2 = hyperbola2d(doa_meters(ii)-perturbation,[xapprox(si),yapprox(si)],[xapprox(sj),yapprox(sj)],t);
        h = [h1,fliplr(h2)];
        [hlats, hlons] = xy2latlong(h(1,:),h(2,:),X0(1),X0(2));
        hyperbolas(2*ii-1:2*ii,:) = [hlats; hlons];
    end
end

%% Output preparation
% Prepare output
tdoa_res = struct();
tdoa_res.hm = hm;
tdoa_res.res_linear = lls;
tdoa_res.res_accurate = nlls;
if generate_hyperbola
    tdoa_res.hyperbolas = hyperbolas;
else
    tdoa_res.hyperbolas = [];
end

fprintf("%f,%f\n",tdoa_res.res_linear(1),tdoa_res.res_linear(2));
fprintf("%f,%f\n",tdoa_res.res_accurate(1),tdoa_res.res_accurate(2));
end