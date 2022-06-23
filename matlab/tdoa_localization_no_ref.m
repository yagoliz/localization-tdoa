function tdoa_res = tdoa_localization(config)

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
if isstring(config) || ischar(config)
    problem_data = json_load(config);
elseif isstruct(config)
    problem_data = config;
else
    error("Unsupported data type for configuration. It must be filename or struct.")
end

% Sensors
NUM_SENSORS = length(problem_data.sensors);

EXP = problem_data.config.filenum;

folder_location = problem_data.config.folder;
if ~isfield(problem_data.config, "folder_date")
    folder_date = [];
else
    folder_date = problem_data.config.folder_date;
end
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

% Add corrections
ell = referenceEllipsoid('WGS84');
correct_offset = problem_data.config.correct;
if ~isfield(problem_data.config, "use_lte")
    use_lte = true;
else
    use_lte = problem_data.config.use_lte;
end

if ~isfield(problem_data.config, "correct_multipath")
    correct_multipath = false;
else
    correct_multipath = problem_data.config.correct_multipath;
end

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

% Reference transmitter
fRS = problem_data.transmitters.reference.freq * 1e6;
sr_ltess = problem_data.config.sample_rate_ltess;
rs_lat = problem_data.transmitters.reference.coord(1);
rs_lon = problem_data.transmitters.reference.coord(2);

% Unknown transmitter
fUS = problem_data.transmitters.unknown.freq * 1e6;
sr_local = problem_data.config.sample_rate;
if (isfield(problem_data.transmitters.unknown, 'coord'))
    us_lat = problem_data.transmitters.unknown.coord(1);
    us_lon = problem_data.transmitters.unknown.coord(2);
    distances_us = zeros(NUM_SENSORS,1);
    theoretical_doa = zeros(N,1);
end

% Opening loop
for ii = 1:NUM_SENSORS
    sensors(ii,1:2) = problem_data.sensors(ii).coordinates;
    sensors(ii,3) = problem_data.sensors(ii).height;

    distances_rs(ii) = distance(rs_lat, rs_lon, sensors(ii,1), sensors(ii,2), ell);

    if (isfield(problem_data.transmitters.unknown, 'coord'))
        distances_us(ii) = distance(us_lat, us_lon, sensors(ii,1), sensors(ii,2), ell);
    end

    localization_struct = dir(join([folder_location, filesep,  folder_date, filesep, problem_data.sensors(ii).name, filesep, 'E', num2str(EXP), '-', '*localization*'],''));
    ltess_struct = dir(join([folder_location, filesep, folder_date, filesep, problem_data.sensors(ii).name, filesep, 'E', num2str(EXP), '-', '*ltess*'],''));
    
    localization_files{ii} = join([localization_struct.folder, filesep, localization_struct.name],'');
    ltess_files{ii} = join([ltess_struct.folder, filesep, ltess_struct.name],'');
end

%% PPM correction
signal_local = cell(NUM_SENSORS,1);
signal_ltess = cell(NUM_SENSORS,1);
for ii = 1:NUM_SENSORS
    signal_sensor = spec_load(localization_files{ii});
    signal_ltess{ii} = spec_load(ltess_files{ii});

    if correct_offset && use_lte
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
if ~isfield(problem_data.config, "bw_us")
    bw_us = 0;
else
    bw_us = problem_data.config.bw_us;
end
interpol_factor = problem_data.config.interp;
ns_freq = problem_data.config.samples_per_slice;
ns_slice = 0.7 * ns_freq;

multipath_delays = zeros(NUM_SENSORS,1);
for ii = 1:N
    % Select the sensors and compute distance differences
    si = combinations(ii,1);
    sj = combinations(ii,2);
    
    % Compute theoretical TDOA in case US signal coordinates are provided
    if (isfield(problem_data.transmitters.unknown, 'coord'))
        theoretical_doa(ii) = distances_us(si) - distances_rs(sj);
    end

    % Compute TDOA and store it
    [doa_meters(ii), doa_samples(ii), doa_meters2(ii), doa_samples2(ii), ~, ~, multipath_delays(sj)] = ...
        tdoa2_no_ref(signal_local{si}, signal_local{sj}, ns_freq, ns_slice, sr_local, ...
              max_lag, corr_type, report_level, bw_us, interpol_factor, correct_multipath, multipath_delays(si));
end

%% Multilateration
% Data conversion
sensors_ecef = llh2ecef(sensors);
center_ecef = mean(sensors_ecef,1);
sensors_ecef_centered = sensors_ecef - center_ecef;

% With ECEF coordinates
if correct_offset && ~use_lte
    tdoa_values = doa_meters;
else
    tdoa_values = doa_meters2;
end

% Non Linear Optimization routines
X0_init = [0, 0, 0];
offset_init = zeros(1,NUM_SENSORS);
for ii = 2:NUM_SENSORS
    offset_init(ii) = -doa_meters(ii-1);
end
options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','Display','none','OptimalityTolerance',1e-18,'StepTolerance',1e-18);

optimfun = @(X) optimmin_ecef_no_offset(X, tdoa_values, sensors_ecef_centered, combinations);
X_ecef = lsqnonlin(optimfun,[X0_init,offset_init],[],[],options);
nlls(:) = ecef2llh(X_ecef(1:3)+center_ecef);

%% Output preparation
% Prepare output
tdoa_res = struct();
tdoa_res.res_linear = lls;
tdoa_res.res_accurate = nlls;

if (isfield(problem_data.transmitters.unknown, 'coord'))
    tdoa_res.error_linear = havdist(lls', [us_lat, us_lon]);
    tdoa_res.error_nonlin = havdist(nlls(1:2)', [us_lat, us_lon]);
end

fprintf("%f,%f\n",tdoa_res.res_linear(1),tdoa_res.res_linear(2));
fprintf("Error LLS: %f\n",havdist(tdoa_res.res_linear(1:2)', [us_lat, us_lon]));
fprintf("%f,%f\n",tdoa_res.res_accurate(1),tdoa_res.res_accurate(2));
fprintf("Error NLLS: %f\n",havdist(tdoa_res.res_accurate(1:2)', [us_lat, us_lon]))
end