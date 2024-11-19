function times = xcorr_bench(config)

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

E = problem_data.config.exp;
times = zeros(1,E);

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

report_level = 0;
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
ns_slice = 0.1 * ns_freq;

for ii = 1:E
    % Select the sensors and compute distance differences
    si = combinations(1,1);
    sj = combinations(1,2);
    rs_diff_ij = distances_rs(si) - distances_rs(sj);

    % Compute TDOA and store it
    tic
    tdoa2(signal_local{si}, signal_local{sj}, ns_freq, ns_slice, sr_local, rs_diff_ij, ...
              max_lag, corr_type, report_level, bw_rs, bw_us, interpol_factor,...
              problem_data.sensors(si).name, problem_data.sensors(sj).name);
    times(ii) = toc;
end
end