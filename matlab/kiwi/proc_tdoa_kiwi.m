% Translation from KiwiSDR TDOA processor
% Source: https://github.com/hcab14/TDoA

function [tdoa, input] = proc_tdoa_kiwi(directory, files, config)
%% Configuration setup
% Paths for all files
addpath([pwd, filesep(), 'm']);
addpath([pwd, filesep(), 'gnss_pos']);
addpath([pwd, filesep(), 'iq']);

for i=1:numel(files)
  input(i).fn = files{i};
end

config.dir       = directory;
config.plot_kiwi = true;
if isfield(config, 'lat_range')
  config.plot_kiwi_json = true;
  config = tdoa_autoresolution(config);
end
if ~isfield(config, 'use_constraints')
  config.use_constraints = false;
end

%% Compute lags
[input,status.input] = tdoa_read_data(config, input, directory);
[tdoa, status.cross_correlations] = tdoa_compute_lags_new(input);

if config.use_constraints
  [tdoa,status.cross_correlations] = tdoa_cluster_lags(config, tdoa, input, status.cross_correlations);
  [tdoa,input,status.constraints]  = tdoa_verify_lags (config, tdoa, input);
end

config.plotname = 'TDoA map';
config.title    = sprintf('%g kHz %s', input(1).freq, input(1).time);

[tdoa,status.position] = tdoa_plot_map(input, tdoa, config);
if config.new
  tdoa = tdoa_plot_dt_new(input, tdoa, config, 1e-2);
else
  tdoa = tdoa_plot_dt (input, tdoa, config, 2.5e-3);
end

%% Save into a .mat file (except the raw IQ samples and times)
for i=1:numel(input)
  input(i).t=[];
  input(i).z=[];
end
save('-mat', sprintf('%s/tdoa_data.mat', config.dir))

end
