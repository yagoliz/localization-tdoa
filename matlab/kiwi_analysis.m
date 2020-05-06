% =========================================================================
%  TDOA Localization with KiwiSDR data
%  Author: Yago Lizarribar
% =========================================================================

clear;
clc;
close all;

% adds subfolder with functions to PATH
[p,n,e] = fileparts(mfilename('fullpath'));
addpath([p '/tdoa']);
addpath([p '/ltess']);

%% Constants
realPos = [50.01556, 9.01083];
c = 299704644.54; %m/s

%% Load the config
load('kiwi/tdoa_data');
filesensor1 = files{1};
filesensor2 = files{2};
filesensor3 = files{3};

%% Position of the sensors
% RX and Ref TX Position
rx1_lat = input(1).coord(1); % RX 1
rx1_long = input(1).coord(2);

rx2_lat = input(2).coord(1); % RX 2
rx2_long = input(2).coord(2);

rx3_lat = input(3).coord(1); % RX 3
rx3_long = input(3).coord(2);

% Reference (midpoint of both transmitters)
geo_ref_lat  = mean([rx1_lat, rx2_lat, rx3_lat]);
geo_ref_long = mean([rx1_long, rx2_long, rx3_long]);

% Convert [latitude,longitude] to [x,y]
[rx1_x, rx1_y] = latlong2xy(rx1_lat, rx1_long, geo_ref_lat, geo_ref_long);
[rx2_x, rx2_y] = latlong2xy(rx2_lat, rx2_long, geo_ref_lat, geo_ref_long);
[rx3_x, rx3_y] = latlong2xy(rx3_lat, rx3_long, geo_ref_lat, geo_ref_long);

%% Read the IQ
% Sensor 1
[signal1Mat, fs1] = audioread(filesensor1);
iqsignal1 = complex(signal1Mat(:,1), signal1Mat(:,2));

% Sensor 2
[signal2Mat, fs2] = audioread(filesensor2);
iqsignal2 = complex(signal2Mat(:,1), signal2Mat(:,2));

% Sensor 3
[signal3Mat, fs3] = audioread(filesensor3);
iqsignal3 = complex(signal3Mat(:,1), signal3Mat(:,2));

if fs1 ~= fs2 || fs1 ~= fs3
    error('Sampling frequencies are not the same');
end

%% Process tdoa for sensor pairs
doa12 = tdoa_kiwi(iqsignal1, iqsignal2, 5, 'dphase');
doa13 = tdoa_kiwi(iqsignal1, iqsignal3, 5, 'dphase');
doa23 = tdoa_kiwi(iqsignal2, iqsignal3, 5, 'dphase');

% Conver to meters
doa12_meters = (doa12/fs1) * c;
doa13_meters = (doa13/fs1) * c;
doa23_meters = (doa23/fs1) * c;

doa_array = [doa12_meters, doa13_meters, doa23_meters];
sensors = [rx1_x, rx1_y; rx2_x, rx2_y; rx3_x, rx3_y];

NUM_SENSORS = 3;
combinations = nchoosek(1:NUM_SENSORS,2);
NUM_HYPERBOLAS = length(combinations);
hyp_array = cell(NUM_HYPERBOLAS,1);
% NUM_HYPERBOLAS = NUM_SENSORS - 1;
% hyp_array = cell(NUM_HYPERBOLAS,1);
% doa_array = zeros(NUM_HYPERBOLAS,1);
% combinations = ones(NUM_HYPERBOLAS,2);
% combinations(:,2) = combinations(:,2) + [1:NUM_SENSORS-1]';

for ii = 1:NUM_HYPERBOLAS
    sensor_1 = combinations(ii,1);
    sensor_2 = combinations(ii,2);
    hyp_array{ii} = hyperbola(doa_array(ii), sensors(sensor_1,:), sensors(sensor_2,:));
    doa_array(ii) = doa_array(ii);
end

%% Plot area
figure();

xmin = min([rx1_x, rx2_x, rx3_x]);
xmax = max([rx1_x, rx2_x, rx3_x]);
ymin = min([rx1_y, rx2_y, rx3_y]);
ymax = max([rx1_y, rx2_y, rx3_y]);

hold on; grid on;
xlim([xmin, xmax]); ylim([ymin, ymax]);
plot(realPos(1,1), realPos(1,2), 'kx', 'LineWidth', 3);
xlabel('X axis (m)');
ylabel('Y axis (m)');

% Plot sensors
for ii = 1:NUM_SENSORS
   plot(sensors(ii,1), sensors(ii,2), 'rx', 'LineWidth', 3); 
end

for ii = 1:length(combinations)
    hyp_points = hyp_array{ii};
    plot(hyp_points(1,:), hyp_points(2,:), 'b.-');
end


