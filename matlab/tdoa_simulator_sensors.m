% =========================================================================
%  TDOA Simulator
%  Author: Yago Lizarribar
% =========================================================================

clear;
clc;
close all;

% Let's add the path to the functions
addpath('tdoa');

%% Constants
c = 299704644.54; % Speed of light in the air (m/s)

%% Geometrical setup
% Area setup
xmin = 0; xmax = 10000;
ymin = 0; ymax = 10000;
r = 10000;

% Positions of the sensors
load('sensors');

%% Fake transmitter simulation
transmitter = [xmax/2, ymax/2];
sensordist = sqrt(sum((sensor - transmitter).^2,2));

%% Simulation
ds = 1;
sampling_rate = 2*1e6;
dt = ds / sampling_rate;

error_array = zeros(length(sensor)-2,1);

for NUM_SENSORS = 3:length(sensor)

    NUM_HYPERBOLAS = NUM_SENSORS - 1;
    hyp_array = cell(NUM_HYPERBOLAS,1);
    doa_array = zeros(NUM_HYPERBOLAS,1);
    combinations = ones(NUM_HYPERBOLAS,2);
    combinations(:,2) = combinations(:,2) + [1:NUM_SENSORS-1]';

    % Calculating the hyperbolas
    for ii = 1:NUM_HYPERBOLAS
        sensor_1 = combinations(ii,1);
        sensor_2 = combinations(ii,2);
        delta_meters = c * dt;
        doa_s1_s2 = sensordist(sensor_1) - sensordist(sensor_2) + delta_meters;
        hyp_array{ii} = hyperbola(doa_s1_s2, sensor(sensor_1,:), sensor(sensor_2,:));
        doa_array(ii) = doa_s1_s2;
    end

    % Get the heatmap
    [heat_x, heat_y, mse_doa] = heatmap(doa_array, sensor(1:NUM_SENSORS,:), [xmin, xmax], [ymin, ymax], combinations, r);

    % Now we obtain the maximum value
    [x, y] = find(mse_doa == max(mse_doa(:)));
    error_array(NUM_SENSORS-2) = norm(transmitter - [x*r/xmax,y*r/ymax]);


end

save('tdoa_sensors_2s.mat','error_array');