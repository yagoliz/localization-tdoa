% =========================================================================
%  TDOA Simulator
%  Author: Yago Lizarribar
% =========================================================================

clear;
clc;
close all;

% Let's add the path to the functions
addpath('tdoa');

%% Geometrical setup
% Area setup
xmin = 0; xmax = 1000;
ymin = 0; ymax = 1000;

% Positions of the sensors
NUM_SENSORS = 5;
sensor = zeros(NUM_SENSORS, 2);

% Let's set the positions manually
sensor(1,:) = [100, 200];
sensor(2,:) = [600, 750];
sensor(3,:) = [800, 500];
sensor(4,:) = [200, 800];
sensor(5,:) = [600, 100];

%% Fake transmitter simulation
transmitter = [500, 500];
sensordist = sqrt(sum((sensor - transmitter).^2,2));

%% Calculate the hyperbolas
% combinations = nchoosek(1:NUM_SENSORS,2);
% NUM_HYPERBOLAS = length(combinations);
% hyp_array = cell(NUM_HYPERBOLAS,1);
% doa_array = zeros(NUM_HYPERBOLAS,1);
NUM_HYPERBOLAS = NUM_SENSORS - 1;
hyp_array = cell(NUM_HYPERBOLAS,1);
doa_array = zeros(NUM_HYPERBOLAS,1);
combinations = ones(NUM_HYPERBOLAS,2);
combinations(:,2) = combinations(:,2) + [1:NUM_SENSORS-1]';

for ii = 1:NUM_HYPERBOLAS
    sensor_1 = combinations(ii,1);
    sensor_2 = combinations(ii,2);
    doa_s1_s2 = sensordist(sensor_1) - sensordist(sensor_2) + 50*rand(1)j;
    hyp_array{ii} = hyperbola(doa_s1_s2, sensor(sensor_1,:), sensor(sensor_2,:));
    doa_array(ii) = doa_s1_s2;
end

%% Get the heatmap
[heat_x, heat_y, mse_doa] = heatmap(doa_array, sensor, [xmin, xmax], [ymin, ymax], combinations, 1000);

%% Get the optimal point
[x, y] = solution2d(doa_array, sensor);

%% Plot area
figure();

hold on; grid on;
xlim([xmin, xmax]); ylim([ymin, ymax]);
plot(transmitter(1,1), transmitter(1,2), 'kx', 'LineWidth', 3);

for ii = 1:length(x)
    plot(x(ii), y(ii), 'gx', 'LineWidth', 3);
end

% Plot sensors
for ii = 1:NUM_SENSORS
   plot(sensor(ii,1), sensor(ii,2), 'rx', 'LineWidth', 3); 
end

for ii = 1:length(combinations)
    hyp_points = hyp_array{ii};
    plot(hyp_points(1,:), hyp_points(2,:), 'b.-');
end

%% Let's plot the heatmap
figure();

% Plot the contours/image
% imagesc(heat_x, heat_y, flipud(rot90(log10(mse_doa)))); colorbar; colormap jet;
% set(gca,'YDir','normal')
h = surf(heat_x, heat_y, rot90(log10(mse_doa))); view(2); colorbar; colormap jet;
set(h, 'edgecolor', 'none');
xlabel('X axis');
ylabel('Y axis');

hold on; grid on;
xlim([xmin, xmax]); ylim([ymin, ymax]);
plot(transmitter(1,1), transmitter(1,2), 'kx', 'LineWidth', 3);

for ii = 1:length(x)
    plot(x(ii), y(ii), 'gx', 'LineWidth', 3);
end

% Plot sensors
for ii = 1:NUM_SENSORS
   plot(sensor(ii,1), sensor(ii,2), 'rx', 'LineWidth', 3); 
end
