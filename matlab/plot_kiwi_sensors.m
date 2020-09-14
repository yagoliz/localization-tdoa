clear;
clc;
close all;

% adds subfolder with functions to PATH
[p,n,e] = fileparts(mfilename('fullpath'));
addpath([p '/tdoa']);
addpath([p '/kiwi/m/']);

% Load the config
load(['kiwi/data/dcf77_', num2str(1), '_pre.mat']);

%% Position of the sensors
% RX and Ref TX Position
rx1_lat = input(1).coord(1); % RX 1
rx1_long = input(1).coord(2);

rx2_lat = input(2).coord(1); % RX 2
rx2_long = input(2).coord(2);

rx3_lat = input(3).coord(1); % RX 3
rx3_long = input(3).coord(2);

dcf77 = [50.01556, 9.01083];

% Plot the points
figure(); hold on; grid on;
geoplot([rx1_lat, rx2_lat, rx3_lat],[rx1_long, rx2_long, rx3_long],'rx',dcf77(1), dcf77(2), 'gx', 'MarkerSize',10,'LineWidth',5);
geolimits([min([rx1_lat, rx2_lat, rx3_lat]) - 5, max([rx1_lat, rx2_lat, rx3_lat]) + 5], [min([rx1_long, rx2_long, rx3_long]) - 5, max([rx1_long, rx2_long, rx3_long]) + 5]);
