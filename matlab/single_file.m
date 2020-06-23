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
addpath([p '/kiwi/m/']);

%% Constants
realPos = [50.01556, 9.01083];
c = 299704644.54; % m/s
R = 6371; % km

%% Load the config
load(['kiwi/data/dcf77_bw.mat']);
input = alignInputs(input);

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
iqsignal1 = input(1).z;
fs1 = input(1).fs;

% Sensor 2
iqsignal2 = input(2).z;
fs2 = input(2).fs;

% Sensor 3
iqsignal3 = input(3).z;
fs3 = input(3).fs;

if abs(fs1-fs2) > 1 || abs(fs1-fs3) > 1
    error('Sampling frequencies are not the same');
end

%% Process tdoa for sensor pairs
doa12 = tdoa_kiwi(iqsignal1, iqsignal2, fs1, 1, 'dphase');
doa13 = tdoa_kiwi(iqsignal1, iqsignal3, fs1, 1, 'dphase');
doa23 = tdoa_kiwi(iqsignal2, iqsignal3, fs1, 1, 'dphase');

% Conver to meters
doa12_meters = (doa12/fs1) * c;
doa13_meters = (doa13/fs1) * c;
doa23_meters = (doa23/fs1) * c;

% Multilateration
doa_array = [doa12_meters, doa13_meters, doa23_meters];
sensors = [rx1_x, rx1_y; rx2_x, rx2_y; rx3_x, rx3_y];

[x, y] = solution2d(doa_array(1:2), sensors);
[lat, long] = xy2latlong(x, y, geo_ref_lat, geo_ref_long);

% Error estimation
angle = distance(lat(1), long(1), realPos(1), realPos(2));
errors = angle * (pi/180) * R;

%% Alignment function
function inputAli = alignInputs(input)

  t0 = max([input(1).t(1), input(2).t(1), input(3).t(1)]);
  t1 = min([input(1).t(end), input(2).t(end), input(3).t(end)]);
  
  n = numel(input);
  for i=1:n
    if ~input(i).use
      continue
    end
    b = input(i).t<t0 | input(i).t>t1;
    input(i).t(b) = [];
    input(i).z(b) = [];
    input(i).use  = numel(input(i).z)/input(i).fs > 10;
    if ~input(i).use
      printf('tdoa_read_data: %-40s excluded (%.2f sec < %g sec overlap)\n', ...
             input(i).fn, numel(input(i).z)/input(i).fs, 10);
    end
  end
  
  inputAli = input;
end
  