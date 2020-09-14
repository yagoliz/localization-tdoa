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

% Load the config
load(['kiwi/data/dcf77_', num2str(1), '_pre.mat']);
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
%     geo_ref_lat = rx3_lat;
%     geo_ref_long = rx3_long;

% Convert [latitude,longitude] to [x,y]
[rx1_x, rx1_y] = latlong2xy(rx1_lat, rx1_long, geo_ref_lat, geo_ref_long);
[rx2_x, rx2_y] = latlong2xy(rx2_lat, rx2_long, geo_ref_lat, geo_ref_long);
[rx3_x, rx3_y] = latlong2xy(rx3_lat, rx3_long, geo_ref_lat, geo_ref_long);

%% Read the IQ
% Sensor 1
iqsignal1 = input(1).z;
fs1 = input(1).fs;
t1 = input(1).t;

% Sensor 2
iqsignal2 = input(2).z;
fs2 = input(2).fs;
t2 = input(2).t;

% Sensor 3
iqsignal3 = input(3).z;
fs3 = input(3).fs;
t3 = input(3).t;

minLength = min([length(iqsignal1), length(iqsignal2), length(iqsignal3)]);
iqsignal1 = double(iqsignal1(1:minLength));
iqsignal2 = double(iqsignal2(1:minLength));
iqsignal3 = double(iqsignal3(1:minLength));

if abs(fs1-fs2) > 1 || abs(fs1-fs3) > 1
    error('Sampling frequencies are not the same');
end

%% Process tdoa for sensor pairs
[doa13, iqcorrelate13, corrfactor13, dt13] = tdoa_kiwi(iqsignal1, t1, iqsignal3, t3, mean([fs1,fs3]), 1, 'iq');
[doa13, iqcorrelate13_5, corrfactor13, dt13] = tdoa_kiwi(iqsignal1, t1, iqsignal3, t3, mean([fs1,fs3]), 5, 'iq');
[doa13, iqcorrelate13_10, corrfactor13, dt13] = tdoa_kiwi(iqsignal1, t1, iqsignal3, t3, mean([fs1,fs3]), 10, 'iq');

figure(); grid on; hold on; plot(iqcorrelate13,'.--','MarkerSize',10); plot(iqcorrelate13_5,'.--','MarkerSize',10); plot(iqcorrelate13_10,'.--','MarkerSize',10);
xlabel('Chunk number'); ylabel('TDOA (samples)'); xlim([1,20]); legend('No upsampling', 'Upsampling x5', 'Upsampling x10','Location','SouthEast');


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