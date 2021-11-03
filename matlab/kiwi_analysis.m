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

%% Data
load(['kiwi/data/dcf77_1_pre.mat']);

%% Constants
% Earth Radius
R = 6371;
c = 299792.458 * 1e3;

% Sensor coordinates
% sensors = [52.4862, 13.4761;...
%            46.1700, 14.3000;...
%            45.7793,  0.6146];
       
sensors = [input(1).coord; input(2).coord; input(3).coord];
       
dcf77 = [50.01556, 9.01083];
% dcf77 = [52.296667, -2.105278];


% Get distances from transmitter to sensor
sensorDistances = zeros(3,1);
for ii = 1:3
    sensorDistances(ii) = R * (pi/180) * distance(sensors(ii,1), sensors(ii,2), dcf77(1), dcf77(2));
end

tdoaTheory = [sensorDistances(1) - sensorDistances(2);...
              sensorDistances(1) - sensorDistances(3);...
              sensorDistances(2) - sensorDistances(3)];
          
tdoaTheory = tdoaTheory ./ c;
tdoaTheorySamples = tdoaTheory * 12e3;

%% Loop
numfiles = 30;
correlations = zeros(numfiles,3);
errors = zeros(numfiles,1);
errors1 = zeros(numfiles,1);
errors2 = zeros(numfiles,1);
errors3 = zeros(numfiles,1);

for ii = 1:numfiles
    % Load the config
%     load(['kiwi/data/dcf77_', num2str(ii), '_pre.mat']);
    load(['kiwi/data/tdf_',num2str(ii),'.mat']);
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
    [doa12, iqcorrelate12, corrfactor12, dt12] = tdoa_kiwi(iqsignal1, t1, iqsignal2, t2, mean([fs1,fs2]), 1, 'iq');
    [doa13, iqcorrelate13, corrfactor13, dt13] = tdoa_kiwi(iqsignal1, t1, iqsignal3, t3, mean([fs1,fs3]), 1, 'iq');
    [doa23, iqcorrelate23, corrfactor23, dt23] = tdoa_kiwi(iqsignal2, t2, iqsignal3, t3, mean([fs2,fs3]), 1, 'iq');
    
    correlations(ii,1) = doa12 / 12;
    correlations(ii,2) = doa13 / 12;
    correlations(ii,3) = doa23 / 12;

%     % Convert to meters
    doa12_meters = (doa12/mean([fs1,fs2])) * c;
    doa13_meters = (doa13/mean([fs1,fs3])) * c;
    doa23_meters = (doa23/mean([fs2,fs3])) * c;
    doa12_meters2 = (mean(dt12 + iqcorrelate12/mean([fs1, fs2]))) * c;
    doa13_meters2 = (mean(dt13 + iqcorrelate13/mean([fs1, fs3]))) * c;
    doa23_meters2 = (mean(dt23 + iqcorrelate23/mean([fs2, fs3]))) * c;

    % Multilateration
    doa_array = [doa12_meters, doa13_meters, doa23_meters];
    doa_array2 = [doa12_meters2, doa13_meters2, doa23_meters2];
    sensors = [rx1_x, rx1_y; rx2_x, rx2_y; rx3_x, rx3_y];
    
    % Let's do all the combinations
    % Sensor 1 as ref
    [x, y] = solution2d(doa_array(1:2), sensors);
    [x2, y2] = solution2d(doa_array2(1:2), sensors);
    [lat, long] = xy2latlong(x, y, geo_ref_lat, geo_ref_long);
    [lat2, long2] = xy2latlong(x2, y2, geo_ref_lat, geo_ref_long);
% 
%     % Error estimation
%     angle1 = distance(lat1(1), long1(1), dcf77(1), dcf77(2));
%     errors1(ii) = angle1 * (pi/180) * R;
%     
%     % Sensor2 as ref
%     doa2 = [-doa12_meters, doa23_meters];
%     sens2 = [sensors(2,:); sensors(1,:); sensors(3,:)];
%     [x2, y2] = solution2d(doa2, sens2);
%     [lat2, long2] = xy2latlong(x2, y2, geo_ref_lat, geo_ref_long);
%     angle2 = distance(lat2(1), long2(1), dcf77(1), dcf77(2));
%     errors2(ii) = angle2 * (pi/180) * R;
%     
%     % Sensor 3 as ref
%     doa3 = [-doa13_meters, -doa23_meters];
%     sens3 = [sensors(3,:); sensors(1,:); sensors(2,:)];
%     [x3, y3] = solution2d(doa3, sens3);
%     [lat3, long3] = xy2latlong(x3, y3, geo_ref_lat, geo_ref_long);
%     angle3 = distance(lat3(1), long3(1), dcf77(1), dcf77(2));
%     errors3(ii) = angle3 * (pi/180) * R;
%     
%     
%     xavg = mean([x(1), x2(1), x3(1)]);
%     yavg = mean([y(1), y2(1), y3(1)]);
%     [lat, long] = xy2latlong(xavg, yavg, geo_ref_lat, geo_ref_long);
    
    angle = distance(lat(1), long(1), dcf77(1), dcf77(2));
    angle2 = distance(lat2(1), long2(1), dcf77(1), dcf77(2));
    
    errors(ii) = angle * (pi/180) * R;
    errors2(ii) = angle2 * (pi/180) * R;
end

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


