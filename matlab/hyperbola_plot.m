% =========================================================================
%  Plotting hyperbolae
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
% Earth Radius
R = 6371 * 1e3;
c = 3e8;

% Sensor coordinates
sensors = [52.4862, 13.4761;...
           46.1700, 14.3000;...
           45.7793,  0.6146];
       
dcf77 = [50.01556, 9.01083];


%% Loop
correlations = zeros(30,3);
errors = zeros(30,1);
for ii = 2
    % Load the config
    load(['kiwi/data/dcf77_', num2str(ii), '_pre.mat']);
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
    
    [dcf_x, dcf_y] = latlong2xy(dcf77(1), dcf77(2), geo_ref_lat, geo_ref_long);

    %% Read the IQ
    % Sensor 1
    iqsignal1 = input(1).z;
    t1 = input(1).t;
    fs1 = input(1).fs;

    % Sensor 2
    iqsignal2 = input(2).z;
    t2 = input(2).t;
    fs2 = input(2).fs;

    % Sensor 3
    iqsignal3 = input(3).z;
    t3 = input(3).t;
    fs3 = input(3).fs;
    
    minLength = min([length(iqsignal1), length(iqsignal2), length(iqsignal3)]);
    iqsignal1 = double(iqsignal1(1:minLength));
    iqsignal2 = double(iqsignal2(1:minLength));
    iqsignal3 = double(iqsignal3(1:minLength));

    if abs(fs1-fs2) > 1 || abs(fs1-fs3) > 1
        error('Sampling frequencies are not the same');
    end

    %% Process tdoa for sensor pairs
%     iq1,t1,iq2,t2,fs,interpol_factor,corr_type
    [doa12, iqcorrelate12, corrfactor12] = tdoa_kiwi(iqsignal1, t1, iqsignal2, t2, fs1, 1, 'iq');
    [doa13, iqcorrelate13, corrfactor13] = tdoa_kiwi(iqsignal1, t1, iqsignal3, t3, fs1, 1, 'iq');
    [doa23, iqcorrelate23, corrfactor23] = tdoa_kiwi(iqsignal2, t2, iqsignal3, t3, fs1, 1, 'iq');
    
    correlations(ii,1) = doa12 / 12;
    correlations(ii,2) = doa13 / 12;
    correlations(ii,3) = doa23 / 12;

%     % Conver to meters
    doa12_meters = (doa12/fs1) * c;
    doa13_meters = (doa13/fs1) * c;
    doa23_meters = (doa23/fs1) * c;

    % Multilateration
    doa_array = [doa12_meters, doa13_meters, doa23_meters];
    sensors = [rx1_x, rx1_y; rx2_x, rx2_y; rx3_x, rx3_y];
    
    % Let's do all the combinations
    % Sensor 1 as ref
    [x, y] = solution2d(doa_array(1:2), sensors);
%     [lat, long] = xy2latlong(x, y, geo_ref_lat, geo_ref_long);

    % Error estimation
%     angle = distance(lat(1), long(1), dcf77(1), dcf77(2));
    
    % Sensor2 as ref
    doa2 = [-doa12_meters, doa23_meters];
    sens2 = [sensors(2,:); sensors(1,:); sensors(3,:)];
    [x2, y2] = solution2d(doa2, sens2);
%     [lat2, long2] = xy2latlong(x2, y2, geo_ref_lat, geo_ref_long);
    
    % Sensor 3 as ref
    doa3 = [-doa13_meters, -doa23_meters];
    sens3 = [sensors(3,:); sensors(1,:); sensors(2,:)];
    [x3, y3] = solution2d(doa3, sens3);
%     [lat3, long3] = xy2latlong(x3, y3, geo_ref_lat, geo_ref_long);
    
    xavg = mean([x(1), x2(1), x3(1)]);
    yavg = mean([y(1), y2(1), y3(1)]);
    [lat, long] = xy2latlong(xavg, yavg, geo_ref_lat, geo_ref_long);
    
    combinations = nchoosek(1:3,2);
    NUM_HYPERBOLAS = length(combinations);
    hyp_array = cell(NUM_HYPERBOLAS,1);
%     doa_array = zeros(NUM_HYPERBOLAS,1);

    for jj = 1:NUM_HYPERBOLAS
        sensor_1 = combinations(jj,1);
        sensor_2 = combinations(jj,2);
        hyp_array{jj} = hyperbola(doa_array(jj), sensors(sensor_1,:), sensors(sensor_2,:));
    end
    
    %% Plot area
    figure();

    hold on; grid on;
%     xlim([xmin, xmax]); ylim([ymin, ymax]);
%     plot(dcf_x, dcf_y, 'kx', 'LineWidth', 3);
    xlabel('X axis (m)');
    ylabel('Y axis (m)');

    % Plot sensors
    for jj = 1:3
       plot(sensors(jj,1), sensors(jj,2), 'rx', 'LineWidth', 3); 
    end

    for jj = 1:length(combinations)
        hyp_points = hyp_array{jj};
        plot(hyp_points(1,:), hyp_points(2,:), 'b.-');
    end
    
    plot(x , y , 'gx', 'LineWidth', 3);
    plot(x2, y2, 'gx', 'LineWidth', 3);
    plot(x3, y3, 'gx', 'LineWidth', 3);
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