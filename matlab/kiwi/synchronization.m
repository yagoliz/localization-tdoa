%% Here we will analyze the correlation offsets for all the obtained data
clear; clc; close all;

% Variables
avgCorrAli = zeros(3,30);
offsets = zeros(3,30);

% Earth Radius
R = 6371 * 1e3;
c = 3e8;

% Sensor coordinates
sensors = [52.4862, 13.4761;...
           46.1700, 14.3000;...
           45.7793,  0.6146];
       
dcf77 = [50.01556, 9.01083];

% Get distances from transmitter to sensor
sensorDistances = zeros(3,1);
for ii = 1:3
    sensorDistances(ii) = R * (pi/180) * distance(sensors(ii,1), sensors(ii,2), dcf77(1), dcf77(2));
end

tdoaTheory = [sensorDistances(1) - sensorDistances(2);...
              sensorDistances(1) - sensorDistances(3);...
              sensorDistances(2) - sensorDistances(3)];
          
tdoaTheory = tdoaTheory ./ c;


%% Main loop
for ii = 1:30
    % We first do everything for unaligned samples
    filepathPre = ['data/dcf77_', num2str(ii), '_pre.mat'];
    load(filepathPre);
    
    fs = input(1).fs;
    
    % Now it's the turn for aligned samples
    inputAli = alignInputs(input);
    
    iq1Ali = inputAli(1).z;
    iq2Ali = inputAli(2).z;
    iq3Ali = inputAli(3).z;
    
    avgCorrAli(1,ii) = correlate(iq1Ali, iq2Ali, fs);
    avgCorrAli(2,ii) = correlate(iq1Ali, iq3Ali, fs);
    avgCorrAli(3,ii) = correlate(iq2Ali, iq3Ali, fs);
        
end

% Plotting
avgCorrAli = avgCorrAli / fs;
offsets = avgCorrAli - tdoaTheory;
offsets = abs(offsets) * 1e6;

figure(); hold on; grid on;
title('Synchronization errors between sensors');
boxplot(offsets', 'labels', {'Sensors 1 & 2', 'Sensors 1 & 3', 'Sensors 2 & 3'});
ylabel('Delay (\mus)');

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

%% Correlation function
function meanCorr = correlate(iq1, iq2, fs)
    jump = 10000;
    maxValue = 200000;
    iqcorrelate12 = zeros(maxValue/jump,1);
    
    data = 1;
    for i = 1:jump:maxValue-jump+1
        [r, lags] = xcorr(iq1(i:i + (jump-1)), iq2(i:i + (jump-1)), round(6371e3*pi/fs), 'coeff');
        [~,idx] = max(abs(r));
        iqcorrelate12(data) = lags(idx);
        data = data + 1;
    end
    meanCorr = mean(iqcorrelate12);
end