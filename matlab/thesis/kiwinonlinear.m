% =========================================================================
%  TDOA Localization with KiwiSDR data
%  Author: Yago Lizarribar
% =========================================================================

clear all;
clc;
close all;

global sensors doa_meters2 combinations R;

% adds subfolder with functions to PATH
[p,n,e] = fileparts(mfilename('fullpath'));
addpath([p '/tdoa']);
addpath([p '/kiwi/m/']);

% Data

%% Config
tx_type = 'dcf2';

% DCF77
if strcmp(tx_type, 'dcf')
    load('kiwi/data/dcf77_1_pre.mat');
    transmitter = [50.01556, 9.01083];

elseif strcmp(tx_type, 'dcf2')
    load('kiwi/data/dcf2_1.mat');
    transmitter = [50.01556, 9.01083];

elseif strcmp(tx_type, 'bbc')
    load('kiwi/data/bbc_1.mat');
    transmitter = [52.296667, -2.105278];
    
elseif strcmp(tx_type, 'bbc2')
    load('kiwi/data/bbc2_1.mat');
    transmitter = [52.296667, -2.105278];
    
elseif strcmp(tx_type, 'tdf')
    load('kiwi/data/tdf_1.mat');
    transmitter = [47.1695, 2.2046];
end

   
%% Constants
% Earth Radius
R = 6371 * 1e3;
c = 299792.458 * 1e3;
       

% Get distances from transmitter to sensor
NUM_SENSORS = 5;

if NUM_SENSORS > numel(input)
    error('There are not as many sensors!');
end

%% Sensors might be disordered
sensors = zeros(NUM_SENSORS,2);
for jj = 1:NUM_SENSORS
    sensors(jj,:) = input(jj).coord;
end    


%% Optimization parameters
optionsfmin = optimoptions('fminunc','Algorithm','quasi-newton');
X0 = [mean(sensors(:,1)), mean(sensors(:,2))];

% We can use non linear least squares
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');

%% Loop
files = dir(['kiwi/data/', tx_type, '_*']);
numfiles = length(files);
errors = zeros(numfiles,1);
errorsfmin = zeros(numfiles,1);

% Combinations
combinations = nchoosek(1:NUM_SENSORS,2);
doa = zeros(size(combinations,1),1);
doa_meters = doa;
doa_meters2 = doa;
iqcorrelate = doa;
corrfactor = doa;
correlations = zeros(numfiles,size(combinations,1));

for ii = 1:numfiles
    fprintf('----------------------------------------------------------\n');
    fprintf('Loop: %d\n', ii);
    fprintf('----------------------------------------------------------\n');
    
    %% Load data
    load([files(ii).folder, filesep, files(ii).name]);
    [input, minLength] = alignInputs(input);

    %% Read the IQ
    iqsignal = zeros(minLength,NUM_SENSORS);
    t = zeros(minLength, NUM_SENSORS);
    fs = zeros(1, NUM_SENSORS);
    
    for jj = 1:NUM_SENSORS
        iqsignal(:,jj) = input(jj).z(1:minLength);
        fs(:,jj) = input(jj).fs;
        t(:,jj) = input(jj).t(1:minLength);
    end
    
    fsdist = dist(fs);
    if any(fsdist(fsdist > 10))
        error('Sampling frequencies are not the same');
    end
    
    %% Sensors might be disordered
    sensors = zeros(NUM_SENSORS,2);
    for jj = 1:NUM_SENSORS
        sensors(jj,:) = input(jj).coord;
    end    

    %% Process tdoa for sensor pairs
    
    for jj = 1:size(combinations,1)
        si = combinations(jj,1);
        sj = combinations(jj,2);
        [doa(jj), iqcorrelate, corrfactor, dt] = tdoa_kiwi(iqsignal(:,si), t(:,si),...
                                                           iqsignal(:,sj), t(:,sj),...
                                                           mean([fs(si),fs(sj)]), 1, 'iq');
        
        correlations(ii,jj) = doa(jj)/12;
        doa_meters(jj) = (doa(jj)/mean([fs(si),fs(sj)])) * c;
        doa_meters2(jj) = (mean(dt + iqcorrelate/mean([fs(si),fs(sj)]))) * c;
        
    end
    
    % Let's do non linear optimization
    [Xfmin,fval,exitflag,outputmin] = fminunc(@optim,X0,optionsfmin);
    X = lsqnonlin(@optimnlls,X0,[],[],options);
    
    
    angle = distance(X(1), X(2), transmitter(1), transmitter(2));
    anglefmin = distance(Xfmin(1), Xfmin(2), transmitter(1), transmitter(2));
    
    errors(ii) = deg2km(angle);
    errorsfmin(ii) = deg2km(anglefmin);
end

%% Alignment function
function [inputAli, minLength] = alignInputs(input)
    numelem = numel(input);
    startt = zeros(numelem,1);
    endt = zeros(numelem,1);
    for ii = 1:numelem
        startt(ii) = input(ii).t(1);
        endt(ii) = input(ii).t(end);
        
    end
    t0 max(startt);
    t1 = min(endt);

    minLength = inf;
    for ii=1:numelem
        if ~input(ii).use
          continue
        end
        b = input(ii).t<t0 | input(ii).t>t1;
        input(ii).t(b) = [];
        input(ii).z(b) = [];
        input(ii).use  = numel(input(ii).z)/input(ii).fs > 10;
        if ~input(ii).use
          printf('tdoa_read_data: %-40s excluded (%.2f sec < %g sec overlap)\n', ...
                 input(ii).fn, numel(input(ii).z)/input(ii).fs, 10);
        end
        
        if length(input(ii).t) < minLength
            minLength = length(input(ii).t);
        end
        
        if abs(input(ii).fs - 12e3) > 1e3
            input(ii).fs = numel(input(ii).t)/(input(ii).t(end) - input(ii).t(1));
        end
    end

    inputAli = input;
  
end

%% BFGS approach
function F = optim(X)
    global sensors doa_meters2 combinations;
    lat = X(1);
    lon = X(2);
    
    d = zeros(size(sensors,1),1);
    for ii = 1:size(sensors,1)
        d(ii) = deg2km(distance(lat, lon, sensors(ii,1), sensors(ii,2)));
    end
    
    t = zeros(size(combinations,1),1);
    for ii = 1:size(t,1)
        si = combinations(ii,1);
        sj = combinations(ii,2);
        t(ii) = (d(si) - d(sj))*1000;
    end
    
    F = 0.5 * sum((doa_meters2 - t).^2);
end


%% Non linear LeastSquares function
function F = optimnlls(X)
    global sensors doa_meters2 combinations;
    lat = X(1);
    lon = X(2);
    
    d = zeros(size(sensors,1),1);
    for ii = 1:size(sensors,1)
        d(ii) = deg2km(distance(lat, lon, sensors(ii,1), sensors(ii,2)));
    end
    
    t = zeros(size(combinations,1),1);
    for ii = 1:size(t,1)
        si = combinations(ii,1);
        sj = combinations(ii,2);
        t(ii) = (d(si) - d(sj))*1000;
    end
    
    F = (doa_meters2 - t);
end
