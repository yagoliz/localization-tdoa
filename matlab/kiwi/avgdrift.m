%% In this script we'll analyze the drift among pairs of sensors
% Selected sensors (in order):
% - DL0ABT
% - S57BIT
% - F1JEK
clear; close all;
% Average drift in microseconds
drift = zeros(3,30); % We have 30 different campaigns and 3 pairs

plotdrift = false;
for i = 1:30
    filepath = ['data/dcf77_', num2str(i), '.mat'];
    load(filepath);
    
    maxLength = min([length(input(1).t), length(input(2).t), length(input(3).t)]);
    fs = input(1).fs;
    
    % Get drifts
    td12 = input(1).t(1:maxLength) - input(2).t(1:maxLength); td12 = td12 - td12(1);
    td13 = input(1).t(1:maxLength) - input(3).t(1:maxLength); td13 = td13 - td13(1);
    td23 = input(2).t(1:maxLength) - input(3).t(1:maxLength); td23 = td23 - td23(1);
    
    % We introduce these drifts in the drift vector
    drift(1,i) = (td12(1) - td12(end))*1e6/maxLength*fs;
    drift(2,i) = (td13(1) - td13(end))*1e6/maxLength*fs;
    drift(3,i) = (td23(1) - td23(end))*1e6/maxLength*fs;
    
    % Plot
    if plotdrift
        hold on; grid on;
        x = (0:maxLength-1)/12000;
        plot(x, td12); plot(x, td13); plot(x, td23);
        xlim([x(1), x(end)]);
        xlabel('Time (s)'); ylabel('Drift (\mus)');
        legend('Sensors 1 & 2', 'Sensors 1 & 3', 'Sensors 2 & 3');
    end
end

% Plot the drifts
figure();
hold on; grid on; 
boxplot(abs(drift'), 'labels', {'Sensors 1 & 2', 'Sensors 1 & 3', 'Sensors 2 & 3'}); 
title('Drift between sensors per second'); ylabel('Drift per second (\mus/s)');