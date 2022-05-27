% =========================================================================
%  Results for Madrid Localization
%  Author: Yago Lizarribar
% =========================================================================
clear all; close all; clc;
warning ('off','all');

% Speed of light
global c; 
c = 299792458; % m/s

% adds subfolder with functions to PATH
[p,~,~] = fileparts(mfilename('fullpath'));
addpath([p, '../']);
addpath([p, '../tdoa']);
addpath([p, '../ltess']);
addpath([p, '../ltess/lte/']);
addpath([p, '../geodesy']);
addpath([p, '../optimization']);


% Load main data
NUM_SENSORS=3;

problem_data = json_load(sprintf("madrid_s%i.json",NUM_SENSORS));
% 6 sensors
folders = ["20211028", "20211029", "20211129", "20211130", "20211214"]; % 6 & 3 (good) sensors
files = [5, 5, 5, 5, 5];

% 5 sensors
% folders = ["20211024", "20211025", "20211027", "20211028", "20211029", "20211129", "20211130", "20211214"];
% files = [5, 5, 5, 5, 5, 5, 5, 5];

% 4 sensors
% folders = ["20211012", "20211024", "20211025", "20211027" "20211028"];
% files = [17,5,5,5,5];

corr_methods = ["abs", "dphase", "iq"];
interp = [1, 5, 10, 20, 50];
correct = [true, false];
lte = [true, false];

lls = zeros(sum(files),2);
error_lls = zeros(sum(files),1);
nlls = zeros(sum(files),3);
error_nlls = zeros(sum(files),1);

% Main loop
for corr_id = 1:length(corr_methods)
    problem_data.config.corr_method = corr_methods(corr_id);
    for interp_id = 1:length(interp)
        problem_data.config.interp = interp(interp_id);
        for ii = 1:length(correct)
            problem_data.config.correct = correct(ii);
            for jj = 1:length(lte)
            problem_data.config.use_lte = lte(jj);
                processed = 1;
                for foldernum = 1:length(folders)
                    problem_data.config.folder_date = folders(foldernum);
                    for filenum = 0:files(foldernum)-1
                        problem_data.config.filenum = filenum;

                        % Localization routine
                        res = tdoa_localization(problem_data);
                        lls(processed+filenum,:) = res.res_linear';
                        error_lls(processed+filenum) = res.error_linear;
                        nlls(processed+filenum,:) = res.res_accurate';
                        error_nlls(processed+filenum) = res.error_nonlin;
                    end % filenum
                    processed = processed + filenum + 1;
                end % folders
                
                % Saving our beloved file
                filename = sprintf("sensys/madrid/s%i/m%s_i%i_c%i_l%i.mat", NUM_SENSORS, corr_methods(corr_id), interp(interp_id), correct(ii), lte(jj));
                save(filename,'lls', 'error_lls', 'nlls', 'error_nlls');
            end % lte
        end % correct
    end % interp
end % corr_methods
