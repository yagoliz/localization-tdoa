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
NUM_SENSORS=6;

problem_data = json_load(sprintf("madrid_s%i.json",NUM_SENSORS));

% 6 sensors
folders = ["20211028", "20211029", "20211129", "20211130", "20211214"]; % 6 & 3 (good) sensors
files = [5, 5, 5, 5, 5];

corr_methods = ["dphase"];
interp = [20];
correct = [true];
lte = [true];
multipath = [false,true];
% samples = [50000, 100000, 200000, 300000, 400000, 500000, 600000];
samples = [100000, 200000, 500000];

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
                for ll = 1:length(samples)
                    problem_data.config.samples_slice = samples(ll);
                    for kk = 1:length(multipath)
                        problem_data.config.correct_multipath = multipath(kk);
                        processed = 1;
                        for foldernum = 1:length(folders)
                            problem_data.config.folder_date = folders(foldernum);
                            for filenum = 0:files(foldernum)-1
                                problem_data.config.filenum = filenum;
        
                                % Localization routine
                                res = tdoa_localization_samples(problem_data);
                                lls(processed+filenum,:) = res.res_linear';
                                error_lls(processed+filenum) = res.error_linear;
                                nlls(processed+filenum,:) = res.res_accurate';
                                error_nlls(processed+filenum) = res.error_nonlin;
                            end % filenum
                            processed = processed + filenum + 1;
                            fprintf("- Processed: %i\n", processed);
                        end % folders
                        fprintf("Saving...\n");
                        if multipath(kk)
                            % Saving our beloved file
                            filename = sprintf("sensys/madrid2/samples/s%i_m%s_i%i_c%i_l%i_n%i_mp.mat", NUM_SENSORS, corr_methods(corr_id), interp(interp_id), correct(ii), lte(jj), samples(ll));
                            save(filename,'lls', 'error_lls', 'nlls', 'error_nlls');
                        else
                            % Saving our beloved file
                            filename = sprintf("sensys/madrid2/samples/s%i_m%s_i%i_c%i_l%i_n%i.mat", NUM_SENSORS, corr_methods(corr_id), interp(interp_id), correct(ii), lte(jj), samples(ll));
                            save(filename,'lls', 'error_lls', 'nlls', 'error_nlls');
                        end
                    end
                end
            end % lte
        end % correct
    end % interp
end % corr_methods

