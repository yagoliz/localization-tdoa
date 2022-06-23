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
problem_data = json_load(sprintf("zurich.json"));

% 4 sensors
folders = ["20220612", "20220613", "20220614", "20220615" "20220616","20220617"];
% folders = ["20220615"];
files = [5,10,10,10,10,10];
% files = [3];

corr_methods = ["dphase"];
interp = [1];
correct = [true];
lte = [true];
multipath = false;

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
                for kk = 1:length(multipath)
                    problem_data.config.correct_multipath = multipath(kk);
                    processed = 1;
                    for foldernum = 1:length(folders)
                        problem_data.config.folder_date = folders(foldernum);
                        for filenum = 0:files(foldernum)-1
                            problem_data.config.filenum = filenum;
    
                            % Localization routine
                            res = tdoa_localization_zurich(problem_data);
                            lls(processed+filenum,:) = res.res_linear';
                            error_lls(processed+filenum) = res.error_linear;
                            nlls(processed+filenum,:) = res.res_accurate';
                            error_nlls(processed+filenum) = res.error_nonlin;
                        end % filenum
                        processed = processed + filenum + 1;
                    end % folders
                    if multipath(kk)
                        % Saving our beloved file
                        filename = sprintf("sensys/zurich/m%s_i%i_c%i_l%i_.mat", corr_methods(corr_id), interp(interp_id), correct(ii), lte(jj));
                        save(filename,'lls', 'error_lls', 'nlls', 'error_nlls');
                    else
                        % Saving our beloved file
                        filename = sprintf("sensys/zurich/m%s_i%i_c%i_l%i.mat", corr_methods(corr_id), interp(interp_id), correct(ii), lte(jj));
                        save(filename,'lls', 'error_lls', 'nlls', 'error_nlls');
                    end
                end
            end % lte
        end % correct
    end % interp
end % corr_methods
