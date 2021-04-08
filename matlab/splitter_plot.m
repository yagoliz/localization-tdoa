clear; close all;

%% Plot DAB different values
load('results/splitter_dab_original_dphase.mat');
dab_orig = abs(doa_meters_2) + 150; 
dab_mean = abs(doa_meters) + 150;
load('results/splitter_dab_fo_correction_no_interp_dphase');
dab_noin = abs(doa_meters_2) + 150;
load('results/splitter_dab_fo_correction_10_interp_dphase');
dab_10in = abs(doa_meters_2) + 15;

% DAB
indicesdab = min(correlation_value(1,:,:),[],2) > 0.1;
dab_corr_12 = squeeze(min(correlation_value(1,:,1:51), [], 2));

figure(); grid on; hold on;
title('DAB');
boxplot([dab_orig(indicesdab), dab_mean(indicesdab), dab_noin(indicesdab), dab_10in(indicesdab)], 'labels', {'Original', 'RS Mean', 'FOC', 'FOC + Upsampling x10'});
ylabel('TDOA (m)'); ylim([0,1000]);

%% Plot DVB-T different values
load('results/splitter_dvbt_original_dphase');
dvb_orig = abs(doa_meters_2) + 150;
dvb_mean = abs(doa_meters) + 150;
load('results/splitter_dvbt_fo_correction_no_interp_dphase');
dvb_noin = abs(doa_meters_2) + 150;
load('results/splitter_dvbt_fo_correction_10_interp_dphase');
dvb_10in = abs(doa_meters_2) + 15;

% DVB-T
indicesdvb = min(correlation_value,[],2) > 0.1;
dvb_corr_12 = min(correlation_value, [], 2);


figure(); grid on; hold on;
title('DVB-T');
boxplot([dvb_orig(indicesdvb), dvb_mean(indicesdvb), dvb_noin(indicesdvb), dvb_10in(indicesdvb)], 'labels', {'Original', 'RS Mean', 'FOC', 'FOC + Upsampling x10'});
ylabel('TDOA (m)'); ylim([0,1000]);

%% Plot GSM different values
load('results/splitter_gsm_original_dphase');
gsm_orig = abs(doa_meters_2) + 150;
gsm_mean = abs(doa_meters) + 150;
load('results/splitter_gsm_fo_correction_no_interp_dphase');
gsm_noin = abs(doa_meters_2) + 150;
load('results/splitter_gsm_fo_correction_10_interp_dphase');
gsm_10in = abs(doa_meters_2) + 15;

indicesgsm = min(correlation_value,[],2) > 0.1;
gsm_corr_12 = min(correlation_value, [], 2);


figure(); grid on; hold on;
title('GSM');
boxplot([gsm_orig(indicesgsm), gsm_mean(indicesgsm), gsm_noin(indicesgsm), gsm_10in(indicesgsm)], 'labels', {'Original', 'RS Mean', 'FOC', 'FOC + Upsampling x10'});
ylabel('TDOA (m)'); ylim([0,1000]);

%% Compare different DAB and DVB-T
figure(); grid on; hold on;
g = [zeros(sum(indicesdab),1); ones(sum(indicesdvb),1); 2*ones(sum(indicesgsm),1)];
boxplot([dab_10in(indicesdab); dvb_10in(indicesdvb); gsm_10in(indicesgsm)], g,'labels',{'DAB','DVB-T','GSM'});
ylabel('TDOA (m)'); ylim([0,150]);

%% Correlation vs Result
figure(); grid on; hold on;
scatter(dab_corr_12(indicesdab(1:51)), dab_noin(indicesdab(1:51)));
scatter(dvb_corr_12(indicesdvb), dvb_noin(indicesdvb));
scatter(gsm_corr_12(indicesgsm), gsm_noin(indicesgsm));
legend('DAB','DVB-T','GSM');
xlabel('Correlation Value'); ylabel('TDOA Error (m)');
ylim([0,500]);
% DAB
% figure(); grid on; hold on;
% scatter(dab_corr_12(:), dab_noin(1:3:end));
% scatter(dab_corr_13(:), dab_noin(2:3:end));
% scatter(dab_corr_23(:), dab_noin(3:3:end));
% legend('1&2','1&3','2&3'); title('DAB');
% xlabel('Correlation value'); ylabel('TDOA Error');
% 
% % DVB-T
% figure(); grid on; hold on;
% scatter(dvb_corr_12(:), dvb_noin(1:3:end));
% scatter(dvb_corr_13(:), dvb_noin(2:3:end));
% scatter(dvb_corr_23(:), dvb_noin(3:3:end));
% legend('1&2','1&3','2&3'); title('DVB-T');
% xlabel('Correlation value'); ylabel('TDOA Error');
