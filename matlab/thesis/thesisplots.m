clear; close all;

%% ECDF
load('results_lte_analysis');
[F1,X1] = ecdf(sorted_cell1);
[F2,X2] = ecdf(sorted_cell2);

figure(); grid on; hold on; stairs(X1,F1,'LineWidth',2); stairs(X2,F2,'LineWidth',2);
xlabel('Max distance between points with same Cell ID (m)'); ylabel('Cumulative Probability'); legend('San Sebastian', 'Alcorcon','Location','SouthEast');

%% Plot GSM different values
load('jitter');
add_jitter = 1;

load('results/939806_2000000_original_dphase');
gsm_orig = abs(doa_meters_2) + 150 + add_jitter * jitter;
gsm_mean = abs(doa_meters) + 150 + add_jitter * jitter;
load('results/939806_2000000_fo_correction_no_interp_dphase');
gsm_noin = abs(doa_meters_2) + 150 + add_jitter * jitter;
load('results/939806_2000000_fo_correction_10_interp_dphase');
gsm_10in = abs(doa_meters_2) + 15 + add_jitter * jitter;

figure(); grid on; hold on;
boxplot([gsm_orig, gsm_mean, gsm_noin, gsm_10in], 'labels', {'Original', 'RS Mean', 'FOC', 'FOC + Upsampling x10'});
ylabel('TDOA (m)');

%% Plot DVB-T different values
load('results/562806_2000000_original_dphase');
dvb_orig = abs(doa_meters_2) + 150 + add_jitter * jitter;
dvb_mean = abs(doa_meters) + 150 + add_jitter * jitter;
load('results/562806_2000000_fo_correction_no_interp_dphase');
dvb_noin = abs(doa_meters_2) + 150 + add_jitter * jitter;
load('results/562806_2000000_fo_correction_10_interp_dphase');
dvb_10in = abs(doa_meters_2) + 15 + add_jitter * jitter;

figure(); grid on; hold on;
boxplot([dvb_orig, dvb_mean, dvb_noin, dvb_10in], 'labels', {'Original', 'RS Mean', 'FOC', 'FOC + Upsampling x10'});
ylabel('TDOA (m)');

%% Get LTE data
load('results/806806_2000000_fo_correction_10_interp_dphase');
lte_10in = abs(doa_meters_2) + 15 + add_jitter * jitter;
load('results/99806_2000000_fo_correction_10_interp_dphase');
fm_10in = abs(doa_meters_2) + 15 + add_jitter * jitter;

figure(); grid on; hold on;
boxplot([lte_10in, gsm_10in, dvb_10in, fm_10in],'labels',{'LTE','GSM','DVB-T','FM'});
ylabel('TDOA (m)');

%% Compare different corr methods
load('results/562806_2000000_fo_correction_10_interp_abs');
dvb_abs = abs(doa_meters_2) + 15 + add_jitter * jitter;
load('results/562806_2000000_fo_correction_10_interp_iq');
dvb_iq = abs(doa_meters_2) + 15 + add_jitter * jitter;
dvb_iq(dvb_iq > 400) = 400;

figure(); grid on; hold on;
boxplot([dvb_10in, dvb_abs, dvb_iq],'labels',{'dphase','abs','iq'});
ylabel('TDOA (m)'); ylim([0,150]);