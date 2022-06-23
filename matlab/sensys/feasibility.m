clear; close all;

c = 299792458; % m/s

load('results/jitter');
add_jitter = 1;
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [0 0 8.86 7.8]);

%% Plot GSM different values
load('results/splitter_gsm_lte_original_dphase');
gsm_orig = abs(doa_meters_2) + 150 + add_jitter * jit50; gsm_orig(gsm_orig > 1000) = 1000;
gsm_mean = abs(doa_meters) + 150 + add_jitter * jit50; gsm_mean(gsm_mean > 1000) = 1000;
load('results/splitter_gsm_lte_fo_correction_no_interp_dphase');
gsm_noin = abs(doa_meters_2) + 150 + add_jitter * jit50; gsm_noin(gsm_noin > 1000) = 1000;
load('results/splitter_gsm_lte_fo_correction_10_interp_dphase');
gsm_10in = abs(doa_meters_2) + 15 + add_jitter * jit50; gsm_10in(gsm_10in > 1000) = 1000;

data = [gsm_orig; gsm_mean; gsm_noin; gsm_10in]/c*1e6;
gsm = [ones(50,1); 2*ones(50,1); 3*ones(50,1); 4*ones(50,1)];
gsms = array2table([data,gsm],'VariableNames',{'data','type'});
names = ["orig", "mean", "foc", "ups"];
namedNames= categorical(gsms.type,1:4,names);

figure(); grid on; hold on; axis square;
% boxplot([gsm_orig, gsm_mean, gsm_noin, gsm_10in]/c*1e6, 'labels', {'orig', 'mean', 'foc', 'ups'});
boxchart(namedNames,gsms.data,'Notch','on');
ylabel('TDOA (μs)'); ylim([0,3]); yticks([1,2]);
box on;
set(gca,'FontSize',40,'XTickLabelRotation',0);

%%  Plot LTE
load('results/splitter_lte_lte_original_dphase');
lte_orig = abs(doa_meters_2) + 150 + add_jitter * jit600;
lte_mean = abs(doa_meters) + 150 + add_jitter * jit600;
load('results/splitter_lte_lte_fo_correction_no_interp_dphase');
lte_noin = abs(doa_meters_2) + 150 + add_jitter * jit600;
load('results/splitter_lte_lte_fo_correction_10_interp_dphase');
lte_10in = abs(doa_meters_2) + 15 + add_jitter * jit600;

figure(); grid on; hold on; axis square;
boxplot([lte_orig, lte_mean, lte_noin, lte_10in]/c * 1e6, 'labels', {'orig', 'mean', 'foc', 'ups'});
ylabel('TDOA (μs)'); ylim([0,3]);
lte_10in = randsample(lte_10in,50);

%% Plot DVB-T different values
load('results/splitter_dvbt_lte_original_dphase');
dvb_orig = abs(doa_meters_2) + 150 + add_jitter * jit50; dvb_orig(dvb_orig > 1000) = 1000;
dvb_mean = abs(doa_meters) + 150 + add_jitter * jit50; dvb_mean(dvb_mean > 1000) = 1000;
load('results/splitter_dvbt_lte_fo_correction_no_interp_dphase');
dvb_noin = abs(doa_meters_2) + 150 + add_jitter * jit50; dvb_noin(dvb_noin > 1000) = 1000;
load('results/splitter_dvbt_lte_fo_correction_10_interp_dphase');
dvb_10in = abs(doa_meters_2) + 15 + add_jitter * jit50; dvb_10in(dvb_10in > 1000) = 1000;

data = [dvb_orig; dvb_mean; dvb_noin; dvb_10in]/c*1e6;
dvb = [ones(50,1); 2*ones(50,1); 3*ones(50,1); 4*ones(50,1)];
dvbs = array2table([data,dvb],'VariableNames',{'data','type'});
names = ["orig", "mean", "foc", "ups"];
namedNames= categorical(dvbs.type,1:4,names);

figure(); grid on; hold on; axis square;
% boxplot([dvb_orig, dvb_mean, dvb_noin, dvb_10in]/c * 1e6, 'labels', {'orig', 'mean', 'foc', 'ups'});
boxchart(namedNames,dvbs.data,'Notch','on');
ylabel(['TDOA (μs)']); ylim([0,3]); yticks([1,2]);
box on;
set(gca,'FontSize',40,'XTickLabelRotation',0)

%% Get all methods
load('results/splitter_fm_lte_fo_correction_10_interp_dphase.mat');
fm_10in = abs(randsample(doa_meters_2,50)) + 15 + add_jitter * jit50;
load('results/splitter_dab_lte_fo_correction_10_interp_dphase.mat')
dab_10in = abs(doa_meters_2) + 15 + add_jitter * jit50;

data = [lte_10in; gsm_10in; dvb_10in; dab_10in; fm_10in]/c*1e6;
type = [ones(50,1); 2*ones(50,1); 3*ones(50,1); 4*ones(50,1); 5*ones(50,1)];
techs = array2table([data,type],'VariableNames',{'data','type'});
names = ["lte","gsm","dvbt","dab","fm"];
namedNames= categorical(techs.type,1:5,names);

figure(); grid on; hold on; axis square
% boxplot([lte_10in, gsm_10in, dvb_10in, dab_10in, fm_10in]/c*1e6,'labels',{'lte','gsm','dvb-t', 'dab', 'fm'});
boxchart(namedNames, techs.data,'Notch','on');
ylim([0,0.3]);
ylabel(['TDOA (μs)']); ylim([0,0.3]); yticks([0.1,0.2]);
box on;
set(gca,'FontSize',40,'XTickLabelRotation',0);

%% Compare different corr methods
load('results/562806_2000000_fo_correction_10_interp_abs.mat');
dvb_abs = abs(doa_meters_2) + 15 + add_jitter * jit600;
dvb_abs = randsample(dvb_abs,50);
load('results/562806_2000000_fo_correction_10_interp_iq.mat');
dvb_iq = abs(doa_meters_2) + 15 + add_jitter * jit600;
dvb_iq(dvb_iq > 400) = 400;
dvb_iq = randsample(dvb_iq,50);

data = [dvb_10in; dvb_abs; dvb_iq]/c*1e6;
method = [ones(50,1); 2*ones(50,1); 3*ones(50,1)];
methods = array2table([data,method],'VariableNames',{'data','type'});
names = ["dphase","abs","iq"];
namedNames= categorical(methods.type,1:3,names);

figure(); grid on; hold on; axis square;
% boxplot([dvb_10in, dvb_abs, dvb_iq]/c*1e6,'labels',{'dphase','abs','iq'});
boxchart(namedNames, methods.data,'Notch','on');
ylabel(['TDOA (μs)']); ylim([0,0.3]); yticks([0.1,0.2]);
box on;
set(gca,'FontSize',40,'XTickLabelRotation',0);
%% Subplot
x0=10;
y0=10;
width=550;
height=400;

subplot(1,4,1);
axis equal;
boxplot([gsm_orig, gsm_mean, gsm_noin, gsm_10in]/c * 1e6, 'labels', {'Orig', 'Mean', 'FOC', 'Ups x10'});
ylabel('TDOA (\mus)'); ylim([0,3]); yticks([1,2])
set(gca,'FontSize',30,'XTickLabelRotation',0);
set(gcf,'position',[x0,y0,width,height]);

subplot(1,4,2);
axis equal;
boxplot([dvb_orig, dvb_mean, dvb_noin, dvb_10in]/c * 1e6, 'labels', {'Orig', 'Mean', 'FOC', 'Ups x10'});
ylabel(['TDOA (\mus)']); ylim([0,3]); yticks([1,2])
set(gca,'FontSize',30,'XTickLabelRotation',0);
set(gcf,'position',[x0,y0,width,height]);

subplot(1,4,3);
axis equal;
boxplot([dvb_10in, dvb_abs, dvb_iq]/c*1e6,'labels',{'dphase','abs','iq'});
ylabel(['TDOA (\mus)']); ylim([0,0.3]); yticks([0.1,0.2])
set(gca,'FontSize',30,'XTickLabelRotation',0);
set(gcf,'position',[x0,y0,width,height]);

subplot(1,4,4);
axis equal;
boxplot([lte_10in, gsm_10in, dvb_10in, dab_10in, fm_10in]/c*1e6,'labels',{'LTE','GSM','DVB-T', 'DAB', 'FM'});
ylabel(['TDOA (\mus)']); ylim([0,0.3]); yticks([0.1,0.2]);
set(gca,'FontSize',30,'XTickLabelRotation',0);
set(gcf,'position',[x0,y0,width,height]);