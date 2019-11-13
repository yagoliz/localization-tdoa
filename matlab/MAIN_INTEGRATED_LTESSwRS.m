%% STEP 0 -- INIT

clear
close all

saveall = 0;

N_receiver = 3;
Test_for_dfo = 3;
Test_for_y1y2 = 10;

fRS = 806e6;
fUS = 244e6;

fs_y1y2 = 2.048e6;
% SAMPLING_RATE_DFO = 1.92e6;
SAMPLING_RATE_DFO = fs_y1y2;
num_samples_y1y2 = 3e6;

addpath('./common/')

RESAMPLE_FACTOR = 60;
PSS_STEP = 9600;
SEARCH_WINDOW = 150;
CORRELATION_FACTOR = 0.1;
CORRELATION_FACTOR_PREAMBLE = 0.1;
PREAMBLE=20; % Number of analyzed PSS before start to jump
FLIP=0; % 0 disabled 1: enabled  (testing purposes)
POLYNOMIAL_DEGREE=1;

[Z, Z_t] = get_Zadoof();

%% STEP 2 -- SIGNAL ANALYSIS

addpath('./functions/')

alignment = true;

Ts = 1/fs_y1y2;
t = 0:1/fs_y1y2:((1/fs_y1y2)*(num_samples_y1y2-1));
alpha = exp(-Ts/0.5e-6);
r_filt_b = 1-alpha;
r_filt_a = [1, -alpha];
r_filt_z = [];

analyzed_devices = nchoosek(1:N_receiver,2);

delay_time_RS = zeros(size(analyzed_devices,1),Test_for_y1y2);
delay_time_US = zeros(size(analyzed_devices,1),Test_for_y1y2);
delay_time_RS_ceck = zeros(size(analyzed_devices,1),Test_for_y1y2);

for index1 = 1:size(analyzed_devices,1)
    for index2 = 1:Test_for_y1y2
        
        fprintf('index1:%d --- index2: %d\n',index1,index2)
        filename_for_corr_1 = sprintf('RTL_test1_LTE244LTE\\RTL_D%d_T%d.dat',analyzed_devices(index1,1),index2);
        filename_for_corr_2 = sprintf('RTL_test1_LTE244LTE\\RTL_D%d_T%d.dat',analyzed_devices(index1,2),index2);
        
        fileID1 = fopen(filename_for_corr_1);
        fileID2 = fopen(filename_for_corr_2);
        A = fread(fileID1);
        B = fread(fileID2);
        fclose(fileID1);
        fclose(fileID2);
        
        inphase1 = A(1:2:end) -128;
        quadrature1 = A(2:2:end) -128;
        inphase2 = B(1:2:end) -128;
        quadrature2 = B(2:2:end) -128;
        
        r1 = (inphase1 + 1i*quadrature1);
        [r1_filt_B,~] = filter(r_filt_b,r_filt_a,r1,r_filt_z);
        r1_filt_B = r1_filt_B - mean(r1_filt_B);
        r2 = (inphase2 + 1i*quadrature2);
        [r2_filt_B,r_filt_z_A] = filter(r_filt_b,r_filt_a,r2,r_filt_z);
        r2_filt_B = r2_filt_B - mean(r2_filt_B);
        
        [PPM_1,~,~,~,~,~] = getDrift(r1(1:num_samples_y1y2/3)./128,SAMPLING_RATE_DFO, Z, Z_t, PREAMBLE, ...
            PSS_STEP, SEARCH_WINDOW, CORRELATION_FACTOR, RESAMPLE_FACTOR, FLIP, POLYNOMIAL_DEGREE);
        T_s = 1/SAMPLING_RATE_DFO;
        delta_f=(PPM_1*1e-6)*806e6;
        Z_t_rotated = {};
        Z_t_rotated{1} = Z_t{1}.*exp(-1i*2*pi*T_s*delta_f*(1:length(Z_t{1}))');
        Z_t_rotated{2} = Z_t{2}.*exp(-1i*2*pi*T_s*delta_f*(1:length(Z_t{2}))');
        Z_t_rotated{3} = Z_t{3}.*exp(-1i*2*pi*T_s*delta_f*(1:length(Z_t{3}))');
        [PPM_1,~,~,~,~,~,~,~,~,~] = getDrift(r1(1:num_samples_y1y2/3)./128,SAMPLING_RATE_DFO, Z, Z_t_rotated, ...
            PREAMBLE, PSS_STEP, SEARCH_WINDOW, CORRELATION_FACTOR, RESAMPLE_FACTOR, FLIP, POLYNOMIAL_DEGREE);

        [PPM_2,~,~,~,~,~] = getDrift(r2(1:num_samples_y1y2/3)./128,SAMPLING_RATE_DFO, Z, Z_t, PREAMBLE, ...
            PSS_STEP, SEARCH_WINDOW, CORRELATION_FACTOR, RESAMPLE_FACTOR, FLIP, POLYNOMIAL_DEGREE);
        T_s = 1/SAMPLING_RATE_DFO;
        delta_f=(PPM_2*1e-6)*806e6;
        Z_t_rotated = {};
        Z_t_rotated{1} = Z_t{1}.*exp(-1i*2*pi*T_s*delta_f*(1:length(Z_t{1}))');
        Z_t_rotated{2} = Z_t{2}.*exp(-1i*2*pi*T_s*delta_f*(1:length(Z_t{2}))');
        Z_t_rotated{3} = Z_t{3}.*exp(-1i*2*pi*T_s*delta_f*(1:length(Z_t{3}))');
        [PPM_2,~,~,~,~,~,~,~,~,~] = getDrift(r2(1:num_samples_y1y2/3)./128,SAMPLING_RATE_DFO, Z, Z_t_rotated, ...
            PREAMBLE, PSS_STEP, SEARCH_WINDOW, CORRELATION_FACTOR, RESAMPLE_FACTOR, FLIP, POLYNOMIAL_DEGREE);
        
        index_samples = [0.1e6 0.9e6 1.1e6 1.9e6 2.1e6 2.9e6];
        
        % FO correction
        r1_filt_B(1:num_samples_y1y2/3) = r1_filt_B(1:num_samples_y1y2/3).*(exp(-1i*2*pi*t(1:num_samples_y1y2/3)*(-PPM_1)*(1e-6)*fUS).');
        r2_filt_B(1:num_samples_y1y2/3) = r2_filt_B(1:num_samples_y1y2/3).*(exp(-1i*2*pi*t(1:num_samples_y1y2/3)*(-PPM_2)*(1e-6)*fUS).');
        r1_filt_B((1+num_samples_y1y2/3):(2*num_samples_y1y2/3)) = r1_filt_B((1+num_samples_y1y2/3):(2*num_samples_y1y2/3)).*(exp(-1i*2*pi*t((1+num_samples_y1y2/3):(2*num_samples_y1y2/3))*(-PPM_1)*(1e-6)*fRS).');
        r2_filt_B((1+num_samples_y1y2/3):(2*num_samples_y1y2/3)) = r2_filt_B((1+num_samples_y1y2/3):(2*num_samples_y1y2/3)).*(exp(-1i*2*pi*t((1+num_samples_y1y2/3):(2*num_samples_y1y2/3))*(-PPM_2)*(1e-6)*fRS).');
        r1_filt_B((1+2*num_samples_y1y2/3):end) = r1_filt_B((1+2*num_samples_y1y2/3):end).*(exp(-1i*2*pi*t((1+2*num_samples_y1y2/3):end)*(-PPM_1)*(1e-6)*fUS).');
        r2_filt_B((1+2*num_samples_y1y2/3):end) = r2_filt_B((1+2*num_samples_y1y2/3):end).*(exp(-1i*2*pi*t((1+2*num_samples_y1y2/3):end)*(-PPM_2)*(1e-6)*fUS).');
        
        % time correction
% % % % % % %         r1_filt_B = interp(r1_filt_B,4);
% % % % % % %         r2_filt_B = interp(r2_filt_B,4);
        Sampling_vector = 1:numel(r1_filt_B);
        r1_filt_B = interp1(Sampling_vector,r1_filt_B,Sampling_vector*(1+((1e-6)*(PPM_1)))).';
        r2_filt_B = interp1(Sampling_vector,r2_filt_B,Sampling_vector*(1+((1e-6)*(PPM_2)))).';
        
        % first corr -- corr1_reliability = corr_reliability(corr_signal);
        corr_signal = correlate_iq(r1_filt_B(index_samples(1):index_samples(2)), r2_filt_B(index_samples(1):index_samples(2)), 'dphase', 0);
        [~, idx1] = max(corr_signal);
        delay_time_RS(index1,index2) = idx1 - (index_samples(2)-index_samples(1));
        % time alignment
        if alignment && delay_time_RS(index1,index2) > 0
            r1_alignment = r1_filt_B(delay_time_RS(index1,index2)+1:end);
            r2_alignment = r2_filt_B;
        elseif alignment && delay_time_RS(index1,index2) <= 0
            r1_alignment = r1_filt_B;
            r2_alignment = r2_filt_B(-delay_time_RS(index1,index2)+1:end);
        else
            r1_alignment = r1_filt_B;
            r2_alignment = r2_filt_B;
        end
        % second corr
        corr_signal = correlate_iq(r1_alignment(index_samples(3):index_samples(4)), r2_alignment((index_samples(3):index_samples(4))), 'dphase', 0);
        [~, idx1] = max(corr_signal);
        delay_time_US(index1,index2) = idx1 - (index_samples(4)-index_samples(3));
        % third corr
        corr_signal = correlate_iq(r1_alignment(index_samples(5):min([index_samples(6) numel(r1_alignment) numel(r2_alignment)])), r2_alignment((index_samples(5):min([index_samples(6) numel(r1_alignment) numel(r2_alignment)]))), 'dphase', 0);
        [~, idx1] = max(corr_signal);
        delay_time_RS_ceck(index1,index2) = idx1 - (min([index_samples(6) numel(r1_alignment) numel(r2_alignment)])-index_samples(5));
        
        fprintf('%d %d %d\n',delay_time_RS(index1,index2),delay_time_US(index1,index2),delay_time_RS_ceck(index1,index2))
    end
end
