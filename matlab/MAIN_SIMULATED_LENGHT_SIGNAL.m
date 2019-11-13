%% STEP 0 -- INIT

clear
close all

step_1 = 1;
step_2 = 1;

N_receiver = 3;
Test_for_dfo = 3;
Test_for_y1y2 = 25;
fRS = 806e6;
fUS = 244e6;

SAMPLING_RATE_DFO = 1.92e6;
% num_sample_dfo minimun 1s of samples
fs_y1y2 = 2.048e6;
num_samples_y1y2 = 3e6;
                

%% STEP 1 -- LTESS-TRACK

if step_1
    addpath('./common/')
    
    PPM_final = zeros(Test_for_dfo,N_receiver);
    PPM2_final = zeros(Test_for_dfo,N_receiver);
    
    RESAMPLE_FACTOR = 60;
    PSS_STEP = 9600;
    SEARCH_WINDOW = 150;
    CORRELATION_FACTOR = 0.1;
    CORRELATION_FACTOR_PREAMBLE = 0.1;
    PREAMBLE=20; % Number of analyzed PSS before start to jump
    FLIP=0; % 0 disabled 1: enabled  (testing purposes)
    POLYNOMIAL_DEGREE=1;
    
    [Z, Z_t] = get_Zadoof();
    
    for index1 = 1:N_receiver
        for index2 = 1:Test_for_dfo
            
            filename_for_LTESS = sprintf('LTESS_test1\\test_%d_device_%d.dat',index2,index1);
            
            capbuf = spec_load(filename_for_LTESS);
            chunk = capbuf(1:SAMPLING_RATE_DFO*1);
            [PPM,PSS_percent,data_Z,Y_Z,p_Z,PPM2] = getDrift(chunk,SAMPLING_RATE_DFO, Z, Z_t, PREAMBLE, ...
                PSS_STEP, SEARCH_WINDOW, CORRELATION_FACTOR, RESAMPLE_FACTOR, FLIP, POLYNOMIAL_DEGREE);
            fprintf('Zadoof adaptation --> PPM: %f [%f] - PSS detected: %f\n', PPM, PPM2, PSS_percent)
            T_s = 1/SAMPLING_RATE_DFO;
            delta_f=(PPM*1e-6)*806e6;
            Z_t_rotated = {};
            Z_t_rotated{1} = Z_t{1}.*exp(-1i*2*pi*T_s*delta_f*(1:length(Z_t{1}))');
            Z_t_rotated{2} = Z_t{2}.*exp(-1i*2*pi*T_s*delta_f*(1:length(Z_t{2}))');
            Z_t_rotated{3} = Z_t{3}.*exp(-1i*2*pi*T_s*delta_f*(1:length(Z_t{3}))');
            [PPM,PSS_percent,data,Y,p,PPM2, th_learned, Z_sequence, last, p_loc20 ] = getDrift(chunk,SAMPLING_RATE_DFO, Z, Z_t_rotated, ...
                PREAMBLE, PSS_STEP, SEARCH_WINDOW, CORRELATION_FACTOR, RESAMPLE_FACTOR, FLIP, POLYNOMIAL_DEGREE);
            fprintf('Final result --> PPM: %f [%f] - PSS detected: %f\n', PPM, PPM2, PSS_percent)
            PPM_final(index2,index1) = PPM;
            PPM2_final(index2,index1) = PPM2;
        end
    end
    fprintf('<------------------------------------------->\n')
    for index1 = 1:N_receiver
        for index2 = 1:Test_for_dfo
            fprintf('Device %d --> PPM: %f - PPM2: %f\n', index1, PPM_final(index2,index1), PPM2_final(index2,index1))
        end
    end
    fprintf('<------------------------------------------->\n')
    
end

%% STEP 2 -- SIGNAL ANALYSIS

index_samples_T = [0.2e6 0.8e6 1.2e6 1.8e6 2.2e6 2.8e6; ...
    0.2e6 0.8e6 1.3e6 1.7e6 2.3e6 2.7e6; ...
    0.2e6 0.8e6 1.4e6 1.6e6 2.4e6 2.6e6; ...
    0.2e6 0.8e6 1.45e6 1.55e6 2.45e6 2.55e6; ...
    0.2e6 0.8e6 1.49e6 1.51e6 2.49e6 2.51e6; ...
    0.2e6 0.8e6 1.499e6 1.501e6 2.499e6 2.501e6; ...
    0.2e6 0.8e6 1.4999e6 1.5001e6 2.4999e6 2.5001e6; ...
    0.2e6 0.8e6 1.49995e6 1.50005e6 2.49995e6 2.50005e6; ...
    0.2e6 0.8e6 1.499975e6 1.500025e6 2.499975e6 2.500025e6; ...
    0.2e6 0.8e6 1.499990e6 1.500010e6 2.499990e6 2.500010e6];

for varius_signal_lenght = 1:10
    if step_2
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
        
        if Test_for_dfo > 1
            PPM_final_def = mean(PPM_final);
        else
            PPM_final_def = PPM_final;
        end
        
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
                
                index_samples = index_samples_T(varius_signal_lenght,:);
                
                % FO correction
                r1_filt_B(1:num_samples_y1y2/3) = r1_filt_B(1:num_samples_y1y2/3).*(exp(-1i*2*pi*t(1:num_samples_y1y2/3)*(-PPM_final_def(analyzed_devices(index1,1)))*(1e-6)*fUS).');
                r2_filt_B(1:num_samples_y1y2/3) = r2_filt_B(1:num_samples_y1y2/3).*(exp(-1i*2*pi*t(1:num_samples_y1y2/3)*(-PPM_final_def(analyzed_devices(index1,2)))*(1e-6)*fUS).');
                r1_filt_B((1+num_samples_y1y2/3):(2*num_samples_y1y2/3)) = r1_filt_B((1+num_samples_y1y2/3):(2*num_samples_y1y2/3)).*(exp(-1i*2*pi*t((1+num_samples_y1y2/3):(2*num_samples_y1y2/3))*(-PPM_final_def(analyzed_devices(index1,1)))*(1e-6)*fRS).');
                r2_filt_B((1+num_samples_y1y2/3):(2*num_samples_y1y2/3)) = r2_filt_B((1+num_samples_y1y2/3):(2*num_samples_y1y2/3)).*(exp(-1i*2*pi*t((1+num_samples_y1y2/3):(2*num_samples_y1y2/3))*(-PPM_final_def(analyzed_devices(index1,2)))*(1e-6)*fRS).');
                r1_filt_B((1+2*num_samples_y1y2/3):end) = r1_filt_B((1+2*num_samples_y1y2/3):end).*(exp(-1i*2*pi*t((1+2*num_samples_y1y2/3):end)*(-PPM_final_def(analyzed_devices(index1,1)))*(1e-6)*fUS).');
                r2_filt_B((1+2*num_samples_y1y2/3):end) = r2_filt_B((1+2*num_samples_y1y2/3):end).*(exp(-1i*2*pi*t((1+2*num_samples_y1y2/3):end)*(-PPM_final_def(analyzed_devices(index1,2)))*(1e-6)*fUS).');
                % time correction
% % % % % % %                 r1_filt_B = interp(r1_filt_B,4);
% % % % % % %                 r2_filt_B = interp(r2_filt_B,4);
                Sampling_vector = 1:numel(r1_filt_B);
                r1_filt_B = interp1(Sampling_vector,r1_filt_B,Sampling_vector*(1+((1e-6)*(PPM_final_def(analyzed_devices(index1,1)))))).';
                r2_filt_B = interp1(Sampling_vector,r2_filt_B,Sampling_vector*(1+((1e-6)*(PPM_final_def(analyzed_devices(index1,2)))))).';
                
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
        
    end
    
    delay_time_US_GLOBAL(varius_signal_lenght,:) = sum(abs(delay_time_US'));
    delay_time_RS_ceck_GLOBAL(varius_signal_lenght,:) = sum(abs(delay_time_RS_ceck'));
    
end

X_scale = [0.6e6 0.4e6 0.2e6 0.1e6 0.02e6 2e3 200 100 50 20];

subplot(2,1,1)
semilogx(X_scale,delay_time_US_GLOBAL/25,'linewidth',3);
title("Corr - Unknown signal")
legend('RTL1 - RTL2','RTL1 - RTL3','RTL2 - RTL3','Location','NorthEast')
xlabel('# Samples','FontSize', 18)
ylabel('Average TDOA [#S]','FontSize', 18)
grid; 
axis([10 630000 0 1.5]);
set(gca,'fontsize',18)

subplot(2,1,2)
semilogx(X_scale,delay_time_RS_ceck_GLOBAL/25,'linewidth',3);
title("Corr - Reference signal check")
legend('RTL1 - RTL2','RTL1 - RTL3','RTL2 - RTL3','Location','NorthEast') 
xlabel('# Samples','FontSize', 18)
ylabel('Average TDOA [#S]','FontSize', 18)
grid; 
axis([10 630000 0 1.5]);
set(gca,'fontsize',18)

sgtitle('50 tests - With corr synchronization','fontweight','bold','FontSize', 18)

