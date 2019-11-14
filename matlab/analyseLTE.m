%% Initialization
clear
close all

% Load the configuration
config;

% Get some file specific parameters
NUM_SAMPLES_LTESS = 1.92e6 * 2;
NUM_SAMPLES_TDOA = 1.92e6 * 3;

% Variables that contain the data
data_rtl = zeros(NUM_SAMPLES_LTESS + NUM_SAMPLES_TDOA, RECEIVERS);

% Get the files
% RTL-SDR 0
file_rtl{1} = sprintf('D%d-E%d-ltess.dat', 0, index_to_analyze);
data_rtl(:,1) = spec_load([data_path, file_rtl{1}]);
% RTL-SDR 1
file_rtl{2} = sprintf('D%d-E%d-ltess.dat', 1, index_to_analyze);
data_rtl(:,2) = spec_load([data_path, file_rtl{2}]);
% RTL-SDR 2
file_rtl{3} = sprintf('D%d-E%d-ltess.dat', 2, index_to_analyze);
data_rtl(:,3) = spec_load([data_path, file_rtl{3}]);

%% LO Offset correction
PPM = zeros(1,RECEIVERS);
PPM2 = zeros(1,RECEIVERS);

if lo_correction    
    for receiver = 1:RECEIVERS
        [PPM(receiver), PPM2(receiver)] = ltess(data_rtl(1:NUM_SAMPLES_LTESS, receiver), SAMPLING_RATE_LTESS);
    end
    fprintf('<------------------------------------------->\n')
    for receiver = 1:RECEIVERS
        fprintf('Device %d --> PPM: %f - PPM2: %f\n', receiver, PPM(receiver), PPM2(receiver))
    end
    fprintf('<------------------------------------------->\n')
end

%% TDOA calculation
if tdoa_calculation == 1    
    alignment = true;
    
    Ts = 1/SAMPLING_RATE_TDOA;
    t = 0:1/SAMPLING_RATE_TDOA:((1/SAMPLING_RATE_TDOA)*(NUM_SAMPLES_TDOA-1));
    
    alpha = exp(-Ts/0.5e-6);
    r_filt_b = 1-alpha;
    r_filt_a = [1, -alpha];
    r_filt_z = [];
    
    analyzed_devices = nchoosek(1:RECEIVERS,2);
    
    delay_time_RS = zeros(size(analyzed_devices,1),1);
    delay_time_US = zeros(size(analyzed_devices,1),1);
    
    % First chunk of data (Reference signal)
    index_samples(1) = GUARD_INTERVAL;
    index_samples(2) = index_samples(1) + NUM_SAMPLES_TDOA/3 - GUARD_INTERVAL;
    
    % Second chunk of data (Unknown signal)
    index_samples(3) = index_samples(2) + GUARD_INTERVAL;
    index_samples(4) = index_samples(3) + NUM_SAMPLES_TDOA/3 - GUARD_INTERVAL;
    
    for receiver = 1:RECEIVERS
        for method_index = 1:length(methods)
            
            fprintf('Experiment:%d --- Method: %s\n', receiver, methods{method_index})
            
            r1 = data_rtl(NUM_SAMPLES_LTESS+1:end, analyzed_devices(receiver,1));
            r2 = data_rtl(NUM_SAMPLES_LTESS+1:end, analyzed_devices(receiver,2));

            [r1_filt_B,~] = filter(r_filt_b,r_filt_a,r1,r_filt_z);
            r1_filt_B = r1_filt_B - mean(r1_filt_B);

            [r2_filt_B,r_filt_z_A] = filter(r_filt_b,r_filt_a,r2,r_filt_z);
            r2_filt_B = r2_filt_B - mean(r2_filt_B);
                        
            % FO correction
            % Correction for reference signal
            r1_filt_B(1:NUM_SAMPLES_TDOA/3) = r1_filt_B(1:NUM_SAMPLES_TDOA/3).* ...
                (exp(-1i*2*pi*t(1:NUM_SAMPLES_TDOA/3)*(-PPM(analyzed_devices(receiver,1)))*(1e-6)*fUS).');
            
            r2_filt_B(1:NUM_SAMPLES_TDOA/3) = r2_filt_B(1:NUM_SAMPLES_TDOA/3).* ...
                (exp(-1i*2*pi*t(1:NUM_SAMPLES_TDOA/3)*(-PPM(analyzed_devices(receiver,2)))*(1e-6)*fUS).');
            
            % Correction for unknown signal
            r1_filt_B((1+NUM_SAMPLES_TDOA/3):(2*NUM_SAMPLES_TDOA/3)) = ... 
                r1_filt_B((1+NUM_SAMPLES_TDOA/3):(2*NUM_SAMPLES_TDOA/3)).* ... 
                (exp(-1i*2*pi*t((1+NUM_SAMPLES_TDOA/3):(2*NUM_SAMPLES_TDOA/3))* ...
                (-PPM(analyzed_devices(receiver,1)))*(1e-6)*fRS).');
            
            r2_filt_B((1+NUM_SAMPLES_TDOA/3):(2*NUM_SAMPLES_TDOA/3)) = ...
                r2_filt_B((1+NUM_SAMPLES_TDOA/3):(2*NUM_SAMPLES_TDOA/3)).* ...
                (exp(-1i*2*pi*t((1+NUM_SAMPLES_TDOA/3):(2*NUM_SAMPLES_TDOA/3))* ...
                (-PPM(analyzed_devices(receiver,2)))*(1e-6)*fRS).');
            
            % Time correction
            Sampling_vector = 1:numel(r1_filt_B);
            r1_filt_B = interp1(Sampling_vector,r1_filt_B,Sampling_vector*(1+((1e-6)*(PPM(analyzed_devices(receiver,1)))))).';
            r2_filt_B = interp1(Sampling_vector,r2_filt_B,Sampling_vector*(1+((1e-6)*(PPM(analyzed_devices(receiver,2)))))).';
                        
            % First correlation - To obtain the delay 
            [corr_signal, idx1] = correlate_iq(r1_filt_B(index_samples(1):index_samples(2)), r2_filt_B(index_samples(1):index_samples(2)), methods{method_index}, 0);
            delay_time_RS(receiver,method_index) = idx1;
            
            % Time alignment
            if alignment && delay_time_RS(receiver,method_index) > 0
                r1_alignment = r1_filt_B(delay_time_RS(receiver,method_index)+1:end);
                r2_alignment = r2_filt_B;
            elseif alignment && delay_time_RS(receiver,method_index) <= 0
                r1_alignment = r1_filt_B;
                r2_alignment = r2_filt_B(-delay_time_RS(receiver,method_index)+1:end);
            else
                r1_alignment = r1_filt_B;
                r2_alignment = r2_filt_B;
            end
            
            % Second correlation - Unknown signal
            chunk_indices = linspace(index_samples(3), index_samples(4), chunks + 1);
            idx2 = zeros(1,length(chunk_indices)-1);
            
            for i = 1:length(chunk_indices) - 1
                signal1 = r1_alignment(chunk_indices(i):chunk_indices(i+1));
                signal1_res = resample(signal1, resampling_factor, 1);
                
                signal2 = r2_alignment(chunk_indices(i):chunk_indices(i+1));
                signal2_res = resample(signal2, resampling_factor, 1);
                
                [corr_signal, idx2(i)] = correlate_iq(signal1_res, signal2_res, methods{method_index}, 0);
                idx2(i) = idx2(i) / resampling_factor;
            end
            
            delay_time_US(receiver,method_index) = median(idx2);
            
            % Print results
            fprintf('---\nResults:\n');
            fprintf('Delay on first chunk: %d | TDOA: %d\n',delay_time_RS(receiver,method_index),delay_time_US(receiver,method_index));
            fprintf('---------------------------------------------------\n');
        end
    end
    
end