function [doa_meters, doa_samples, doa_meters_2, doa_samples_2, correlation_value, correlation_value_interp, keep] = ... 
    tdoa2_test(signal1_complex, signal2_complex, num_samples_per_freq, num_samples_per_slice, sample_rate, rx_distance_diff, ...
    max_lag, corr_type,  report_level, signal_bandwidth_khz, interpol)

    % tdoa2 calculates the TDOA of two signals captured by two RXs
    %	output:
    %   doa_meters: 		delay in meters (how much signal1 is later than signal 2)
    %	doa_samples: 		delay in samples
    %   reliability: 		reliab. of the correlation, 0(bad)..1(good)
    %   
    %   input:
    %	signal1: 			signal with length 3.6e6 from RX 1
    %	signal2: 			signal with length 3.6e6 from RX 2
    %   rx_distance_diff: 	difference in distance in meters between two RX to Ref (sign matters)
    %   rx_distance:        distance between RX1 and RX2 in meters (always positive)
    %	smooting_factor: 	for wideband signals
    %	corr_type: 			switch between abs and delta phase (abs: 0, delta phase: 1);
    %	report_level:		no reports: 0, show figures >0
    %   signal_bandwidth_khz bandwidth of FIR filter for signal filtering applied to meas signal (400, 200, 40, 12, 0)
    %   interpol            interpolation factor (0 or 1 = no interpolation)
    %
    %   requirement: signal capture with 2 Msps:
    
    %   still open: correct valid signal generation for interpolation 
    %              (currently: native valid signal is used, which gives an approximation, which should be ok)
    
    
    % slice signal into three parts
    % 1111111111111111111111111xxxxxxxxxxxxx2222222222222222222222222xxx3333..
    % |-num_samples_per_slice-|
    % |-num_samples_per_freq+guard_interval-|-num_samples_per_slice-|
    % |-------2*num_samples_per_freq + % guard_interval----------------|..

%     num_samples_per_freq = 4e6;
%     num_samples_per_slice = 1.2e6;
    guard_interval = 2e5; % time to switch to a new frequency, fixed, empirically determined
%     sample_rate = 2e6;  % in Hz

    signal11_complex = signal1_complex(1                                       : num_samples_per_slice*interpol);
    signal12_complex = signal1_complex((num_samples_per_freq   + guard_interval)*interpol : (num_samples_per_freq   + guard_interval + num_samples_per_slice - 1)*interpol);
    signal13_complex = signal1_complex((2*num_samples_per_freq + guard_interval)*interpol : (2*num_samples_per_freq + guard_interval + num_samples_per_slice - 1)*interpol);

    signal21_complex = signal2_complex(1                                       : num_samples_per_slice*interpol);
    signal22_complex = signal2_complex((num_samples_per_freq   + guard_interval)*interpol : (num_samples_per_freq   + guard_interval + num_samples_per_slice - 1)*interpol);
    signal23_complex = signal2_complex((2*num_samples_per_freq + guard_interval)*interpol : (2*num_samples_per_freq + guard_interval + num_samples_per_slice - 1)*interpol);
    
    % This variable holds the output from xcorr
    correlation_value = [0, 0, 0];
    correlation_value_interp = [0, 0, 0];

    %% Filter measurement to signal bandwidth
       
    disp('Filter measurement signal to actual bandwidth');

%     signal12_complex_unfiltered = signal12_complex; % copy unfiltered signals for display
%     signal22_complex_unfiltered = signal22_complex;
%     
%     signal12_complex = filter_iq(signal12_complex, signal_bandwidth_khz);
%     signal22_complex = filter_iq(signal22_complex, signal_bandwidth_khz);
    
    disp(' ');
    
    
    if report_level > 2
        % display spectrum before and after filtering
         
        figure;
        subplot(2,1,1);
        spectrum1 = 10*log10(abs(fftshift(fft(signal12_complex_unfiltered))));
        spectrum1 = smooth(spectrum1, 201);
        
        spectrum2 = 10*log10(abs(fftshift(fft(signal22_complex_unfiltered))));
        spectrum2 = smooth(spectrum2, 201);
        
        freq_axis = -(length(spectrum1)/2) : 1 : ((length(spectrum1)/2)-1);
        plot(freq_axis, spectrum1, freq_axis, spectrum2);
        title('Measurement Signal 1 & 2, unfiltered');
        grid;
        
        subplot(2,1,2);
        spectrum1 = 10*log10(abs(fftshift(fft(signal12_complex))));
        spectrum1 = smooth(spectrum1, 201);
        
        spectrum2 = 10*log10(abs(fftshift(fft(signal22_complex))));
        spectrum2 = smooth(spectrum2, 201);
        
        freq_axis = -(length(spectrum1)/2) : 1 : ((length(spectrum1)/2)-1);
        plot(freq_axis, spectrum1, freq_axis, spectrum2);
        title('Measurement Signal 1 & 2, filtered');
        grid;
    end
    

    
    %% Correlation for slice 1 (ref)
    signal11_complex = filter_iq(signal11_complex, signal_bandwidth_khz);
    signal21_complex = filter_iq(signal21_complex, signal_bandwidth_khz);
    
    % native
    [corr_signal_1, lags1] = correlate_iq(signal11_complex, signal21_complex, corr_type, max_lag);
%     corr1_reliability = corr_reliability( corr_signal_1 );

    [correlation_value(1), idx1] = max(abs(corr_signal_1));
    
    delay1_native = lags1(idx1)/interpol; % >0: signal1 later, <0 signal2 later
        
    %% Correlation for slice 2 (measure)
    [corr_signal_2, lags2] = correlate_iq(signal12_complex, signal22_complex, corr_type, max_lag);
        
    [correlation_value(2), idx2] = max(abs(corr_signal_2));
    delay2_native = lags2(idx2)/interpol; % >0: signal1 later, <0 signal2 later    
    
    %% Correlation for slice 3 (ref check)
    %native

    signal13_complex = filter_iq(signal13_complex, signal_bandwidth_khz);
    signal23_complex = filter_iq(signal23_complex, signal_bandwidth_khz);
    
    
    [corr_signal_3, lags3] = correlate_iq(signal13_complex, signal23_complex, corr_type, max_lag);
%     corr3_reliability = corr_reliability(corr_signal_3);
%     corr_signal_3_valid = zeros(length(corr_signal_3),1) - 1;
%     corr_signal_3_valid(idx1-valid_samples_left-3:idx1+valid_samples_right+3) = corr_signal_3(idx1-valid_samples_left-3:idx1+valid_samples_right+3);
    corr_signal_3_valid = corr_signal_3;

    [correlation_value(3), idx3] = max(abs(corr_signal_3_valid));
    delay3_native = lags3(idx3)/interpol;
        
    %% Display Correlation Signals

    % display reference and reference check signal
    if (report_level > 1)
        figure('units','normalized','outerposition',[0 0 1 1]);
        
        if (interpol > 1) 
            x_corr_interp = 1:(2*interpol*length(signal11_complex) - 1);
            x_corr_interp_plot = (x_corr_interp-1)./interpol +(1/interpol);
        end

        win_len = 10;
        idx = idx1;

        idx_left = idx-win_len;
        if idx_left < 1 
            idx_left = 1;
        end
        idx_right = idx+win_len;
        if idx_right > length(corr_signal_1) 
            idx_right = length(corr_signal_1);
        end

        % First reference signal
        subplot(3,1,1); hold on;
        stem(lags1(idx_left:idx_right), corr_signal_1(idx_left:idx_right));
        
        if (interpol > 1)
            plot(lags1_interp((interpol*idx_left):(interpol*idx_right))./interpol, corr_signal_1_interp((interpol*idx_left):(interpol*idx_right)), 'm.', 'MarkerSize',7);
            plot(delay1_interp, corr_signal_1_interp(idx1_interp_i), 'dg');
        end
        plot(lags1(idx), corr_signal_1(idx), 'dr');
        title('Reference Signal (Slice 1)'); legend('Original','Upsampled'); ylim([-0.1 1.1]); xlim([lags1(idx_left+8) lags1(idx_right)-8]); grid;


        % Measurement signal
        subplot(3,1,2); hold on;
        idx = idx2;

        idx_left = idx-win_len;
        if idx_left < 1 
            idx_left = 1;
        end
        idx_right = idx+win_len;
        if idx_right > length(corr_signal_2) 
            idx_right = length(corr_signal_2);
        end
        
        stem(lags2(idx_left:idx_right) - delay1_native, corr_signal_2(idx_left:idx_right));
        if interpol > 1
            plot(lags2_interp((interpol*idx_left):(interpol*idx_right))./interpol - delay1_native, corr_signal_2_interp((interpol*idx_left):(interpol*idx_right)), 'm.');
            plot(delay2_interp - delay1_native, corr_signal_2_interp(idx2_interp_i), 'dg');     
        end
        plot(lags2(idx) - delay1_native, corr_signal_2(idx), 'dr');
        title('measurement Signal (Slice 2)'); legend('Original','Upsampled'); ylim([-0.1 1.1]); xlim([lags2(idx_left) - delay1_native lags2(idx_right) - delay1_native]); grid;
        
        % Reference signal check
        subplot(3,1,3); hold on;
        idx = idx3;

        idx_left = idx-win_len;
        if idx_left < 1 
            idx_left = 1;
        end
        idx_right = idx+win_len;
        if idx_right > length(corr_signal_3) 
            idx_right = length(corr_signal_3);
        end
        
        stem(lags3(idx_left:idx_right) - delay1_native, corr_signal_3(idx_left:idx_right));
        if interpol > 1
            plot(lags3_interp((interpol*idx_left):(interpol*idx_right))./interpol - delay1_native, corr_signal_3_interp((interpol*idx_left):(interpol*idx_right)), 'm.');
            plot(delay3_interp - delay1_native, corr_signal_3_interp(idx3_interp_i), 'dg');
        end
        plot(lags3(idx) - delay1_native, corr_signal_3(idx), 'dr');
        title('Reference Signal Check (Slice 3)'); legend('Original','Upsampled'); ylim([-0.1 1.1]); xlim([lags3(idx_left) - delay1_native lags3(idx_right) - delay1_native]); grid;
    end

  

    %% Calculate Correlation Results
    
    delay1 = delay1_native;
    delay2 = delay2_native;
    delay3 = delay3_native;
    
    if abs(delay1 - delay3) > 2 * interpol % A small delay can be tolerated (interpol always needs to be 1)
        disp('<strong>WARNING: BAD REFERENCE SIGNALS: ref delays differ by more than 2 samples! </strong>');
%         if corr1_reliability > corr3_reliability
%             avg_delay13 = delay1;
%             disp(['<strong>taking ref with higher reliability, i.e. ref (reliability: ' num2str(corr1_reliability) '(ref) > ' num2str(corr3_reliability) '(ref check) </strong>']);
%         else
%             avg_delay13 = delay3;
%             disp(['<strong>taking ref with higher reliability, i.e. ref check (reliability: ' num2str(corr1_reliability) '(ref) < ' num2str(corr3_reliability) '(ref check) </strong>']);
%         end
        keep = false;
    else
        keep = true;
    end
    avg_delay13 = (delay1 + delay3) / 2; % Assuming delay1 is 0
    
    
    ref_signal_diff_samples = (rx_distance_diff / 3e8) * sample_rate; % known ref signal delay in samples

    % doa_samples/_meters specifies how much signal1 is later than signal2
    doa_samples = delay2 - avg_delay13 + ref_signal_diff_samples; % time difference of arrival without delays due to reception start time and ref transmitter (desc. above)
    doa_samples_2 = delay2 - delay1 + ref_signal_diff_samples;
    doa_meters = (doa_samples / sample_rate) * 3e8;
    doa_meters_2 = (doa_samples_2 / sample_rate) * 3e8;
    

%     reliability = min([corr3_reliability, corr2_reliability, corr1_reliability]);

    if report_level > 0
        disp(' ');
        disp('CORRELATION RESULTS');

        disp(['raw delay1 (ref) (nativ): ' num2str(delay1_native) ', reliability nativ (0..1): ' num2str(correlation_value(1))]);
        disp(['raw delay2 (measure) (nativ): ' num2str(delay2_native) ', reliability nativ: ' num2str(correlation_value(2))]);
        disp(['raw delay3 (ref check) (nativ): ' num2str(delay3_native) ', reliability nativ: ' num2str(correlation_value(3))]);
        disp(['merged delay of ref and ref check: ' num2str(avg_delay13) ]);
        disp(' ');

        disp(['specified distance difference to ref tx [m]: ' int2str(rx_distance_diff)]);
        disp(['specified distance difference to ref tx [samples]: ' num2str(ref_signal_diff_samples)]);
        %disp(['specified distance between two RXes [m]: ' num2str(rx_distance)]);
        disp(' ');

        disp('FINAL RESULT');
        disp(['TDOA in samples: ' num2str(doa_samples) '(how much is signal1 later than signal2)']);
        disp(['TDOA in distance [m]: ' num2str(doa_meters) ]);
        disp(['Total Reliability (min of all 3): ' num2str(min(correlation_value)) ]);
        disp(' ');
    end
end

