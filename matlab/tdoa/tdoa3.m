function [doa_meters, doa_samples, doa_meters_2, doa_samples_2, correlation_value, correlation_value_interp, mp_delay] = ...
    tdoa3(signal1_complex, signal2_complex, num_samples_per_freq, num_samples_per_slice, sample_rate, rx_distance_diff, ...
    max_lag, corr_type,  report_level, signal_bandwidth_khz_rs, signal_bandwidth_khz_us, interpol, correct_multipath, mp_si)

global c;

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

%   Constants for our tdoa estimation
guard_interval = 5e4; % time to switch to a new frequency, fixed, empirically determined
peak_find = round(max_lag);

%     Signal preparation
signal11_complex = signal1_complex(guard_interval                          : guard_interval + num_samples_per_slice - 1);
signal13_complex = signal1_complex(2*num_samples_per_freq + guard_interval : 2*num_samples_per_freq + guard_interval + num_samples_per_slice - 1);

signal21_complex = signal2_complex(guard_interval                          : guard_interval + num_samples_per_slice - 1);
signal23_complex = signal2_complex(2*num_samples_per_freq + guard_interval : 2*num_samples_per_freq + guard_interval + num_samples_per_slice - 1);

% This variable holds the output from xcorr
correlation_value = [0, 0, 0];
correlation_value_interp = [0, 0, 0];

%% Correlation for slice 1 (ref)
signal11_complex = filter_iq(signal11_complex, signal_bandwidth_khz_rs);
signal21_complex = filter_iq(signal21_complex, signal_bandwidth_khz_rs);

% native
[corr_signal_1, lags1] = correlate_iq(signal11_complex, signal21_complex, corr_type);
[correlation_value(1), idx1] = max(corr_signal_1);
delay1_native = lags1(idx1); % >0: signal1 later, <0 signal2 later

% with interpolation
if (interpol > 1)
    lags1_interp = lags1(idx1-peak_find):1/interpol:lags1(idx1+peak_find);
    corr_signal_1_interp = interp1(lags1(idx1-peak_find:idx1+peak_find),...
        corr_signal_1(idx1-peak_find:idx1+peak_find),...
        lags1_interp,...
        "makima");
    %        [correlation_value_interp(1),idx1_interp] = max(abs(corr_signal_1_interp));
    [correlation_value_interp(1),idx1_interp] = max(abs(corr_signal_1_interp));
    delay1_interp = lags1_interp(idx1_interp);
else
    delay1_interp = 0;
end

%% Correlation for slice 3 (ref check)
signal13_complex = filter_iq(signal13_complex, signal_bandwidth_khz_rs);
signal23_complex = filter_iq(signal23_complex, signal_bandwidth_khz_rs);


[corr_signal_3, lags3] = correlate_iq(signal13_complex, signal23_complex, corr_type);
corr_signal_3_valid = zeros(length(corr_signal_1),1);
corr_signal_3_valid(idx1-round(max_lag)/2:idx1+round(max_lag)/2) = corr_signal_3(idx1-round(max_lag)/2:idx1+round(max_lag)/2);

[correlation_value(3), idx3] = max(corr_signal_3_valid);
delay3_native = lags3(idx3);

% with interpolation
if (interpol > 1)
    lags3_interp = lags3(idx3-peak_find):1/interpol:lags3(idx3+peak_find);
    corr_signal_3_interp = interp1(lags3(idx3-peak_find:idx3+peak_find),...
        corr_signal_3(idx3-peak_find:idx3+peak_find),...
        lags3_interp,...
        "makima");
    [correlation_value_interp(3),idx3_interp] = max(corr_signal_3_interp);
    delay3_interp = lags3_interp(idx3_interp);
else
    delay3_interp = 0;
end

%% Correlation for slice 2 (measure)
if interpol > 1
    signal12_complex = signal1_complex(num_samples_per_freq + delay1_interp + guard_interval : num_samples_per_freq + guard_interval + delay1_interp + num_samples_per_slice - 1);
else
    signal12_complex = signal1_complex(num_samples_per_freq + delay1_native + guard_interval : num_samples_per_freq + guard_interval + delay1_native + num_samples_per_slice - 1);
end
signal22_complex = signal2_complex(num_samples_per_freq + guard_interval : num_samples_per_freq + guard_interval   + num_samples_per_slice - 1);

signal12_complex = filter_iq(signal12_complex, signal_bandwidth_khz_us);
signal22_complex = filter_iq(signal22_complex, signal_bandwidth_khz_us);

[corr_signal_2, lags2] = correlate_iq(signal12_complex, signal22_complex, corr_type);

% Multipath correction
mp_delay = 0;
if correct_multipath
    [~, idx2] = max(corr_signal_2);
    lags_normalized = lags2(idx2-max_lag:idx2+max_lag);
    corr_normalized = corr_signal_2(idx2-max_lag:idx2+max_lag);
    corr_normalized = sign(cos(angle(corr_normalized))).*abs(corr_normalized); 
    [peaks, locations] = findpeaks(corr_normalized/max(abs(corr_signal_2)),'MinPeakHeight',0.3,'MinPeakProminence',0.01);

    % If we find multipath we have to correct it
    if length(peaks) > 1
        max_peak = max(locations) - max(mp_si); % Peak is at maximum multipath component and then to the left
        distances = abs(locations - max_peak); % We select the closest peak (hopefully it lands on one peak)

        [min_distance,real_loc] = min(distances);
        if min_distance > 0.5
            fprintf('WARNING: Peak appeared at a larger distance than expected: %f',min_distance);
        end
        real_peak = locations(real_loc);

        correlation_value(2) = corr_normalized(real_peak);
        delay2_native = lags_normalized(real_peak);
        mp_delay = abs(min(locations) - real_peak); 
    else
        [correlation_value(2), idx2] = max(corr_signal_2);
        delay2_native = lags2(idx2); % >0: signal1 later, <0 signal2 later
    end
else
    [correlation_value(2), idx2] = max(corr_signal_2);
    delay2_native = lags2(idx2); % >0: signal1 later, <0 signal2 later
end

% with interpolation, only one valid
if (interpol > 1)
    lags2_interp = lags2(idx2-peak_find):1/interpol:lags2(idx2+peak_find);
    corr_signal_2_interp = interp1(lags2(idx2-peak_find:idx2+peak_find),...
        corr_signal_2(idx2-peak_find:idx2+peak_find),...
        lags2_interp,...
        "makima");
    
    if correct_multipath
        index = find(lags2_interp == delay2_native);
        [correlation_value_interp(2), ~] = max(corr_signal_2_interp(index-interpol-1:index+interpol+1));
        delay2_interp = lags2_interp(corr_signal_2_interp == correlation_value_interp(2));
        
    else
        [correlation_value_interp(2),idx2_interp] = max(corr_signal_2_interp);
        delay2_interp = lags2_interp(idx2_interp);
    end
else
    delay2_interp = 0;
end


%% Calculate Correlation Results
if (interpol <= 1)
    delay2 = delay2_native;
else
    delay2 = delay2_interp;
end

ref_signal_diff_samples = (rx_distance_diff / c) * sample_rate; % known ref signal delay in samples

% doa_samples/_meters specifies how much signal1 is later than signal2
doa_samples = delay2 + ref_signal_diff_samples; % time difference of arrival without delays due to reception start time and ref transmitter (desc. above)
doa_samples_2 = doa_samples;
doa_meters = (doa_samples / sample_rate) * c;
doa_meters_2 = (doa_samples_2 / sample_rate) * c;

if report_level > 0
    disp(' ');
    disp('CORRELATION RESULTS');

    disp(['raw delay1 (ref) (nativ/interp): ' int2str(delay1_native) ' / ' num2str(delay1_interp) ', reliability nativ (0..1): ' num2str(correlation_value(1))]);
    disp(['raw delay2 (measure) (nativ/interp): ' int2str(delay2_native) ' / ' num2str(delay2_interp) ', reliability nativ: ' num2str(correlation_value(2))]);
    disp(['raw delay3 (ref check) (nativ/interp): ' int2str(delay3_native) ' / ' num2str(delay3_interp) ', reliability nativ: ' num2str(correlation_value(3))]);
    disp(' ');

    disp(['specified distance difference to ref tx [m]: ' int2str(rx_distance_diff)]);
    disp(['specified distance difference to ref tx [samples]: ' num2str(ref_signal_diff_samples)]);
    disp(' ');

    disp('FINAL RESULT');
    disp(['TDOA in samples: ' num2str(doa_samples) '(how much is signal1 later than signal2)']);
    disp(['TDOA in distance [m]: ' num2str(doa_meters) ]);
    disp(['Total Reliability (min of all 3): ' num2str(min(correlation_value)) ]);
    disp(' ');
end
end

