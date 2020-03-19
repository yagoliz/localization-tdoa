function [ doa_meters, delay ] = tdoa(signal1_complex, signal2_complex, sample_rate, num_samples_per_slice, smoothing_factor, corr_type)

signal1_complex = signal1_complex(1:num_samples_per_slice);
signal2_complex = signal2_complex(1:num_samples_per_slice);

%% Correlation 
[corr_signal, lags] = correlate_iq(signal1_complex, signal2_complex, corr_type, smoothing_factor);
[~, idx1] = max(corr_signal);

delay = lags(idx1); % >0: signal1 later, <0 signal2 later
doa_meters = (delay / sample_rate) * 3e8;

end

