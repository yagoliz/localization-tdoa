function signal_corrected = correct_fo(signal, PPM, samplingRate, fRS, fUS)
%CORRECT_FO Correct local oscillator offset of RTL-SDR
% -signal: vector with complex IQ
% -PPM: LO offset in parts per million
% -fRS: Frequency of the reference signal (Hz)
% -fUS: Frequency of the unknown signal (Hz)
% Signal is assumed to be split in 3 chunks of equal length

Ts = 1/samplingRate;
phi = PPM * 1e-6;
t = 0:Ts/(1+phi):Ts/(1+phi)*(length(signal) - 1);
% t = 0:Ts:Ts*(length(signal) - 1);
chunk_length = length(signal)/3;

% First chunk correction
reference_chunk = signal(1 : chunk_length);
reference_time = t(1 : chunk_length);
reference_chunk = reference_chunk .* (exp(-1i * 2 * pi .* reference_time * PPM * (1e-6) * fRS)');

% Second chunk correction
unknown_chunk = signal(chunk_length + 1 : 2 * chunk_length);
unknown_time = t(chunk_length + 1 : 2 * chunk_length);
unknown_chunk = unknown_chunk .* (exp(-1i * 2 * pi .* unknown_time * PPM * (1e-6) * fUS)');

% Third chunk correction
secondary_reference_chunk = signal(2*chunk_length + 1 : end);
secondary_reference_time = t(2 * chunk_length + 1 : end);
secondary_reference_chunk = secondary_reference_chunk .* (exp(-1i * 2 * pi .* secondary_reference_time * PPM * (1e-6) * fRS)');

% Time correction
signal_corrected = [reference_chunk; unknown_chunk; secondary_reference_chunk];
% sampling_vector = 1:numel(signal_corrected);
% signal_corrected = interp1(sampling_vector, signal_corrected, sampling_vector.*(1+(1e-6) * PPM))';
signal_corrected = resample(signal_corrected, t, samplingRate);

% Eliminate al odd numbers (Nan)
signal_corrected(isnan(signal_corrected)) = 0;

end

