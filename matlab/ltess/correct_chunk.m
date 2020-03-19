function signal_corrected = correct_chunk(signal, PPM, samplingRate, f)
%CORRECT_FO Correct local oscillator offset of RTL-SDR
% -signal: vector with complex IQ
% -PPM: LO offset in parts per million
% -fRS: Frequency of the reference signal (Hz)
% -fUS: Frequency of the unknown signal (Hz)
% Signal is assumed to be split in 3 chunks of equal length

Ts = 1/samplingRate;
phi = PPM * 1e-6;
t = 0:Ts/(1+phi):Ts/(1+phi)*(length(signal) - 1);

% Chunk correction
signal = signal .* (exp(-1i * 2 * pi .* t * PPM * (1e-6) * f)');

% Time correction
signal_corrected = resample(signal, t, samplingRate);

% Eliminate al odd numbers (Nan)
signal_corrected(isnan(signal_corrected)) = 0;

end

