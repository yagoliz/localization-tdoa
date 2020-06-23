%% Analyze DCF77 and GYN2BR

clear all; close all;
% DCF77
load('data/dcf77_2_pre.mat');
input = alignInputs(input);

iq1 = input(1).z;
t1 = input(1).t - input(1).t(1);
phase1 = unwrap(angle(iq1));

iq2 = input(2).z;
t2 = input(2).t - input(2).t(1);

iq3 = input(3).z;
t3 = input(3).t - input(3).t(1);
phase3 = unwrap(angle(iq3));

% FFTs
freqs1 = linspace(-6, 6, length(iq1));
fft1 = fftshift(fft(iq1));
fft1 = fft1(freqs1 >= -0.5 & freqs1 <= 0.5);
iq1mod = ifft(ifftshift(fft1));
phase1mod = unwrap(angle(iq1mod));
t1mod = linspace(t1(1), t1(end), length(iq1mod));

fft3 = fftshift(fft(iq3));
fft3 = fft3(freqs1 >= -0.5 & freqs1 <= 0.5);
iq3mod = ifft(ifftshift(fft3));
phase3mod = unwrap(angle(iq3mod));
t3mod = linspace(t3(1), t3(end), length(iq3mod));


% Plotting
figure(1); grid on; hold on; plot(t1, phase1); plot(t3, phase3); legend('Sensor 1', 'Sensor 3'); 
xlabel('Time (s)'); ylabel('Phase (º)'); title('Phase using all data');

figure(2); grid on; hold on; plot(t1mod, phase1mod); plot(t3mod, phase3mod); legend('Sensor 1', 'Sensor 3');
xlabel('Time (s)'); ylabel('Phase (º)'); title('Phase when cropping');

%% Alignment function
function inputAli = alignInputs(input)

  t0 = max([input(1).t(1), input(2).t(1), input(3).t(1)]);
  t1 = min([input(1).t(end), input(2).t(end), input(3).t(end)]);
  
  n = numel(input);
  for i=1:n
    if ~input(i).use
      continue
    end
    b = input(i).t<t0 | input(i).t>t1;
    input(i).t(b) = [];
    input(i).z(b) = [];
    input(i).use  = numel(input(i).z)/input(i).fs > 10;
    if ~input(i).use
      printf('tdoa_read_data: %-40s excluded (%.2f sec < %g sec overlap)\n', ...
             input(i).fn, numel(input(i).z)/input(i).fs, 10);
    end
  end
  
  inputAli = input;
  
end


