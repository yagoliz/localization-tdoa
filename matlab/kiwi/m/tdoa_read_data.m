%% Read data from tdoa files

function [input,status]=tdoa_read_data(plot_info, input, directory)
  if nargin == 1
    directory = 'gnss_pos';
    printf('using default dir="gnss_pos"\n');
  end

  n = numel(input);
  for ii=1:n
    tic;
    [input(ii).name, ...
     input(ii).vname, ...
     input(ii).fname, ...
     input(ii).time, ...
     input(ii).freq]  = parse_iq_filename(input(ii).fn);
    input(ii).coord   = get_coord(input(ii).vname, directory);
    [x,xx,fs,gpsfix] = proc_kiwi_iq_wav(input(ii).fn, 255);
    input(ii).gpsfix  = gpsfix;
    input(ii).use     = false;
    tmin(ii)          = NaN;
    tmax(ii)          = NaN;
    status.per_file(ii).name      = input(ii).name;
    status.per_file(ii).file_name = input(ii).fn;
    status.per_file(ii).time_sec  = toc;
    status.per_file(ii).message   = '';
    status.per_file(ii).last_gnss_fix = gpsfix;
    if gpsfix == 255
      printf('tdoa_read_data: %-40s no GPS timestamps\n', input(ii).fn);
      status.per_file(ii).message = 'no GNSS timestamps';
      continue
    end

    if gpsfix == 254
      printf('tdoa_read_data: %-40s no recent GPS timestamps\n', input(ii).fn);
      status.per_file(ii).message = 'no recent GNSS timestamps';
      continue
    end

    input(ii).t      = cat(1,xx.t);
    input(ii).z      = cat(1,xx.z);
    input(ii).gpssec = cat(1,x.gpssec)+1e-9*cat(1,x.gpsnsec);
    if numel(input(ii).t) == 0 || numel(input(ii).gpssec) <= 2
      printf('tdoa_read_data: %-40s number of samples = %d == 0 || number of blocks = %d <= 2\n', ...
             input(ii).fn, numel(input(ii).t), numel(input(ii).gpssec));
      status.per_file(ii).message = sprintf('number of samples = %d == 0 || number of blocks = %d <= 2', ...
                                           numel(input(ii).t), numel(input(ii).gpssec));
      continue
    end
    if max(input(ii).z) == 0
      printf('tdoa_read_data: %-40s max(z)==0\n', input(ii).fn);
      status.per_file(ii).message = 'max(z)==0';
      continue;
    end
    if max(abs(diff(input(ii).t))) > 2/fs
      printf('tdoa_read_data: %-40s max(abs(diff(input(i).t))) = %f > %f\n', ...
             input(ii).fn, max(abs(diff(input(ii).t))), 2/fs);
      status.per_file(ii).message = sprintf('max(abs(diff(input(i).t))) = %f > %f', ...
                                           max(abs(diff(input(ii).t))), 2/fs);
      continue
    end
    tmin(ii)      = min(input(ii).t);
    tmax(ii)      = max(input(ii).t);
    
%     input(ii).fs  = 512/mean(diff(input(ii).gpssec)(2:end));
    input(ii).use = true;
    printf('tdoa_read_data: %-40s %s last_gnss_fix=%3d\n', ...
           input(ii).fn, input(ii).name, gpsfix);
  end

  % Exclude bad stations
  b_use = vertcat(input.use);
  tmin  = tmin(b_use);
  tmax  = tmax(b_use);
  n     = numel(input);

  % Truncate all data to the same common time interval
  t0 = max(tmin);
  t1 = min(tmax);
  for ii=1:n
    if ~input(ii).use
      continue
    end
    b = input(ii).t<t0 | input(ii).t>t1;
    input(ii).t(b) = [];
    input(ii).z(b) = [];
    input(ii).use  = numel(input(ii).z)/input(ii).fs > 10;
    if ~input(ii).use
      printf('tdoa_read_data: %-40s excluded (%.2f sec < %g sec overlap)\n', ...
             input(ii).fn, numel(input(ii).z)/input(ii).fs, 10);
      status.per_file(ii).message = sprintf('excluded (%.2f sec < %g sec overlap)', numel(input(ii).z)/12000, 10);
    end
  end

  % Generate json info for each station
  counter=1;
  stn_status = {'BAD', 'GOOD'};
  for ii=1:n
    if ~input(ii).use
      status.per_file(ii).idx    = -1;
      status.per_file(ii).status = 'BAD';
    else
      status.per_file(ii).idx    = counter;
      status.per_file(ii).status = 'GOOD';
      counter = counter + 1;
    end
  end

  % Exlude bad stations
  input = input(vertcat(input.use));

  if numel(input) < 2
    printf('tdoa_read_data: n=%d < 2 good stations found\n', numel(input));
    status.result = struct('status', 'BAD',...
                           'message', sprintf('%d/%d good stations < 2', numel(input), n));
  else
    status.result = struct('status', 'GOOD',...
                           'message', sprintf('%d/%d good stations', numel(input), n));
  end

  % Resample if necessary (20.25 kHz vs. 12 kHz modes)
  [input,status] = resample_ifneeded(input,status);
end

%% Resampling
function [input,status]=resample_ifneeded(input, status)
  % Round sampling frequencies to nearest multiple of 10 Hz
  fs        = round(cat(1,input.fs)/10)*10;
  [fs0,idx] = min(fs);
  for ii=1:numel(input)
    if ii==idx
      continue
    end
    if abs(fs(ii)/fs0-1) > 0.1
      status.per_file(ii).message = sprintf('resampled %g kHz to %g kHz', 1e-3*[fs(ii) fs0]);
      input(ii).z   = resample(input(ii).z, fs0, fs(ii)); % Factor fs0/fs(i)
      dt = mean(diff(input(ii).t)) *fs(ii)/fs0;
      input(ii).t   = input(ii).t(1) + (0:numel(input(ii).z))*dt;
      input(ii).fs  = 1/dt;
    end
  end
end

%% Parsing IQ files
% e.g. fn = '../files/02697/20180707T211018Z_77500_F1JEK-P_iq.wav'
function [name, vname, fname, time, freq] = parse_iq_filename(fn)
  [~, filename, ext] = fileparts(fn);
  if ~strcmp(ext, '.wav')
    error('Wrong extension: %s', fn);
  end
  tokens = strsplit(filename, '_');
  if numel(tokens) < 4
    error('Malformed filename: %s', fn);
  end
  if ~strcmp(tokens{4}, 'iq')
    error('Filename does not indicate an IQ recording: %s', fn);
  end
  time  = tokens{1};
  freq  = 1e-3 * str2double(tokens{2});
  fname = tokens{3};
  name  = strrep(fname, '-', '/'); % Recover encoded slashes
  vname = strrep(fname, '-', '_');
end
