%% Process data to localize dcf77

% Filenames
files = { ...
'20200420T084345Z_77500_DL0ABT_iq.wav', ...
'20200420T084345Z_77500_S57BIT_iq.wav', ...
'20200420T084346Z_77500_F1JEK_iq.wav', ...
};

% Main directory
directory = [pwd, filesep, 'gnss_pos'];

% Configuration (empty struct)
config = struct();

% Process data
[tdoa, input] = proc_tdoa_kiwi(directory, files, config);