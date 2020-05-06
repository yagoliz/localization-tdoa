%% Kiwi Information

function sensors = kiwi_info(directory)
  files = dir([directory filesep() '*.m']);
  for ii=1:length(files)
    try
      run(files(ii).name);
    catch err
      error('%s: %s', files(ii), err.message);
    end
  end
end
