%% Get coordinates

function latlon=get_coord(vname, directory)
  sensors = kiwi_info(directory);
  latlon=getfield(sensors, vname).coord;
end

