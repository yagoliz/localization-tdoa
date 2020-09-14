function data = json_load(filename)
% Simple helper function to load json into a struct
fid = fopen(filename); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
data = jsondecode(str);
end