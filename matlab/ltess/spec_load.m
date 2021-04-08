function [ Y ] = spec_load( filename, data_type )

    if nargin == 1
        data_type = 'unsigned char';
    end

    fd = fopen(filename,'r');
    
    if strcmp(data_type, 'unsigned char')

        %data = fread(fd, 9600*200*2, 'unsigned char');
        data = fread(fd, 'unsigned char');
        real = data(1:2:size(data,1));
        imag = data(2:2:size(data,1));

        min_size = min(size(real,1), size(imag,1));
        real = real(1:min_size);
        imag = imag(1:min_size);

        real = (real-127.0)/128.0;
        imag = (imag-127.0)/128.0;

        complex_signal = complex(real,imag);
    elseif strcmp(data_type, 'float')
        data = fread(fd, 'float');
        real = data(1:2:size(data,1));
        imag = data(2:2:size(data,1));

        min_size = min(size(real,1), size(imag,1));
        real = real(1:min_size);
        imag = imag(1:min_size);
        complex_signal = complex(real,imag);
        
    end

    fclose(fd);

    Y = complex_signal;

end

