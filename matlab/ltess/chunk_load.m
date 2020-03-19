function [ Y ] = chunk_load( fd, num_samples )

data = fread(fd, num_samples*2, 'unsigned char');
real = data(1:2:size(data,1));
imag = data(2:2:size(data,1));

min_size = min(size(real,1), size(imag,1));
real = real(1:min_size);
imag = imag(1:min_size);

real = (real-127.0)/128.0;
imag = (imag-127.0)/128.0;

complex_signal = complex(real,imag);

Y = complex_signal;

end

