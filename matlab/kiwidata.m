load('results_mlat/iq1.mat');

dcf77 = [50.01556, 9.01083];
% dcf77 = [52.296667, -2.105278];

kiwicoords = [50, 9] .* ones(24, 2);
kiwicoords = [kiwicoords; 50.1, 8.9; 50.1, 8.9; 50, 8.9; 50, 8.9; 50, 8.9; 50, 8.9];
% kiwicoords = [52.25, -2.20; 52.30, -1.90; 52.25, -2.20; 52.32, -2.12; 52.34, -2.12];

filenum = length(kiwicoords);
errorsk = zeros(filenum,1);
for ii = 1:filenum
   [angle, r] = distance(kiwicoords(ii,1), kiwicoords(ii,2), dcf77(1), dcf77(2));
   errorsk(ii) = angle * 6371 * (pi/180);
end