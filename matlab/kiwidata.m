load('results_mlat/iq1.mat');

dcf77 = [50.01556, 9.01083];

kiwicoords = [50, 9] .* ones(24, 2);
kiwicoords = [kiwicoords; 50.1, 8.9; 50.1, 8.9; 50, 8.9; 50, 8.9; 50, 8.9; 50, 8.9];

errorsk = zeros(30,1);
for ii = 1:30
   [angle, r] = distance(kiwicoords(ii,1), kiwicoords(ii,2), dcf77(1), dcf77(2));
   errorsk(ii) = angle * 6371 * (pi/180);
end