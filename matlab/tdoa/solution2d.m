function [x, y] = solution2d(doa_array,sensors)
%SOLUTION2D Compute the TDOA algorithm and return the optimal position
%   Input:
%   - doa_array = array with the TDOAs between different pairs of sensors.
%   It is assumed that the first sensor is the reference and the rest are
%   TDOA calculations between the first and the corresponding sensor from 2
%   to M, being M the number of sensors
%   - sensors = matrix Mx2 with the positions of the sensors
%   
%   Output:
%   x, y = Estimation of the position of the unknown transmitter


% If we have 3 sensors, solution is exact (or there is none)
if length(sensors) == 3
    % This method is based on Fang's method:
    % https://ieeexplore.ieee.org/document/102710
    [x, y] = exactcomputation(doa_array, sensors);
   
    
% If we have 4 sensors, we can add a new variable and compute an exact
% solution
elseif length(sensors) == 4
    [A, b] = getMatrices(doa_array, sensors);
    result = A\b;
    x = result(2);
    y = result(3);
 

% Same as with 4 sensors, but we compute the least squares
elseif length(sensors) > 4
    x = 0;
    y = 0;
    
else
    error('Not enough input sensors');
end

end


%% Fang's method
function [x, y] = exactcomputation(doa_array, sensors)

% The core of this method is that reference is located in (0, 0)
% The second sensor is located on (b,0)
% For that we need to rotate the vectors
vectortoref = sensors(2:end,:) - sensors(1,:);

% With the basis line
basisline = vectortoref(1,:);
angle = atan2(basisline(2), basisline(1));
rotmat_0a = [cos(angle), -sin(angle); sin(angle), cos(angle)]; % Matrix to rotate from a to 0

% Now we put the sensors in the reference basis
sensorb = rotmat_0a' * basisline';
sensorc = rotmat_0a' * vectortoref(2,:)';

% We extract the coefficients
b = norm(sensorb);
c = norm(sensorc);
cx = sensorc(1);
cy = sensorc(2);

% We get the TDOA values
R12 = doa_array(1);
R13 = doa_array(2);

% We obtain the equation coefficients
% y = g*x + h
g = ((R13/R12)*b - cx)/cy;
h = (c^2 - R13^2 + R12*R13*(1 - (b/R12)^2))/(2*cy);

% 0 = d*x^2 + e*x + f
d = -(1 + g^2 - (b/R12)^2);
e = b*(1-(b/R12)^2) - 2*g*h;
f = (R12^2/4)*(1-(b/R12)^2)^2 - h^2;

% Now we obtain the positive term x and y
xrotplus = (-e + sqrt(e^2 - 4*d*f))/(2*d);
yrotplus = g*xrotplus + h;

% We rotate the values
resultplus = rotmat_0a * [xrotplus; yrotplus];

% We add translate by the location of the sensor
xplus = resultplus(1) + sensors(1,1);
yplus = resultplus(2) + sensors(1,2);

% We check the doa and see if the sign is equal to the doa_array
plusvector = [xplus, yplus];
plusvectorsign = sign(norm(plusvector - sensors(1,:)) - norm(plusvector - sensors(2,:)));
if plusvectorsign == sign(doa_array(1))
    x = xplus;
    y = yplus;
end

% Now we obtain the negative term x and y
xrotminus = (-e - sqrt(e^2 - 4*d*f))/(2*d);
yrotminus = g*xrotminus + h;

% We rotate the values
resultminus = rotmat_0a * [xrotminus; yrotminus];

% We add translate by the location of the sensor
xminus = resultminus(1) + sensors(1,1);
yminus = resultminus(2) + sensors(1,2);

% We check the doa and see if the sign is equal to the doa_array
minusvector = [xminus, yminus];
minusvectorsign = sign(norm(minusvector - sensors(1,:)) - norm(minusvector - sensors(2,:)));
if minusvectorsign == sign(doa_array(1))
    if exist('x','var') == 1
        x = [x, xminus];
        y = [y, yminus];
    else
        x = xminus;
        y = yminus;
    end
end

end

%% Get matrices

