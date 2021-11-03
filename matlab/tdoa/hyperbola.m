function [ hyp_array ] = hyperbola(doa_meters, rx1, rx2)
%hyperbola: calculates the points of a hyperbola from receiver positions
%               and the doa in meters
  
    
    % Get the coordinates
    rx1_x = rx1(1); rx1_y = rx1(2);
    rx2_x = rx2(1); rx2_y = rx2(2);


    % mit Kosinussatz Dreieck berechnen
    rx_x_dist = rx2_x - rx1_x;
    rx_y_dist = rx2_y - rx1_y;

    dist_12 = abs (rx_x_dist+1i*rx_y_dist);  % positions in complex plane
    angle_12 = angle (rx_x_dist+1i*rx_y_dist); % -pi to +pi

    hyp_x = zeros(1,1);
    hyp_y = zeros(1,1);

    hyp_x_leg1 = zeros(1,1);
    hyp_y_leg1 = zeros(1,1);
    hyp_x_leg2 = zeros(1,1);
    hyp_y_leg2 = zeros(1,1);
    hyp_point_counter = 0;
    
    if abs(doa_meters) > dist_12
        disp(['<strong>TDOA delay (' num2str(doa_meters) ' meters) larger than RX distance (' num2str(dist_12) ' meters) -> no solution possible </strong>']);
        doa_meters = sign(doa_meters) * 0.995 * dist_12;
        disp(['<strong>ATTENTION: Correcting TODA delay to 0.995 * RX distance (maximum possible value) = ' num2str(0.995*doa_meters) '</strong>']);
    end
        
        
    if abs(doa_meters) <= dist_12

        %for r_1 = (exp(0:0.05:4)-1) / 5
        for r_1 = 0:5:1000000
            r_2 = r_1 - doa_meters;
            %disp(['r_1 = ' num2str(r_1) ', r_2 = ' num2str(r_2)]);

            if ((r_2 + r_1) > dist_12)  % checks if triangle can be created
                
				acos_argument = (r_2^2 - r_1^2 - dist_12^2) / (-2*r_1*dist_12);
				
				if (acos_argument >= -1) && (acos_argument <= +1) % checks if triangle can be created
				
					hyp_point_counter = hyp_point_counter + 1;

					hyp_angle = acos(acos_argument); % inner angle of triangle at RX1
                
					abs_angle1 = wrap2pi(angle_12 + hyp_angle);  % 1st solution: hyperbola leg 1
					hyp_x_leg1(hyp_point_counter) = rx1_x + r_1 * cos(abs_angle1);
					hyp_y_leg1(hyp_point_counter) = rx1_y + r_1 * sin(abs_angle1);

					abs_angle2 = wrap2pi(angle_12 - hyp_angle);  % 2nd solution: hyperbola leg 2 
					hyp_x_leg2(hyp_point_counter) = rx1_x + r_1 * cos(abs_angle2);
					hyp_y_leg2(hyp_point_counter) = rx1_y + r_1 * sin(abs_angle2);
                else
                    %disp(['acos argument ' num2str(acos_argument)]);
				end
            end

        end
    else
        disp('TDOA delay larger than RX distance -> no solution possible');
    end
    
    if (hyp_point_counter == 0)
        disp('Hyperbola could not be constructed');
    end

    hyp_x = [fliplr(hyp_x_leg1) hyp_x_leg2];
    hyp_y = [fliplr(hyp_y_leg1) hyp_y_leg2];
    hyp_points = 2* hyp_point_counter;
    
    hyp_array = [hyp_x; hyp_y];
    
    disp(['Hyperbola with totally ' num2str(hyp_points) ' points generated.']);
end