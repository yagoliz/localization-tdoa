function [position] = trilateration_nlls(n, pos, starting_position, dist, W)
    
    size_dist = size(dist);
    if size_dist(1) == 1
        dist = dist';
    end
    position = starting_position;
    min_norm = 1e-3;
    max_it = 500;
    frac = 1/10;
    stop = 0;
    it = 0;
    while stop == 0
        f = zeros(n,1);
        J = zeros(size(pos,1),2);
        for i = 1:n
            temp_val = sqrt((position(1)-pos(i,1))^2 + (position(2)-pos(i,2))^2);
            f(i) = temp_val - dist(i);
            J(i,1) = (position(1) - pos(i,1)) / temp_val;
            J(i,2) = (position(2) - pos(i,2)) / temp_val;
        end
%         # Use linsolve which is more stable, use correct formula for weights
%         #delta_est = (J'*W*J)^-1*(J'*W*f);
        delta_est = linsolve(J'*W*W*J, J'*W*f);
        position = position - frac*delta_est;
        it = it + 1;
        if isnan(position)
            stop = 1;
            position = starting_position;
%             disp(' isnan');
        elseif norm(delta_est) < min_norm
            stop = 1;
%             disp('norm<min_norm');
        elseif it == max_it
            stop = 1;
%             disp('it==max_it');
        end
        
    end
end
