function [position] = calculate_position(vector_distances, anchor_positions, num_least_anchors, nlls_weight, nlls_scale, previous_position)

    if length(anchor_positions) >= 3

        num_available_anchors = length(vector_distances);

        if num_available_anchors >= num_least_anchors
            if nlls_weight == 0
                w = ones(1,length(vector_distances));
            elseif nlls_weight == 1
                w = 1./(1 + vector_distances);
            elseif nlls_weight == 2
                factor = .09;
                w = 1 - factor * (vector_distances - min(vector_distances));
                w = 1./w;
            elseif nlls_weight == 5
                if isempty(previous_position)
                    w = ones(1, length(vector_distances));
                else
                    matrix_pdist = pdist2(anchor_positions, repmat(previous_position', length(anchor_positions), 1));
                    w = matrix_pdist(:, 1);
                    w = w - min(w);
                    if max(w) > .01
                        w = w/max(w);
                    end
                    w = (1./(1+w))';
                end

            end
            
            W = diag(w) * nlls_scale;

            if isempty(previous_position)
                position = trilateration_lls(num_available_anchors, anchor_positions, vector_distances, W);
            else
                position = trilateration_nlls(num_available_anchors, anchor_positions, previous_position, vector_distances, W);
            end

        else
            disp(' Too few anchors to make a position estimate');
            position = NaN;
        end
    else
        disp(' At least 3 anchors please!');
        position = NaN;

    end

end

