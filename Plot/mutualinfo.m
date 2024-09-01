function MI = mutualinfo(selec_feat_discrt, Combined_Labels)
    % This function calculates the mutual information between a set of 
    % discretized features and the labels.
    
    n = size(selec_feat_discrt, 1);
    Y = unique(Combined_Labels);
    values_of_feat = [];
    P = [];
    feat_ind = 0;

    % Calculate joint probabilities
    for r = 1:n
        current_feat = selec_feat_discrt(r, :);
        if ~is_x_in_values_of_x(current_feat, values_of_feat)
            feat_ind = feat_ind + 1;
            values_of_feat = [values_of_feat; current_feat];
            P = [P; zeros(1, length(Y))];
            for rr = r:n
                if isequal(selec_feat_discrt(rr, :), current_feat)
                    y_ind = find(Y == Combined_Labels(rr));
                    P(feat_ind, y_ind) = P(feat_ind, y_ind) + 1;
                end
            end
        end
    end

    % Calculate mutual information
    P(P == 0) = eps;  % Replace 0s with a small epsilon to avoid log(0)
    sum1 = sum(sum(P .* log2(P)));
    sum2 = sum(sum(P, 2) .* log2(sum(P, 2)));
    sum3 = sum(sum(P, 1) .* log2(sum(P, 1)));
    MI = log2(n) + (sum1 - sum2 - sum3) / n;
end

function j = is_x_in_values_of_x(x, values_of_x)
    % This helper function checks if x is part of the set values_of_x.
    % j = true if x is found in values_of_x, false otherwise.
    
    j = false;
    if isempty(values_of_x)
        return;
    end
    
    for r = 1:size(values_of_x, 1)
        if isequal(x, values_of_x(r, :))
            j = true;
            break;
        end
    end
end