function hyp = calculate_hypervolume(pf, repoint)
    % Hypervolume calculation in MATLAB
    % Inputs:
    %   pf      : Pareto front (N x 2 matrix, where N is the number of points)
    %   repoint : Reference point (1 x 2 vector)
    % Outputs:
    %   hyp     : Hypervolume value

    % Sort the Pareto front based on the first objective
    [~, temp_index] = sort(pf(:, 1));
    sorted_pf = pf(temp_index, :);

    % Combine the reference point and sorted Pareto front
    pointset = [repoint; sorted_pf];

    % Initialize hypervolume
    hyp = 0;

    % Compute hypervolume
    for i = 1:size(pf, 1)
        cubei = (pointset(1, 1) - pointset(i + 1, 1)) * (pointset(i, 2) - pointset(i + 1, 2));
        hyp = hyp + cubei;
    end
end
