function state = evaluate_population_evolution(pareto_front, decision_vars, ref_pareto, ref_decision, localFE, t)
    % Evaluate the dynamic evolution state of a population by combining
    % Pareto front change rate and decision variable change rate, with task-specific history.

    % Inputs:
    %   pareto_front  : Current objective values of the population (N x M)
    %   decision_vars : Current decision variables of the population (N x D)
    %   ref_pareto    : Previous Pareto front (N_ref x M)
    %   ref_decision  : Previous decision variables (N_ref x D)
    %   localFE       : Local function evaluations for this task
    %   t             : Task identifier (unique for each task)

    % Persistent maps for task-specific HV and decision variable history
    persistent hv_last_map decision_last_map;

    % Initialize the maps if they are empty
    if isempty(hv_last_map)
        hv_last_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
        decision_last_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
    end

    % Initialize task-specific values if not already present
    if ~isKey(hv_last_map, t)
        hv_last_map(t) = Inf; % Use Inf to ensure proper comparison for the first iteration
    end
    if ~isKey(decision_last_map, t)
        decision_last_map(t) = Inf;
    end

    % Return if the localFE is 0 (initial state)
    if localFE == 0
        state = 'Balanced';
        fprintf('Population Evolution State (Task %d): %s\n', t, state);
        return;
    end

    % Reference point for HV calculation
    refpoint = max(pareto_front, [], 1) + 10 * ones(1, size(pareto_front, 2)); 
    hv_current = calculate_hypervolume(pareto_front, refpoint);  
    hv_ref = calculate_hypervolume(ref_pareto, refpoint);
    hv_change_rate = abs(hv_current - hv_ref);
%     disp(hv_change_rate);

    % Calculate decision variable change rate (using Wasserstein distance)
    wassdistance = computeWasserstein(decision_vars, ref_decision);
    decision_range = sqrt(sum((max(decision_vars) - min(decision_vars)).^2));
    normalized_decision_change = wassdistance / decision_range;
%     disp(normalized_decision_change);

    % Retrieve last generation's history
    hv_last = hv_last_map(t);
%     disp(hv_last);
    decision_last = decision_last_map(t);
%     disp(decision_last);

    % Determine evolution state by comparing with last generation
    if mean(hv_change_rate) <= hv_last && mean(normalized_decision_change) <= decision_last
        state = 'Converging';
    elseif mean(normalized_decision_change) >= decision_last || mean(hv_change_rate) >= hv_last
        state = 'Diverse';
    else
        state = 'Balanced';
    end

    % Update persistent variables to store current generation's values
    hv_last_map(t) = mean(hv_change_rate);
    decision_last_map(t) = mean(normalized_decision_change);

    % Display state
    fprintf('Population Evolution State (Task %d): %s\n', t, state);
end
% function state = evaluate_population_evolution(pareto_front, decision_vars, ref_pareto, ref_decision, localFE, t)
%     % Evaluate the dynamic evolution state of a population by combining
%     % Pareto front change rate and decision variable change rate, with task-specific history.
% 
%     % Inputs:
%     %   pareto_front  : Current objective values of the population (N x M)
%     %   decision_vars : Current decision variables of the population (N x D)
%     %   ref_pareto    : Previous Pareto front (N_ref x M)
%     %   ref_decision  : Previous decision variables (N_ref x D)
%     %   localFE       : Local function evaluations for this task
%     %   t             : Task identifier (unique for each task)
% 
%     % Persistent variables for task-specific HV and decision change history
%     persistent hv_change_history;
%     persistent decision_change_history;
% 
%     % Initialize the maps if they are empty
%     if isempty(hv_change_history)
%         hv_change_history = containers.Map('KeyType', 'double', 'ValueType', 'any');
%     end
%     if isempty(decision_change_history)
%         decision_change_history = containers.Map('KeyType', 'double', 'ValueType', 'any');
%     end
% 
%     % Return if the localFE is 0 (initial state)
%     if localFE == 0
%         state = 'Diverse';
%         fprintf('Population Evolution State (Task %d): %s\n', t, state);
%         return;
%     end
% 
%     % Reference point for HV calculation
%     refpoint = max(pareto_front, [], 1) + 10 * ones(1, size(pareto_front, 2)); 
%     hv_current = calculate_hypervolume(pareto_front, refpoint);  
%     hv_ref = calculate_hypervolume(ref_pareto, refpoint);
%     hv_change_rate = abs(hv_current - hv_ref);
%     disp(hv_change_rate);
% 
%     % Calculate decision variable change rate (using Wasserstein distance)
%     wassdistance = computeWasserstein(decision_vars, ref_decision);
%     decision_range = sqrt(sum((max(decision_vars) - min(decision_vars)).^2));
%     normalized_decision_change = wassdistance / decision_range;
%     disp(normalized_decision_change);
%     
%     % Update task-specific history (only store the most recent result for each task)
%     hv_change_history(t) = mean(hv_change_rate);
%     decision_change_history(t) = mean(normalized_decision_change);
% 
%     % Retrieve the latest values for this task
%     hv_last = hv_change_history(t);
%     disp(hv_last);
%     decision_last = decision_change_history(t);
%     disp(decision_last);
% 
%     % Determine evolution state
%     if mean(hv_change_rate) <= hv_last && mean(normalized_decision_change) <= decision_last
%         state = 'Converging';
%     elseif mean(normalized_decision_change) >= decision_last || mean(hv_change_rate) >= hv_last
%         state = 'Diverse';
%     else
%         state = 'Balanced';
%     end
% 
%     % Display state
%     fprintf('Population Evolution State (Task %d): %s\n', t, state);
% end


