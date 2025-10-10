function [insert_dec, insert_obj] = selectIndividuals(pop_dec, pop_obj, FrontNo, desiredCount)
    % selectIndividuals - 从个体中选择满足条件的个体和适应度值

    % 初始化，提取 FrontNo 为 1 的个体
    insert_dec = pop_dec(FrontNo == 1, :);
    insert_obj = pop_obj(FrontNo == 1, :);

    if size(insert_dec, 1) >= desiredCount
        insert_dec = insert_dec(1:desiredCount, :);
        insert_obj = insert_obj(1:desiredCount, :);
        return;
    end

    uniqueFrontNos = unique(FrontNo);
    uniqueFrontNos(uniqueFrontNos == 1) = [];

    for i = 1:length(uniqueFrontNos)
        if size(insert_dec, 1) >= desiredCount
            break;
        end

        indices = find(FrontNo == uniqueFrontNos(i));
        candidates = pop_dec(indices, :);
        candidates_obj = pop_obj(indices, :);

        if size(candidates, 1) <= 3
            insert_dec = [insert_dec; candidates];
            insert_obj = [insert_obj; candidates_obj];
            continue;
        end

        P1 = candidates_obj(1, :);
        P2 = candidates_obj(end, :);
        A = P2(2) - P1(2);
        B = P1(1) - P2(1);
        C = P2(1) * P1(2) - P1(1) * P2(2);

        distances = abs(A * candidates_obj(:, 1) + B * candidates_obj(:, 2) + C) / sqrt(A^2 + B^2);
        [~, max_index] = max(distances);

        insert_dec = [insert_dec; candidates(max_index, :)];
        insert_obj = [insert_obj; candidates_obj(max_index, :)];

        if size(insert_dec, 1) >= desiredCount
            break;
        end

        left_indices = find(A * candidates_obj(:, 1) + B * candidates_obj(:, 2) + C > 0);
        right_indices = find(A * candidates_obj(:, 1) + B * candidates_obj(:, 2) + C < 0);
        left_obj = candidates_obj(left_indices, :);
        right_obj = candidates_obj(right_indices, :);

        left_crowding_distance = computeCrowdingDistance(left_obj);
        right_crowding_distance = computeCrowdingDistance(right_obj);

        left_max_point_dec = [];
        right_max_point_dec = [];

        if ~isempty(left_crowding_distance)
            [~, left_max_idx] = max(left_crowding_distance);
            left_max_point_dec = candidates(left_indices(left_max_idx), :);
            left_max_point_obj = left_obj(left_max_idx, :);
        end

        if ~isempty(right_crowding_distance)
            [~, right_max_idx] = max(right_crowding_distance);
            right_max_point_dec = candidates(right_indices(right_max_idx), :);
            right_max_point_obj = right_obj(right_max_idx, :);
        end

        remaining_slots = desiredCount - size(insert_dec, 1);

        if remaining_slots >= 2
            if ~isempty(left_max_point_dec)
                insert_dec = [insert_dec; left_max_point_dec];
                insert_obj = [insert_obj; left_max_point_obj];
            end
            if ~isempty(right_max_point_dec)
                insert_dec = [insert_dec; right_max_point_dec];
                insert_obj = [insert_obj; right_max_point_obj];
            end
        elseif remaining_slots == 1
            if ~isempty(left_max_point_dec) && ~isempty(right_max_point_dec)
                if rand() > 0.5
                    insert_dec = [insert_dec; left_max_point_dec];
                    insert_obj = [insert_obj; left_max_point_obj];
                else
                    insert_dec = [insert_dec; right_max_point_dec];
                    insert_obj = [insert_obj; right_max_point_obj];
                end
            elseif ~isempty(left_max_point_dec)
                insert_dec = [insert_dec; left_max_point_dec];
                insert_obj = [insert_obj; left_max_point_obj];
            elseif ~isempty(right_max_point_dec)
                insert_dec = [insert_dec; right_max_point_dec];
                insert_obj = [insert_obj; right_max_point_obj];
            end
        end
    end

    if size(insert_dec, 1) > desiredCount
        insert_dec = insert_dec(1:desiredCount, :);
        insert_obj = insert_obj(1:desiredCount, :);
    end
end

function crowding_distance = computeCrowdingDistance(obj_values)
    n = size(obj_values, 1);
    if n <= 2
        crowding_distance = inf(n, 1);
        return;
    end

    M = size(obj_values, 2);
    crowding_distance = zeros(n, 1);

    for m = 1:M
        [~, sorted_idx] = sort(obj_values(:, m));
        crowding_distance(sorted_idx(1)) = inf;
        crowding_distance(sorted_idx(end)) = inf;
        for j = 2:n - 1
            crowding_distance(sorted_idx(j)) = crowding_distance(sorted_idx(j)) + ...
                (obj_values(sorted_idx(j + 1), m) - obj_values(sorted_idx(j - 1), m));
        end
    end
end
