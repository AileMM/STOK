classdef SGETO < Algorithm
% <Multi-task> <Multi-objective> <None/Constrained>
% COPY EMT-ET


properties (SetAccess = private)
    G = 8  %人为设定的选择8个待迁移的个体
    P = 0.5
    MuC = 20
    MuM = 15
end

methods
    function Parameter = getParameter(Algo)
        Parameter = {'G: Transfer Solutions Number', num2str(Algo.G), ...
                'P: Distrib Probability', num2str(Algo.P), ...
                'MuC: Simulated Binary Crossover', num2str(Algo.MuC), ...
                'MuM: Polynomial Mutation', num2str(Algo.MuM)};
    end

    function setParameter(Algo, Parameter)
        i = 1;
        Algo.G = str2double(Parameter{i}); i = i + 1;
        Algo.P = str2double(Parameter{i}); i = i + 1;
        Algo.MuC = str2double(Parameter{i}); i = i + 1;
        Algo.MuM = str2double(Parameter{i}); i = i + 1;
    end

    function run(Algo, Prob)
        
%         previous_similarity = inf;
        % Initialize    初始化不同任务的不同种群和外部存档
        population = Initialization(Algo, Prob, Individual_ETO);
        Archive = Initial_arc(Algo, Prob, Individual_ETO);  
        
        convergence_pool = cell(1, Prob.T);
        diversity_pool = cell(1, Prob.T);
        new_pool = cell(1, Prob.T);
        
%         convergence_pool = [];
%         diversity_pool = [];
        
        convergence_pop = cell(1, Prob.T);
        diversity_pop = cell(1, Prob.T);
        
        mind_insert_Arc = cell(1, Prob.T);
        insert_Arc = cell(1, Prob.T);
        
        current_obj = cell(1, Prob.T); % 每个任务的当前目标值
        current_dec = cell(1, Prob.T); % 每个任务的当前决策变量
        ref_obj = cell(1, Prob.T);     % 参考目标值
        ref_dec = cell(1, Prob.T);     % 参考决策变量
        pre_obj = cell(1, Prob.T);     % 中间变量
        pre_dec = cell(1, Prob.T);     % 中间变量
        
        surrogate_pool = cell(1, Prob.T);
        
        localFE = zeros(1, Prob.T); % 每个任务的本地FE计数器
        state = cell(1, Prob.T);  % 每个任务的状态
        
        for t = 1:Prob.T
            [rank, FrontNo] = NSGA2Sort(population{t});
            for i = 1:length(population{t})
                population{t}(i).FrontNo = FrontNo(i); 
            end
            population{t} = population{t}(rank);
        end
        
        for t = 1:Prob.T
            [rank, FrontNo] = NSGA2Sort(Archive{t});
            for i = 1:length(Archive{t})
                Archive{t}(i).FrontNo = FrontNo(i); 
            end
            Archive{t} = Archive{t}(rank);
         
            %% 这里的参考值到底是Archive
            Arc_obj = vertcat(Archive{t}.Obj);
            Arc_decs = vertcat(Archive{t}.Dec);
            
            ref_obj{t} = Arc_obj;  
            ref_dec{t} = Arc_decs;

            pre_obj{t} = ref_obj{t};
            pre_dec{t} = ref_dec{t};
        end
        
        
        while Algo.notTerminated(Prob, Archive)
             %% 插值点个数
            desiredCount = 5;
            
            for t = 1:Prob.T
               % 构建代理模型 这里得到的是每个目标的代理模型，构建应该没有问题
                sur_mtgp = buildMTGP(vertcat(Archive{t}.Dec), vertcat(Archive{t}.Obj), Prob.M(t));
                surrogate_pool{t}=sur_mtgp; 
            end
            for t = 1:Prob.T
                % 模型预测这里应该是有问题的
                X = vertcat(population{t}.Dec);
                % Predict mean and covariance for H
                H_pred = zeros(Prob.N, Prob.M(t)); % Predicted mean for H
                Sigma_H = cell(1, Prob.M(t)); % Predicted covariance for H
                Sigma_H_combined = zeros(Prob.M(t), Prob.M(t));
                C = rand(Prob.M(t), Prob.M(t)) * 2;
                while det(C) == 0
                    C = rand(Prob.M(t), Prob.M(t)) * 2;
                end
                % 3. 确保每一列都是线性独立的，可以使用SVD方法来检查并重构C
                [U, S, V] = svd(C);
                threshold = 1e-6;  % 设置一个阈值来剔除小的奇异值
                S(S < threshold) = 0;
                C = U * S * V';  % 重构C，确保满秩
                for i = 1:Prob.M(t)
                    [H_pred(:, i),Sigma_H{i},~] = predict(surrogate_pool{t}{i}, X, 'Alpha', 0.05);
                end
                %% Step 4: Recover original F from H using C^-1
                C_inv = pinv(C);  % 使用伪逆代替普通求逆
                % Recover F (mean prediction)
                F_pred = (C_inv * H_pred')';
                % Recover F covariance matrix
                F_cov = cell(1, Prob.N); % Covariance matrices for each test point
                for j = 1:Prob.N
                    for i = 1:Prob.M(t)
                        for k = 1:Prob.M(t)
                            % 计算目标 i 和目标 k 之间的协方差
                            cov_matrix = cov(Sigma_H{i}', Sigma_H{k}');  % 计算协方差矩阵
                            % 提取目标 i 和目标 k 之间的协方差
                            Sigma_H_combined(i, k) = cov_matrix(1, 2);  % Cov(X, Y) 即目标 i 和目标 k 的协方差
                        end
                    end
                    % 确保协方差矩阵正定
                    epsilon = 1e-6;  % 微小正数
                    Sigma_H_combined = Sigma_H_combined + epsilon * eye(Prob.M(t));  % 使矩阵正定
                    % Transform covariance for F
                    F_cov{j} = C_inv * Sigma_H_combined * C_inv';
                end
                % Extract standard deviations for F (diagonal of covariance matrix)
                F_std = zeros(Prob.N, Prob.M(t));
                for j = 1:Prob.N
                    F_std(j, :) = sqrt(diag(F_cov{j}));
                end
                LCB{t} = F_pred - 2 * F_std; 
                
                %% 这里给到种群的模型的预测值，其实给不给都无所谓 只要划分按照这个来就行
                for i = 1:length(population{t})
                    population{t}(i).Obj = LCB{t}(i,:); 
                end
              %% 后面依赖LCB预测值划分种群！
                [rank, FrontNo] = NSGA2Sort(population{t});
                [convergence_indices, diversity_indices] = classifyByFrontAndRank(FrontNo, rank);

                convergence_pop{t} = population{t}(convergence_indices);
                diversity_pop{t} = population{t}(diversity_indices);
                % 保留之前的划分的个体
                convergence_pool{t} = [convergence_pool{t},convergence_pop{t}];
                diversity_pool{t} = [diversity_pool{t},diversity_pop{t}];

                X = vertcat(Archive{t}.Dec);
                Y = vertcat(Archive{t}.Obj);
                
                current_dec{t} = X; 
                current_obj{t} = Y;
                
            end
                  
            for t = 1:Prob.T 
                state{t} = evaluate_population_evolution(current_obj{t}, current_dec{t}, pre_obj{t}, pre_dec{t}, localFE(t), t);
            end
%             g = 0;
%             while g<= 20
            for t = 1:Prob.T
               offspring{t} = Algo.Generation(Prob,population, t, convergence_pool,diversity_pool,surrogate_pool{t},state);
               population{t} = [population{t}, offspring{t}];
                [rank, FrontNo] = NSGA2Sort(population{t});
                for i = 1:length(population{t})
                    population{t}(i).FrontNo = FrontNo(i);
                end
                [population{t}, ~] = Selection_SPEA2(population{t}, Prob.N);
            end
%                  g = g + 1;
%             end
            
            for t = 1:Prob.T
       
                [~, FrontNo] = NSGA2Sort(population{t});
                mid_population{t} = population{t}(FrontNo == 1);
                insert_rankD = vertcat(mid_population{t}.Dec);
                insert_rankO = vertcat(mid_population{t}.Obj);
                
                %% 中间过程 填充个体的选择 得到的种群中选择
                temp_individual = 2;
                
               if size(insert_rankD,1)< temp_individual
                   [insert_dec, insert_obj] = selectIndividuals(vertcat(population{t}.Dec), vertcat(population{t}.Obj), FrontNo, temp_individual);
               else
                   index = randperm(size(insert_rankD,1),temp_individual);
                   insert_dec = insert_rankD(index,:);
                   insert_obj = insert_rankO(index,:);
               end

               mind_insert_Arc{t} = Individual_ETO();
               for j = 1:size(insert_dec,1)
                    mind_insert_Arc{t}(j).Dec = insert_dec(j,:);
                    mind_insert_Arc{t}(j).Obj = insert_obj(j,:);
               end
            end
            
            %% 更新存档，更新代理模型的构建种群
            for t = 1:Prob.T
                [nondominated_set, ~] = checkNonDominated(Archive{t},mind_insert_Arc{t});
                insert_Arc{t} = Individual_ETO();
                if isempty(nondominated_set)
                    if strcmp(state{t}, 'Converging')
                        row = length(diversity_pool{t});
                        if row==0 || row<num_additional_individuals
                            selected_individuals = Individual_ETO();
                            for i = 1:num_additional_individuals
                                selected_individuals(i).Dec = (Prob.Ub{Prob.T} - Prob.Lb{Prob.T}) .* rand(1, max(Prob.D)) + Prob.Lb{Prob.T};
                            end
                        else
                            random_indices = randperm(row, desiredCount);  
                            selected_individuals = diversity_pool{t}(random_indices);
%                         diversity_pool{t}(random_indices) = [];
                        end
                    elseif strcmp(state{t}, 'Diverse') 
                        row = length(convergence_pool{t});
                        random_indices = randperm(row, desiredCount);  
                        selected_individuals = convergence_pool{t}(random_indices);
%                         convergence_pool{t}(random_indices) = [];
                    else 
                        new_pool{t} = [convergence_pool{t},diversity_pool{t}];
                        row = length(new_pool{t});
                        random_indices = randperm(row, desiredCount);  
                        selected_individuals = new_pool{t}(random_indices);
%                         new_pool{t}(random_indices) = [];
                    end
                    insert_Arc{t} = selected_individuals;
                  
                 elseif size(nondominated_set, 1) == desiredCount
                    for i = 1:size(nondominated_set, 1)
                        insert_Arc{t}(i).Dec = nondominated_set(i, :);
                    end
                    
                 elseif size(nondominated_set, 1) < desiredCount
                   num_additional_individuals = desiredCount - size(nondominated_set, 1);
                   if strcmp(state{t}, 'Converging')
                        row = length(diversity_pool{t});
                        if row==0 || row<num_additional_individuals
                            selected_individuals = Individual_ETO();
                            for i = 1:num_additional_individuals
                                selected_individuals(i).Dec = (Prob.Ub{Prob.T} - Prob.Lb{Prob.T}) .* rand(1, max(Prob.D)) + Prob.Lb{Prob.T};
                            end
                        else
                            random_indices = randperm(row, num_additional_individuals);  
                            selected_individuals = diversity_pool{t}(random_indices);
    %                         diversity_pool{t}(random_indices) = [];
                        end

                    elseif strcmp(state{t}, 'Diverse') 
                        row = length(convergence_pool{t});
                        random_indices = randperm(row, num_additional_individuals);  
                        selected_individuals = convergence_pool{t}(random_indices);
%                         convergence_pool{t}(random_indices) = [];

                    else 
                        new_pool{t} = [convergence_pool{t}, diversity_pool{t}];
                        row = length(new_pool{t});
                        random_indices = randperm(row, num_additional_individuals);  
                        selected_individuals = new_pool{t}(random_indices);
%                         new_pool{t}(random_indices) = [];

                   end
                   for i = 1:size(nondominated_set, 1)
                        insert_Arc{t}(i).Dec = nondominated_set(i, :);
                   end
                   insert_Arc{t} = [insert_Arc{t},selected_individuals];
                   
                 else 
                   index = randperm(size(nondominated_set, 1), desiredCount);  
                    for i = 1:desiredCount
                        idx = index(i); 
                        insert_Arc{t}(i).Dec = nondominated_set(idx, :); % Assign decision variables
                    end
                end
                insert_Arc{t} = Algo.Evaluation(insert_Arc{t}, Prob, t);
                Archive{t} = [Archive{t}, insert_Arc{t}];
                [~,FrontNo] = NSGA2Sort(Archive{t});
                Archive{t} = Archive{t}(FrontNo == 1);
                
                localFE(t) = localFE(t) + desiredCount;
                
                  %% 这里应该与插值点的个体有关
                if  localFE(t) ~= desiredCount
                    pre_obj{t} = current_obj{t};  % 更新上一代的目标值
                    pre_dec{t} = current_dec{t};  % 更新上一代的决策变量
                end
                current_obj{t} = vertcat(Archive{t}.Obj);
                current_dec{t} = vertcat(Archive{t}.Dec);
 
            end
            
%             for t = 1:Prob.T
% %                 population{t} = [population{t}, offspring{t}];
% %                 [rank, FrontNo] = NSGA2Sort(population{t});
% %                 for i = 1:length(population{t})
% %                     population{t}(i).FrontNo = FrontNo(i);
% %                 end
% %                 population{t} = population{t}(rank(1:Prob.N));
%       
%                 %% 这里应该与插值点的个体有关
%                 if  localFE(t) ~= desiredCount
%                     pre_obj{t} = current_obj{t};  % 更新上一代的目标值
%                     pre_dec{t} = current_dec{t};  % 更新上一代的决策变量
%                 end
%                 current_obj{t} = vertcat(Archive{t}.Obj);
%                 current_dec{t} = vertcat(Archive{t}.Dec);
%                 
%             end
        end
    end



    
    function offspring = Generation(Algo, Prob, population, t, convergence_pool,diversity_pool,surrogate, state)
%         for j = 1:Prob.T
%             if t == j
%                 continue; 
%             end
%         end
%         j = 3 - t;
        candidates = setdiff(1:1:Prob.T, t); % 过滤掉 t，得到可能的 j 候选集
        j = candidates(randi(length(candidates))); % 随机选择一个不同于 t 的任务编号
%         similarity = kl_divergence_gaussian(model{t}.mean, model{t}.sigma, model{j}.mean, model{j}.sigma);
        if rand < 0.5
            %% 迁移个体数=5
            inject_pop = Individual_ETO.empty();
            num_individuals = 5;

            if strcmp(state{t}, 'Converging')
                row = length(diversity_pool{j});
                 num_individuals = min(num_individuals, row);
                random_indices = randperm(row, num_individuals);  
                selected_individuals = diversity_pool{j}(random_indices);
                diversity_pool{j}(random_indices) = [];
                
            elseif strcmp(state{t}, 'Diverse') 
                row = length(convergence_pool{j});
                num_individuals = min(num_individuals, row);
                random_indices = randperm(row, num_individuals);  
                selected_individuals = convergence_pool{j}(random_indices);
                convergence_pool{j}(random_indices) = [];
            else 
                new_pool{j} = [convergence_pool{j}, diversity_pool{j}];
                row = length(new_pool{j});
                random_indices = randperm(row, num_individuals);  
                selected_individuals = new_pool{j}(random_indices);
                new_pool{j}(random_indices) = [];
            end
            % 获得迁移个体
            inject_pop = Individual_ETO.empty();
            selected_individuals = vertcat(selected_individuals.Dec);
            for i = 1:size(selected_individuals, 1)
                c = Individual_ETO();
%                 c.Dec = [selected_individuals(i, :), rand(1, max(Prob.D) - Prob.D(t))];
                c.Dec = selected_individuals(i, :);
                c.Dec(c.Dec > 1) = 1;
                c.Dec(c.Dec < 0) = 0;
                inject_pop = [inject_pop, c];
            end

            count = 1;
            population = [population{t}, inject_pop];
            for i = 1:ceil(length(population) / 2)
                p1 = i; p2 = i + fix(length(population) / 2);
                offspring(count) = population(p1);
                offspring(count + 1) = population(p2);
                [offspring(count).Dec, offspring(count + 1).Dec] = GA_Crossover(population(p1).Dec, population(p2).Dec, Algo.MuC);
                offspring(count).Dec = GA_Mutation(offspring(count).Dec, Algo.MuM);
                offspring(count + 1).Dec = GA_Mutation(offspring(count + 1).Dec, Algo.MuM);
                for x = count:count + 1
                    offspring(x).Dec(offspring(x).Dec > 1) = 1;
                    offspring(x).Dec(offspring(x).Dec < 0) = 0;
                end
                count = count + 2;
            end
            % Predict H for each objective 这里要换成row 迁移的知识合并种群产生子代
           [row,~] = size(vertcat(offspring.Dec));
            H_pred = zeros(row, Prob.M(t)); % Predicted mean for H
            Sigma_H = cell(1, Prob.M(t)); % Predicted covariance for H
            Sigma_H_combined = zeros(Prob.M(t), Prob.M(t));
            C = rand(Prob.M(t), Prob.M(t)) * 2;
            while det(C) == 0
                C = rand(Prob.M(t), Prob.M(t)) * 2;
            end
            % 3. 确保每一列都是线性独立的，可以使用SVD方法来检查并重构C
            [U, S, V] = svd(C);
            threshold = 1e-6;  % 设置一个阈值来剔除小的奇异值
            S(S < threshold) = 0;
            C = U * S * V';  % 重构C，确保满秩
            for i = 1:Prob.M(t)
                [H_pred(:, i),Sigma_H{i},~] = predict(surrogate{i}, vertcat(offspring.Dec), 'Alpha', 0.05);
            end
            %% Step 4: Recover original F from H using C^-1
            C_inv = pinv(C);  % 使用伪逆代替普通求逆
            % Recover F (mean prediction)
            F_pred = (C_inv * H_pred')';
            % Recover F covariance matrix
            F_cov = cell(1, row); % Covariance matrices for each test point
            for j = 1:row
                % Combine covariance matrices for H
                for i = 1:Prob.M(t)
                    for k = 1:Prob.M(t)
                        % 计算目标 i 和目标 k 之间的协方差
                        cov_matrix = cov(Sigma_H{i}', Sigma_H{k}');  % 计算协方差矩阵
                        % 提取目标 i 和目标 k 之间的协方差
                        Sigma_H_combined(i, k) = cov_matrix(1, 2);  % Cov(X, Y) 即目标 i 和目标 k 的协方差
                    end
                end
                % 确保协方差矩阵正定
                epsilon = 1e-6;  % 微小正数
                Sigma_H_combined = Sigma_H_combined + epsilon * eye(Prob.M(t));  % 使矩阵正定
                % Transform covariance for F
                F_cov{j} = C_inv * Sigma_H_combined * C_inv';
            end
            % Extract standard deviations for F (diagonal of covariance matrix)
            F_std = zeros(row, Prob.M(t));
            for j = 1:row
                F_std(j, :) = sqrt(diag(F_cov{j}));
            end
            LCB = F_pred - 2 * F_std; 

           for i = 1:length(offspring)
               offspring(i).Obj = LCB(i,:); 
               offspring(i).CV = 0;
           end
        else
            count = 1;
            population = population{t};
            for i = 1:ceil(length(population) / 2)
                p1 = i; p2 = i + fix(length(population) / 2);
                offspring(count) = population(p1);
                offspring(count + 1) = population(p2);
                [offspring(count).Dec, offspring(count + 1).Dec] = GA_Crossover(population(p1).Dec, population(p2).Dec, Algo.MuC);
                offspring(count).Dec = GA_Mutation(offspring(count).Dec, Algo.MuM);
                offspring(count + 1).Dec = GA_Mutation(offspring(count + 1).Dec, Algo.MuM);
                for x = count:count + 1
                    offspring(x).Dec(offspring(x).Dec > 1) = 1;
                    offspring(x).Dec(offspring(x).Dec < 0) = 0;
                end
                count = count + 2;
            end
            % Predict H for each objective 这里要换成row 迁移的知识合并种群产生子代
            [row,~] = size(vertcat(offspring.Dec));
            H_pred = zeros(row, Prob.M(t)); % Predicted mean for H
            Sigma_H = cell(1, Prob.M(t)); % Predicted covariance for H
            Sigma_H_combined = zeros(Prob.M(t), Prob.M(t));
            C = rand(Prob.M(t), Prob.M(t)) * 2;
            while det(C) == 0
                C = rand(Prob.M(t), Prob.M(t)) * 2;
            end
            % 3. 确保每一列都是线性独立的，可以使用SVD方法来检查并重构C
            [U, S, V] = svd(C);
            threshold = 1e-6;  % 设置一个阈值来剔除小的奇异值
            S(S < threshold) = 0;
            C = U * S * V';  % 重构C，确保满秩
            for i = 1:Prob.M(t)
                [H_pred(:, i),Sigma_H{i},~] = predict(surrogate{i}, vertcat(offspring.Dec), 'Alpha', 0.05);
            end
            %% Step 4: Recover original F from H using C^-1
            C_inv = pinv(C);  % 使用伪逆代替普通求逆
            % Recover F (mean prediction)
            F_pred = (C_inv * H_pred')';
            % Recover F covariance matrix
            F_cov = cell(1, row); % Covariance matrices for each test point
            for j = 1:row
                for i = 1:Prob.M(t)
                    for k = 1:Prob.M(t)
                        % 计算目标 i 和目标 k 之间的协方差
                        cov_matrix = cov(Sigma_H{i}', Sigma_H{k}');  % 计算协方差矩阵
                        % 提取目标 i 和目标 k 之间的协方差
                        Sigma_H_combined(i, k) = cov_matrix(1, 2);  % Cov(X, Y) 即目标 i 和目标 k 的协方差
                    end
                end
                % 确保协方差矩阵正定
                epsilon = 1e-6;  % 微小正数
                Sigma_H_combined = Sigma_H_combined + epsilon * eye(Prob.M(t));  % 使矩阵正定
                % Transform covariance for F
                F_cov{j} = C_inv * Sigma_H_combined * C_inv';
            end
            % Extract standard deviations for F (diagonal of covariance matrix)
            F_std = zeros(row, Prob.M(t));
            for j = 1:row
                F_std(j, :) = sqrt(diag(F_cov{j}));
            end
            LCB = F_pred - 2 * F_std;   
           for i = 1:length(offspring)
               offspring(i).Obj = LCB(i,:); 
               offspring(i).CV = 0;
           end
        end
%         updated_similarity = similarity;
    end
end
end

