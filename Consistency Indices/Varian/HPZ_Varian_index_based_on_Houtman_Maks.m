function [VARIAN, Var_approx_ratio, Var_in_sample_residuals, one_minus_v, Var_out_of_sample_residuals, out_of_sample_one_minus_v] = HPZ_Varian_index_based_on_Houtman_Maks (expenditure, identical_choice, index_threshold, Varian_algorithm_settings, aggregators, out_of_sample_required, main_folder_for_results)

% NOTE! We assume all the observations passed to this function are involved
% in cycles. If you will pass an observation that is not involved in any
% cycle, this function will crash.
% It is the job of the "Varian_Manager" function to make sure no such
% observations are passed to this function.



%% Expanation about this function:
% This function transforms the problem of finding the Varian Index, to a
% problem of finding the Houtman-Maks Index, when each observation has a
% different weight/value (in comparison to the original Houtman-Maks, in
% which all observations have a weight/value of 1).
% This transformation transforms a Varian problem with n observations, to a
% weighted-Houtman-Maks problem with m observations, n <= m <= n(n-1).
%
% After performing this transformation, it calls a recursive function that
% solves the weighted-Houtman-Maks problem.
% 
% The best indicator for the Time Complexity of the weighted-Houtman-Maks
% algorithm is the number of observations.
% Therefore, if after the transformation the new number of observations is
% too big, we perform a relaxation of the Varian problem, that reduces the
% new number of observations. If the relaxation wasn't good enough, we use
% a stronger relaxation, then a stronger one, and so on, until the new
% number of observations is small enough.
%
% When we use a relaxation, we get a value (VARIAN), and we get a
% multiplier (Var_approx_ratio), such that (VARIAN*Var_approx_ratio) is the
% upper bound for the real value, and (VARIAN/Var_approx_ratio) is the
% lower bound for the real value.
%% End of Explanation




% these are used for the approximation option:

% we are allowed "unite"/"combine" values of (1-v) if the ratios between
% them is no more than approx term. e.g. if for some original (old) 
% observation there is 1-v = 1.5 and 1-v = 1.3, then if approx_term = 1.2
% we can combine them to one new observation, but if approx_term = 1.1 we can't. 
% initial value (no approximation; precise):
approx_term = 1.00;
previous_approx_term = approx_term;
% in each step that we try to find an approximation of given percentage difference, 
% increment the differnce by this (it updates (getting bigger) over iterations though)  
approx_step = 0.01;
% num of observations in the original problem
[old_obs_num, ~] = size(expenditure);
% if the new obs_num is more than limit_obs_num, we assume that running 
% time will be too long, then we perform a relaxation that reduces the new
% obs_num to no more than target_obs_num (the relaxation makes it an 
% approximation instead of completely exact calculation).
% why 2 different variables? (limit_obs_num & target_obs_num)
% because if we had only one variable, and say its value was 600, it would 
% be a waste if we calculated approximation for a subject that without
% approximation would have 601 new obs_num.
target_obs_num = Varian_algorithm_settings(1);
limit_obs_num_multiplier = 1.05;   % should be bigger or equal to 1
limit_obs_num = target_obs_num * limit_obs_num_multiplier;   % should be bigger or equal to target_obs_num
% it will be a problem if the target num of obs is lower than the original
% num of obs (infinite loop), also even if it is close to it, it will be a
% problem. hence this:
old_obs_num_multiplier = Varian_algorithm_settings(2);   % should be bigger or equal to 1, preferably at least 2
target_obs_num = max(target_obs_num , old_obs_num_multiplier * old_obs_num);
limit_obs_num = max(limit_obs_num , old_obs_num_multiplier * limit_obs_num_multiplier * old_obs_num);



% find SDRP (it's a bit wasteful to calculate GARP, but it is negligible in
% comparison to the complexity of calculating the Varian index
[~, ~, ~, ~, SDRP] = GARP_based_on_expenditures(expenditure, identical_choice, index_threshold);



% now we transform the Varian problem to a weighted Houtman-Maks problem

% The matrix RATIO has at the cell in the i'th row and the j'th column, the ratio between 
% the value of the bundle that was chosen in observation j given the prices of observation 
% i and the value of the bundle that was chosen in observation i 
RATIO = expenditure ./ (diag(expenditure) * ones(old_obs_num,1)');
one_minus_RATIO = 1 - RATIO;

% new number of observations in the equivalent weighted-HM problem
% (if it is too big, we will reduce it by some relaxation, thus calculating
% an approximation instead of exact value)
new_obs_num = sum(sum(SDRP));


approx_iteration = 1;  % counts iterations of the while loop
while approx_iteration == 1 || (approx_iteration == 2 && new_obs_num > limit_obs_num) || (approx_iteration >= 3 && new_obs_num > target_obs_num)
    
    old_to_new_indices = zeros(1, old_obs_num);
    new_to_old = zeros(1, new_obs_num);
    new_counter = 0;
    for old_counter = 1:old_obs_num   
        % find the relations and their ratios, and sort the ratios in descending order 
        relations = find(SDRP(old_counter , :));
        v_weights = one_minus_RATIO(old_counter , relations);
        [v_weights , indices_for_sorting] = sort(v_weights, 'descend');
        relations = relations(indices_for_sorting);   % we sort the relations as well
        num_of_relations = length(relations);

        % assign to the vector
        old_to_new_indices(old_counter) = new_counter + 1;

        % update counter
        new_counter = new_counter + 1;   % anyway we need to increase the new counter by at least 1
        new_counter_temp = new_counter;
        current_max = v_weights(1);
        for i=1:num_of_relations

            if i ~= num_of_relations
                if v_weights(i+1)*approx_term >= current_max
                    % we want to combine this obs with the next obs
                else
                    % we don't combine this obs with the next one
                    current_max = v_weights(i+1);
                    new_counter = new_counter + 1;
                end
            else
                
            end
        end
        % assign to the vector
        new_to_old(new_counter_temp:new_counter) = old_counter;
    end
    % truncate the vector
    new_to_old = new_to_old(1:new_counter);
    
    % update new_obs_num
    new_obs_num = new_counter;
    
    % update approx term
    previous_approx_term = approx_term;
    approx_term = approx_term + approx_step;
    
    if mod(approx_iteration, 10) == 0 
        approx_step = approx_step * 2;
    end
    approx_iteration = approx_iteration + 1;
end
% in the last iteration we incremented it with no need, now we return it back 
approx_term = previous_approx_term;



% initializing the new "SDRP" matrix that will be used in the new problem,
% and the vector of observation weights
new_SDRP = false(new_obs_num , new_obs_num);   % zeros(new_obs_num , new_obs_num);
weights = nan(1, new_obs_num);


new_counter = 0;
max_ratio = 1;
for old_counter = 1:old_obs_num
    % find the relations and their ratios, and sort the ratios in descending order 
    relations = find(SDRP(old_counter , :));
    v_weights = one_minus_RATIO(old_counter , relations);
    [v_weights , indices_for_sorting] = sort(v_weights, 'descend');
    relations = relations(indices_for_sorting);   % we sort the relations as well
    num_of_relations = length(relations);
    
    new_counter = new_counter + 1;   % anyway we need to increase the new counter by at least 1
    current_max = v_weights(1);
    for i=1:num_of_relations
        new_matrix_relation = old_to_new_indices(relations(i));
        new_SDRP(new_counter, new_matrix_relation) = 1;
        if i ~= num_of_relations
            if v_weights(i+1)*approx_term >= current_max
                % we want to combine this obs with the next obs
            else
                % the made-up relation between duplicated observations
                new_SDRP(new_counter, new_counter + 1) = 1;
                weights(new_counter) = sqrt(current_max * v_weights(i)); % the weight we assign is a geometric mean of the minimum and the maximum. it also equals to: v_weights(i) * sqrt(current_max / v_weights(i));
                new_counter = new_counter + 1;
                % we don't combine this obs with the next one
                max_ratio = max(max_ratio , current_max / v_weights(i));
                current_max = v_weights(i+1);
            end
        else
            % in the last observation, we also update these
            max_ratio = max(max_ratio , current_max / v_weights(i));
            weights(new_counter) = sqrt(current_max * v_weights(i)); % the weight we assign is a geometric mean of the minimum and the maximum. it also equals to: v_weights(i) * sqrt(current_max / v_weights(i));
        end
    end
end

% if it is exact, it is 1. if it is not exact, it equals to sqrt(max_ratio)
% (which is no more than sqrt(approx_term), usually very close to sqrt(approx_term)) 
% we take the square root so that if the final result is x, then the lower
% bound is (x / Var_exact) and the upper bound is (x * Var_exact)
Var_approx_ratio = sqrt(max_ratio);


% % FLIP Part I: 
% % this is for efficiency of finding minimal cycles algorithm.
% % it is needed because the "find_all_minimal_cycles" function that will be
% % called from "HPZ_Houtman_Maks_Weighted_Cycles_Approach" starts from the
% % last vertex and goes backwords, and if old vertex A was splitted into A1,
% % A2,...,Ak such that A1->A2->...->Ak, we want that A1 will be handled
% % first by the "find_all_minimal_cycles" function, because then when we
% % continue to Ai = {A2,...,Ak} after removing A1, the algorithm will 
% % immediately recognize that there is no path back to Ai, since all paths 
% % leading to Ai go through A1.
% % if you use an implementation of "find_all_minimal_cycles" that starts
% % from the 1st vertex and goes forward, you should delete the "FLIP" Parts
% % I and II to regain efficiency.
% new_SDRP = new_SDRP(new_obs_num:(-1):1 , new_obs_num:(-1):1);
% weights = weights(new_obs_num:(-1):1);

% after transforming the problem, we can calculate the index:
if ~out_of_sample_required
    % without out-of-sample
    [~, ~, one_minus_v_new_obs, Var_out_of_sample_residuals, out_of_sample_one_minus_v] = HPZ_Weighted_Houtman_Maks (new_SDRP, weights, aggregators, main_folder_for_results);
else
    % with out-of-sample, we pass old_to_new_indices as varargin{1}
   [~, ~, one_minus_v_new_obs, Var_out_of_sample_residuals, out_of_sample_one_minus_v] = HPZ_Weighted_Houtman_Maks (new_SDRP, weights, aggregators, main_folder_for_results, old_to_new_indices); 
end

% % FLIP Part II:
% % RE-FLIP: turn everything back to the right order
% one_minus_v_new_obs = one_minus_v_new_obs(new_obs_num:(-1):1, :);


% transform 1-v to the original observations
one_minus_v = nan(old_obs_num, 3);
for i=1:new_obs_num
    one_minus_v(new_to_old(i), :) = max(one_minus_v(new_to_old(i), :), one_minus_v_new_obs(i, :));
end



% calculating the Varian Index for this subset
% this part is not necessary, as it is done again in HPZ_Varian_Manager for the full set of observations  
VARIAN = nan(1,3);
VARIAN(1) = sum(one_minus_v(:,1)) / old_obs_num;            % Mean
VARIAN(2) = sqrt(sum(one_minus_v(:,2).^2) / old_obs_num);   % AVG(SSQ)
VARIAN(3) = max(one_minus_v(:,3));                          % Max



% calculating in-sample residuals for this subset
% this part is not necessary, as it is done again in HPZ_Varian_Manager for the full set of observations  

% calculate in-sample difference residuals
average_var_residuals = HPZ_Consistency_Indices_In_Sample_Difference_Residuals_Calc (one_minus_v(:,1), @mean);
meanssq_var_residuals = HPZ_Consistency_Indices_In_Sample_Difference_Residuals_Calc (one_minus_v(:,2), @(x) sqrt(meansqr(x)));
max_var_residuals = HPZ_Consistency_Indices_In_Sample_Difference_Residuals_Calc (one_minus_v(:,3), @max); 

% summarize all residuals;
% the first 3 are in-sample component residuals, 
% and the second 3 are in-sample difference residuals.
Var_in_sample_residuals = [one_minus_v(:,1) , one_minus_v(:,2) , one_minus_v(:,3) , average_var_residuals , meanssq_var_residuals , max_var_residuals];



end
