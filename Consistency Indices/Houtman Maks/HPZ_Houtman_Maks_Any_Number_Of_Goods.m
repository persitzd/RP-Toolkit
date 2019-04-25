function [HM, HM_raw, HM_raw_residuals, HM_raw_pseudo_residuals] = HPZ_Houtman_Maks_Any_Number_Of_Goods (DRP, SDRP, identical_choice, varargin)

% for this specific algorithm (finding minimal cycles),
% we want DRP to not contain 1 in (i,j) if i and j are identical choices 
DRP = DRP .* (~identical_choice); % (added 20.02.19)

if isempty(varargin)
    % in the first round, we still haven't created the observations-cycles
    % matrix, so we need to create it. in the following, recursive rounds, 
    % the matrix will be already passed down

    % we always start with approach 1
    approach = 1;
    
    % number of observations
    [obs_num , ~] = size(DRP);
    
    % when we use an algorithm to find all minimal cycles, we don't want that 
    % there would be edges from a vertex to itself, so we make sure it won't occur
    for i=1:obs_num
        DRP(i,i) = 0;
    end
    
    % find all cycles that don't contain sub-cycles
    % ----------------------------------------
    
    % --- according to SARP ---
    % NOTE! When using this function to calculate HM according to SARP, the
    % given DRP must NOT include relation between identical choices, that
    % is: we use DRP & ~identical_choice instead of DRP
    %cycles = find_all_minimal_cycles (DRP & ~identical_choice); 
    %%cycles = find_cycles_without_subcycles (DRP & ~identical_choice);   
    
    % --- according to GARP ---
    cycles = find_all_minimal_GARP_cycles (DRP, SDRP);
    %cycles = find_GARP_cycles_without_GARP_subcycles (DRP, SDRP);
    
    % ----------------------------------------
    
    % number of cycles
    cycles_num = length(cycles);
    
    % initialize the matrix
    obs_cycle_mat = false(obs_num , cycles_num); % = zeros(obs_num , cycles_num);
    
    % we want that in the (i,j) location in the matrix there will be 1,
    % if and only if observation i is part of cycle j
    for j=1:cycles_num
        current_cycle = cycles{j};
        for k=1:length(current_cycle)   %k=1:(length(current_cycle)-1)
            obs_cycle_mat(current_cycle(k) , j) = true; % = 1
        end
    end
    
else
    % in the recursive iterations, we already have obs_cycle_mat, and we
    % don't need DRP and SDRP anymore
    obs_cycle_mat = varargin{1};
    approach = varargin{2};
    
    % number of observations and of cycles
    [obs_num , cycles_num] = size(obs_cycle_mat);
    
end



% if there are no observation, we don't need to perform many calculations
if obs_num == 0
    HM = 0;
    HM_raw = 0;
    HM_raw_residuals = [];
    HM_raw_pseudo_residuals = [];
    return
% if there is only 1 observation, we don't need to perform many calculations
elseif obs_num == 1
    if any(obs_cycle_mat(1,:)) % sum(obs_cycle_mat(1,:)) > 0
        % we need to choose/drop this observation in order to take care of
        % the remaining cycle/s
        HM = 1;
        HM_raw = 1;
        HM_raw_residuals = 1;
        HM_raw_pseudo_residuals = 1;
    else
        % we don't need to choose/drop this observation
        HM = 0;
        HM_raw = 0;
        HM_raw_residuals = 0;
        HM_raw_pseudo_residuals = 0;
    end
    return
end





% initializations (just in case)
HM = 0; %#ok<NASGU>
HM_raw = 0;
HM_raw_residuals = zeros(obs_num, 1);
HM_raw_pseudo_residuals = zeros(obs_num, 1);

% in the last part of this function, there is a case where we "split up"
% and check 2 different scenarios. it could be that we got a raw index of
% say 0.07 in the first, so in the second we know that there is no point in
% having more than 0.07, and this is what this variable is for.
if length(varargin) < 3
    upper_bound = inf;
else
    upper_bound = varargin{3};
end





% sums the matrix per column
% for each cycle, it finds how many observations (from those remaining) are
% involved / part of this cycle
obs_involved_per_cycle = sum(obs_cycle_mat);

% we want to find all cycles that have only 1 observation
% we MUST take these observations in order to "cover" these cycles 
cycles_with_one_obs = find(obs_involved_per_cycle == 1);

if ~isempty(cycles_with_one_obs)

    % we need to find the observations that belong to these cycles
    obs_that_must_be_chosen = zeros(length(cycles_with_one_obs), 1);
    for k=1:length(cycles_with_one_obs)
        obs_that_must_be_chosen(k) = find(obs_cycle_mat(:, cycles_with_one_obs(k)));
    end
    % delete duplicates
    obs_that_must_be_chosen = unique(obs_that_must_be_chosen);

    % increase the HM_raw, since we "choose" (=drop) these observations
    HM_raw = HM_raw + length(obs_that_must_be_chosen);
    
    if upper_bound < HM_raw
        % no need to keep going any further
        HM_raw = inf;
        HM = HM_raw;
        HM_raw_residuals = inf(obs_num, 1);
        HM_raw_pseudo_residuals = inf(obs_num, 1);
        return
    else
        % update the upper bound; we will pass this new upper bound to the
        % recursive call ahead 
        upper_bound = upper_bound - HM_raw;
    end
    
    % initialization
    cycles_involved = [];
    for i=1:length(obs_that_must_be_chosen)
        % for each of these observations, we need to delete the cycles that it
        % is involved with
        obs_cycles = find(obs_cycle_mat(obs_that_must_be_chosen(i), :));
        cycles_involved = unique([cycles_involved , obs_cycles]);
    end

    % finally, we can delete these observations and these cycles from the 
    % obs_cycle_mat matrix
    cycles_to_delete = zeros(cycles_num, 1);
    cycles_to_delete(cycles_involved) = 1;
    cycles_to_remain = (cycles_to_delete == 0);
    obs_to_delete = zeros(obs_num, 1);
    obs_to_delete(obs_that_must_be_chosen) = 1;
    obs_to_delete = (obs_to_delete ~= 0);
    obs_to_remain = (obs_to_delete == 0);
    obs_cycle_mat = obs_cycle_mat(obs_to_remain , cycles_to_remain);

    % update residuals
    HM_raw_residuals(obs_to_delete) = 1;
    HM_raw_pseudo_residuals(obs_to_delete) = 1;
    
    if isempty(obs_cycle_mat)
        % then we finished
    else
        % we need to use another recursive call
        [~, remaining_HM_raw, remaining_HM_raw_residuals, remaining_HM_raw_pseudo_residuals] = HPZ_Houtman_Maks_Any_Number_Of_Goods ([], [], [], obs_cycle_mat, 1, upper_bound);
        HM_raw = HM_raw + remaining_HM_raw;
        HM_raw_residuals(obs_to_remain) = remaining_HM_raw_residuals;
        HM_raw_pseudo_residuals(obs_to_remain) = remaining_HM_raw_pseudo_residuals;
    end
    
else
    % if isempty(cycles_with_one_obs),
    % if we can't tell about any observation that we "must" choose it,
    % then we need to take an observation (e.g. the observation with the
    % biggest number of cycles, or some aggregate of cycles), and check
    % both options: the option of choosing (dropping) it, and the option
    % of not choosing it.
    
    % we apply 1 of 2 algorithms here:
    % Algorithm 1 : delete every observation i such that there is another
    % observation j such that the set of cycles of i is a subset of the set
    % of cycles of j. if the sets are equal, we will arbitrarily delete
    % only one of them. we only want to find some minimal set of
    % observations to drop, not every set, therefore we can always assume
    % that it would be better to drop the observation that is involved with
    % the same cycles and possibly with other cycles.
    % Algorithm 2 : choose an observation, then check 2 options: either to
    % take (drop) this observation or not (each option is checked by a
    % recursive call), then we determine which of the 2 options is better
    % (has lower Houtman-Maks raw index).
    
    if approach == 1
        
        % --- Algorithm 1 ---
        % to improve time complexity, we first create a cell array, in the i'th
        % place it contains a vector with all the indices of cycles that
        % observation i is involved with
        cycles_of_each_obs = cell(obs_num, 1);
        for i=1:obs_num
            cycles_of_each_obs{i} = find(obs_cycle_mat(i,:));
        end
        % to improve time complexity, during the process we keep track over 
        % observations that we already decided to drop 
        % (NOTE: "drop" here means that we DON'T choose to drop this
        % observation, it means that we just "drop" it from our matrix)
        was_dropped = zeros(obs_num, 1);
        % for each couple of observations, we check if one is a subset of the other 
        for i=1:(obs_num-1) 
            for t=(i+1):obs_num
                % if we already dropped one of them, we don't need to check
                if was_dropped(i) || was_dropped(t)
                    continue
                end
                cycles_i = cycles_of_each_obs{i};
                cycles_t = cycles_of_each_obs{t};
                i_not_subset = false;
                t_not_subset = false;
                l_cycle_i = length(cycles_i);
                l_cycle_t = length(cycles_t);
                if l_cycle_i > l_cycle_t
                    i_not_subset = true;
                    difference = l_cycle_i - l_cycle_t;
                    ai = 1;
                    for at=1:l_cycle_t
                        while ai <= l_cycle_i && cycles_t(at) ~= cycles_i(ai)
                            ai = ai + 1;
                        end
                        if ai-at > difference
                            t_not_subset = true;
                            break
                        end
                    end
                elseif l_cycle_i < l_cycle_t
                    t_not_subset = true;
                    difference = l_cycle_t - l_cycle_i;
                    at = 1;
                    for ai=1:l_cycle_i
                        while at <= l_cycle_t && cycles_t(at) ~= cycles_i(ai)
                            at = at + 1;
                        end
                        if at-ai > difference
                            i_not_subset = true;
                            break
                        end
                    end
                else % cycles are of equal length - either identical (both are subsets of each other) or not  
                    % (ahead we use "if... elseif..." so even if both are
                    % false, only one will be dropped)
                    for a=1:l_cycle_i
                        if cycles_t(a) ~= cycles_i(a)
                            i_not_subset = true;
                            t_not_subset = true;
                            break
                        end
                    end
                end

                % if one of them is a subset - we delete it
                % (it must be "if... elseif..." and not "if... if...", so we will never delete both) 
                if ~i_not_subset
                    was_dropped(i) = 1;
                elseif ~t_not_subset
                    was_dropped(t) = 1;
                end
            end
        end

        if ~any(was_dropped)   % sum(was_dropped) == 0
            % we couldn't "drop" any observation
            % we make a recursive call to the 2nd approach (Algorithm 2)
            if ~isempty(obs_cycle_mat)
                [~, HM_raw, HM_raw_residuals, HM_raw_pseudo_residuals] = HPZ_Houtman_Maks_Any_Number_Of_Goods ([], [], [], obs_cycle_mat, 2, upper_bound);
            end
        else
            remaining_obs = find(~was_dropped);   % find(was_dropped == 0);
            if ~isempty(obs_cycle_mat) && ~isempty(remaining_obs)
                % we need to use another recursive call to this approach (1st approach) 
                [~, remaining_HM_raw, remaining_HM_raw_residuals, remaining_HM_raw_pseudo_residuals] = HPZ_Houtman_Maks_Any_Number_Of_Goods ([], [], [], obs_cycle_mat(remaining_obs, :), 1, upper_bound);
                HM_raw = HM_raw + remaining_HM_raw;
                HM_raw_residuals(remaining_obs) = remaining_HM_raw_residuals;
                HM_raw_pseudo_residuals(remaining_obs) = remaining_HM_raw_pseudo_residuals;
            end
        end
    
    elseif approach == 2
        
        % initializations
        assigned_obs = false(1, size(obs_cycle_mat,1));
        obs_subsets = cell(0,0);
        
        while ~all(assigned_obs)
            
            [~, current_obs] = max(~assigned_obs);
            
            % initialization of independent subset of obs
            obs_current_subset = false(1, size(obs_cycle_mat,1));
            previous_obs_current_subset = obs_current_subset;
            obs_current_subset(current_obs) = true;
            
            while any(obs_current_subset & ~previous_obs_current_subset)   % any(obs_to_assign)   % sum(obs_to_assign) > 0
                
                previous_obs_current_subset = obs_current_subset;
                
                % all cycles (columns) that are involved with the obs of 
                % the current version of independent subset of obs 
                cycles_involved = any(obs_cycle_mat(obs_current_subset, :), 1);
                
                % updating the independent subset of obs
                obs_current_subset = any(obs_cycle_mat(: , cycles_involved), 2)';
            end
            
            obs_subsets{end+1} = obs_current_subset; %#ok<AGROW>
            
            % update those who were assigned
            assigned_obs = assigned_obs | obs_current_subset;
        end
        
        % for each subset - calculate independently
        for subset_num = 1:length(obs_subsets)
            
            obs_of_this_subset = obs_subsets{subset_num};
            cycles_of_this_subset = any(obs_cycle_mat(obs_of_this_subset, :), 1);
            
            [~, subset_HM_raw, subset_HM_raw_residuals, subset_HM_raw_pseudo_residuals] = HPZ_Houtman_Maks_Any_Number_Of_Goods ([], [], [], obs_cycle_mat(obs_of_this_subset, cycles_of_this_subset), 3, upper_bound);
            
            HM_raw = HM_raw + subset_HM_raw;
            HM_raw_residuals(obs_of_this_subset, :) = subset_HM_raw_residuals;
            HM_raw_pseudo_residuals(obs_of_this_subset, :) = subset_HM_raw_pseudo_residuals;
            
            % update upper bound
            upper_bound = upper_bound - subset_HM_raw;
            % we do not update upper_bound(1) because here (unlike approach 3) 
            % we don't compare alternatives, but we just calculate
            % separately for different subsets.
            if upper_bound < 0
               % no need to keep going any further
                HM_raw = inf;
                HM = HM_raw;
                HM_raw_residuals = inf(obs_num, 1);
                HM_raw_pseudo_residuals = inf(obs_num, 1);
                return
            end
        end
        
    else % approach == 3
    
        % --- Algorithm 2 ---
        % we prefer to take the observation that is involved with the biggest
        % amonut of cycles, but involvment in a small cycle is of more "value"
        % than involvement in a big cycle, so we give each cycle a weight that
        % equals to 1 / (cycle length)
        cycles_values = 1 ./ obs_involved_per_cycle;
        observations_values = zeros(obs_num, 1);
        for i=1:obs_num
    %         obs_cycles = find(obs_cycle_mat(i, :));
    %         observations_values(i) = sum(cycles_weights(obs_cycles));
            observations_values(i) = sum(obs_cycle_mat(i, :) .* cycles_values);
        end
        [~, obs_to_take] = max(observations_values);


        % option I - choose (drop) this observation
        HM_raw_1 = HM_raw + 1;
        upper_bound_1 = upper_bound - 1;
        if upper_bound_1 < 0
            % no need to keep going any further
            HM_raw_1 = inf;
            remaining_HM_raw_residuals_1 = inf(obs_num, 1);
            remaining_HM_raw_pseudo_residuals_1 = inf(obs_num, 1);
            obs_to_remain_1 = [];
            obs_cycle_mat_1 = [];
        else
            remaining_HM_raw_residuals_1 = [];
            remaining_HM_raw_pseudo_residuals_1 = [];
            cycles_to_remain_1 = ~obs_cycle_mat(obs_to_take, :); % (obs_cycle_mat(obs_to_take, :) == 0);
            obs_to_remain_1 = [1:(obs_to_take-1) , (obs_to_take+1):obs_num];
            obs_cycle_mat_1 = obs_cycle_mat(obs_to_remain_1 , cycles_to_remain_1);
            if ~isempty(obs_cycle_mat_1)
                % we need to use another recursive call
                [~, remaining_HM_raw_1, remaining_HM_raw_residuals_1, remaining_HM_raw_pseudo_residuals_1] = HPZ_Houtman_Maks_Any_Number_Of_Goods ([], [], [], obs_cycle_mat_1, 1, upper_bound_1);
                HM_raw_1 = HM_raw_1 + remaining_HM_raw_1;
            end
        end


        % option II - don't choose (don't drop) this observation
        HM_raw_2 = HM_raw;
        upper_bound_2 = min(HM_raw_1 , upper_bound);
        if upper_bound_2 < 0
            % no need to keep going any further
            HM_raw_2 = inf;
            remaining_HM_raw_residuals_2 = inf(obs_num, 1);
            remaining_HM_raw_pseudo_residuals_2 = inf(obs_num, 1);
            obs_to_remain_2 = [];
            obs_cycle_mat_2 = [];
        else
            remaining_HM_raw_residuals_2 = [];
            remaining_HM_raw_pseudo_residuals_2 = [];
            obs_to_remain_2 = [1:(obs_to_take-1) , (obs_to_take+1):obs_num];
            obs_cycle_mat_2 = obs_cycle_mat(obs_to_remain_2 , :);
            if ~isempty(obs_cycle_mat_2)
                % we need to use another recursive call
                [~, remaining_HM_raw_2, remaining_HM_raw_residuals_2, remaining_HM_raw_pseudo_residuals_2] = HPZ_Houtman_Maks_Any_Number_Of_Goods ([], [], [], obs_cycle_mat_2, 1, upper_bound_2);
                HM_raw_2 = HM_raw_2 + remaining_HM_raw_2;
            end
        end

        if HM_raw_1 < HM_raw_2
            % I (choose / drop) is better
            HM_raw = HM_raw_1;
            HM_raw_residuals(obs_to_take) = 1;
            HM_raw_pseudo_residuals(obs_to_take) = 1;
            if ~isempty(obs_cycle_mat_1)
                HM_raw_residuals(obs_to_remain_1) = remaining_HM_raw_residuals_1;
                HM_raw_pseudo_residuals(obs_to_remain_1) = remaining_HM_raw_pseudo_residuals_1;
            end

        elseif HM_raw_1 > HM_raw_2
            % II (don't choose / don't drop) is better
            HM_raw = HM_raw_2;
            if ~isempty(obs_cycle_mat_2)
                HM_raw_residuals(obs_to_remain_2) = remaining_HM_raw_residuals_2;
                HM_raw_pseudo_residuals(obs_to_remain_2) = remaining_HM_raw_pseudo_residuals_2;
            end

        else % they are equal
            HM_raw = HM_raw_1;  % = HM_raw_2
            % we arbitrarily take the 1st alternative for raw (not pseudo) residuals:  
            HM_raw_residuals(obs_to_take) = 1;
            if ~isempty(obs_cycle_mat_1)
                HM_raw_pseudo_residuals(obs_to_remain_1) = remaining_HM_raw_pseudo_residuals_1;
            end
            % for raw pseudo residuals, we take any vertex that is in any of both  
            HM_raw_pseudo_residuals(obs_to_take) = 1;
            if ~isempty(obs_cycle_mat_1)
                HM_raw_pseudo_residuals(obs_to_remain_1) = max(remaining_HM_raw_pseudo_residuals_1, HM_raw_pseudo_residuals(obs_to_remain_1));
            end
            if ~isempty(obs_cycle_mat_2)
                HM_raw_pseudo_residuals(obs_to_remain_2) = max(remaining_HM_raw_pseudo_residuals_2, HM_raw_pseudo_residuals(obs_to_remain_2));
            end
        end
        
    end
    
end



% calculate the final HM index by dividing by the number of observations
HM = HM_raw / obs_num;



end