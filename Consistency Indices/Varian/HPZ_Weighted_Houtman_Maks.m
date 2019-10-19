function [w_HM, w_HM_raw, w_HM_weight_per_observation, Varian_out_of_sample_residuals, out_of_sample_one_minus_v] = HPZ_Weighted_Houtman_Maks (SDRP, weights, aggregators, varargin)


% This function implements the same algorithm that is implemented in 
% "HPZ_Houtman_Maks_Cycles_Approach", but with some modifications, since
% this function can handle the case that dropping observation i has a
% weight of w_i, instead of being equal and constant (weight of 1) for all
% observations. 
% Since the weights are different, there is meaning in using different
% aggregates, so here we also check 3 different aggregates:
% maximum, average, and average sum of squares (AVGSSQ).
%
% The purpose of this function is currently solely as a mean to calculate
% the Varian Efficiency Index; we proved that any problem of calculating
% the Varian index can be translated to a problem of calculating the HM
% index with weights, just like this function does.



% number of observations
[obs_num , ~] = size(SDRP);

% initializations of outputs
w_HM = [nan, nan, nan];
w_HM_raw = [nan, nan, nan];
w_HM_weight_per_observation = nan(obs_num, 3);



% when we use Jhonson algorithm, we don't want that there would be
% edges from a vertex to itself, so we make sure it won't occur
for i=1:obs_num
    SDRP(i,i) = 0;
end

% takes less space, and more efficient
SDRP = logical(sparse(SDRP));



% FLIP Part I: 
% this is for efficiency of finding minimal cycles algorithm.
% it is needed because the "find_all_minimal_cycles" function that will be
% called from "HPZ_Houtman_Maks_Weighted_Cycles_Approach" starts from the
% last vertex and goes backwords, and if old vertex A was splitted into A1,
% A2,...,Ak such that A1->A2->...->Ak, we want that A1 will be handled
% first by the "find_all_minimal_cycles" function, because then when we
% continue to Ai = {A2,...,Ak} after removing A1, the algorithm will 
% immediately recognize that there is no path back to Ai, since all paths 
% leading to Ai go through A1.
% if you use an implementation of "find_all_minimal_cycles" that starts
% from the 1st vertex and goes forward, you should delete the "FLIP" Parts
% I and II to regain efficiency.
flipped_SDRP = SDRP(obs_num:(-1):1 , obs_num:(-1):1);

% start_time = now;

% find all minimal cycles (cycles that don't contain sub-cycles)
cycles = find_all_minimal_cycles (flipped_SDRP);
cycles_num = length(cycles);

% % USE THIS CODE TO DETERMINE THE RUNNING TIME OF FINDING ALL MINIMAL
% % CYCLES (in comparison to total running time)
% end_time = now;
% running_time = datevec(end_time - start_time);
% months = running_time(2);
% days = running_time(3);
% hours = running_time(4);
% minutes = running_time(5);
% seconds = running_time(6);
% fprintf('\ncycles running time was:\n%d months, %d days, %d hours, %d minutes and %.3f seconds.\n',...
%                                                                 months, days, hours, minutes, seconds);

% FLIP Part II:
% RE-FLIP: turn everything back to the right order
unflip = obs_num:(-1):1;
for c = 1:cycles_num
    cycles{c} = unflip(cycles{c});
end



% now we make preparations for the use of integer programming (using intlinprog)   

if aggregators(1) || aggregators(2)
    
    % in both mean and mean(SSQ) aggregators we use integer programming to perform the calculation  
    
    % weight of each observations
    f = weights';

    % all values are integers that receive value of either 0 or 1
    intcon = 1:obs_num;
    lb = zeros(1, obs_num);
    ub = ones(1, obs_num);

    % the Ax<=b restriction will make sure that each cycle is covered by at least one observation 
    A = zeros(cycles_num, obs_num);  % initialization
    b = - ones(cycles_num, 1);
    for c = 1:cycles_num
        currect_cycle = cycles{c};
        A(c, currect_cycle) = -1;
    end

    % avoid displaying algorithm stages in command window
    intlinprog_options = optimoptions('intlinprog', 'Display','off');
    
    % intlinprog function was changed between R2017a and R2017b,
    % hence we need to address that. we don't bother addressing version before 2014a, 
    % since intlinprog wasn't even introduced before then (there was bintprog instead).  
    is_old_version_intlinprog = all(version('-release') == '2017a') || ...
                    all(version('-release') == '2016b') || all(version('-release') == '2016a') || ...
                    all(version('-release') == '2015b') || all(version('-release') == '2015a') || ...
                    all(version('-release') == '2014b') || all(version('-release') == '2014a');
end


if aggregators(1)   % mean aggregator
    
    % calculate weighted houtman-maks for sum aggregator
    if is_old_version_intlinprog
        [x, fval, exitflag, ~] = intlinprog(f, intcon, A,b, [],[], lb,ub, intlinprog_options);
    else
        [x, fval, exitflag, ~] = intlinprog(f, intcon, A,b, [],[], lb,ub, [], intlinprog_options);
    end

    if exitflag == 1
        % assigning to the matrices
        w_HM_raw(1) = fval;   % sum aggregator
        w_HM_weight_per_observation(:,1) = x .* f;
    else
        % assigning nan values and ending the function
        w_HM_raw(1) = nan;
        w_HM_weight_per_observation(:,1) = nan(1, obs_num);
    end
    
end


if aggregators(2)   % mean(SSQ) aggregator
    
    % calculate weighted houtman-maks for sum-of-squares aggregator (by using f^2 instead of f)  
    if is_old_version_intlinprog
        [x, fval, exitflag, ~] = intlinprog(f.^2, intcon, A,b, [],[], lb,ub, intlinprog_options);
    else
        [x, fval, exitflag, ~] = intlinprog(f.^2, intcon, A,b, [],[], lb,ub, [], intlinprog_options);
    end

    if exitflag == 1
        % assigning to the matrices
        w_HM_raw(2) = fval;   % sum-of-squares aggregator
        w_HM_weight_per_observation(:,2) = x .* f;
    else
        % assigning nan values and ending the function
        w_HM_raw(2) = nan;
        w_HM_weight_per_observation(:,2) = nan(1, obs_num);
    end
    
end


if aggregators(3)   % max aggregator
    
    % here we cannot use integer programming, since the aggregator is not additive  
    
    % initialize the matrix
    obs_cycle_mat = false(obs_num , cycles_num); % zeros(obs_num , cycles_num);

    % we want that in the (i,j) location in the matrix there will be 1,
    % if and only if observation i is part of cycle j
    for j=1:cycles_num
        current_cycle = cycles{j};
        obs_cycle_mat(current_cycle , j) = true;
    end

    % takes less space, and more efficient
    obs_cycle_mat = sparse(logical(obs_cycle_mat));
    
    % Call to Weighted Houtman-Maks recursive algorithm
    [~, w_HM_max_raw, w_HM_max_weight_per_observation] = HPZ_Weighted_Houtman_Maks_Max_Aggregator (obs_cycle_mat, weights);
    w_HM_raw(3) = w_HM_max_raw;
    w_HM_weight_per_observation(:,3) = w_HM_max_weight_per_observation;
end



% calculate the final HM index by dividing by the number of observations (if needed)  
w_HM(1) = w_HM_raw(1) / obs_num;   % mean aggregator
w_HM(2) = sqrt(w_HM_raw(2) / obs_num);   % mean(SSQ) aggregator
w_HM(3) = w_HM_raw(3);   % max aggregator





% out-of-sample (if required)

if ~isempty(varargin)
    
    % this vector will allow us to remove at once all new observations that
    % originated from a specific old/original observation (in the Varian problem)  
    old_to_new_indices = varargin{1};
    % old/original number of observations in the Varian problem
    old_obs_num = length(old_to_new_indices);
    % new obs num
    new_obs_num = obs_num;
    
    % initializations of Varian out-of-sample vectors
    average_var_out_sample = nan(old_obs_num, 1);
    meanssq_var_out_sample = nan(old_obs_num, 1);
    max_var_out_sample = nan(old_obs_num, 1);
    
    % each "row" of this matrix will contain 3 x (old_obs_num-1) matrix,  
    % representing the 1-v values of the out-of-sample solution
    out_of_sample_one_minus_v = nan(old_obs_num, old_obs_num-1, 3);
    
    % calculate out-of-sample for each old observation
    for i=1:old_obs_num
        
%         % initialization of the i'th instance of this cell array
%         out_of_sample_one_minus_v{i} = nan(old_obs_num-1, 3);
        
        % all new observations that originated from old observation i
        if i < old_obs_num
            new_obs_from_old_obs_i = old_to_new_indices(i):(old_to_new_indices(i+1)-1);
        else
            new_obs_from_old_obs_i = old_to_new_indices(i):new_obs_num;
        end
        
        % all new observations except for those that originated from old observation i  
        all_except_i = true(1, new_obs_num);
        all_except_i(new_obs_from_old_obs_i) = false;
        
        % by construction, the first new observation must be included in  
        % all the cycles that the old observations was included in.
        cycles_that_include_observation_i = logical(A(:,old_to_new_indices(i)));
        
        if ~any(~cycles_that_include_observation_i)
            
            % if all cycles include old observation i,
            % then all cycles were "solved" - no need to do anything more
            out_of_sample_one_minus_v(i,:,1) =  zeros(1, old_obs_num-1);
            out_of_sample_one_minus_v(i,:,2) =  zeros(1, old_obs_num-1);
            out_of_sample_one_minus_v(i,:,3) =  zeros(1, old_obs_num-1);
            average_var_out_sample(i) = 0;
            meanssq_var_out_sample(i) = 0;
            max_var_out_sample(i) = 0;
            
        else
            
            % if there are cycles not involved in old observation i
            
            % this vector will translate the new observations without i 
            % to the old observations without i
            old_to_new_indices_without_i = [old_to_new_indices(1:(i-1)) , old_to_new_indices((i+1):end) - length(new_obs_from_old_obs_i)];
            
            if aggregators(1) || aggregators(2)
                % we need to reduce the sizes of these vectors
                new_obs_num_without_i = new_obs_num - length(new_obs_from_old_obs_i);
                intcon_residuals = 1:new_obs_num_without_i;
                lb_residuals = lb(1:new_obs_num_without_i);
                ub_residuals = ub(1:new_obs_num_without_i);

                % a matrix A without observation i and without all minimal
                % cycles that include observation i
                A_without_obs_i = A(~cycles_that_include_observation_i , all_except_i);
                b_without_obs_i = b(~cycles_that_include_observation_i);
                f_without_obs_i = f(all_except_i);
            end

            % Mean
            if aggregators(1)
                % find the raw Varian Index without this observations 
                [x, fval, exitflag, ~] = intlinprog(f_without_obs_i, intcon_residuals, A_without_obs_i,b_without_obs_i, [],[], lb_residuals,ub_residuals, [], intlinprog_options);
                if exitflag == 1
                    % transform 1-v to the original observations
                    one_minus_v = nan(old_obs_num-1, 1);
                    for j=1:(old_obs_num-1)
                        if j < old_obs_num-1
                            current_new_obs = old_to_new_indices_without_i(j) : (old_to_new_indices_without_i(j+1)-1);
                        else
                            current_new_obs = old_to_new_indices_without_i(j) : sum(all_except_i);
                        end
                        x_times_f = x .* f_without_obs_i;
                        one_minus_v(j, :) = max(x_times_f(current_new_obs));
                    end
                    % assigning to the vector and cell array
                    out_of_sample_one_minus_v(i,:,1) = one_minus_v; % x .* f_without_obs_i;
                    average_var_out_sample(i) = mean(one_minus_v); % mean(x .* f_without_obs_i);   % = fval / (old_obs_num-1)
                end
            end

            % AVG(SSQ)
            if aggregators(2)
                % find the raw Varian Index without this observations 
                [x, fval, exitflag, ~] = intlinprog(f_without_obs_i.^2, intcon_residuals, A_without_obs_i,b_without_obs_i, [],[], lb_residuals,ub_residuals, [], intlinprog_options);
                if exitflag == 1
                    % transform 1-v to the original observations
                    one_minus_v = nan(old_obs_num-1, 1);
                    for j=1:(old_obs_num-1)
                        if j < old_obs_num-1
                            current_new_obs = old_to_new_indices_without_i(j) : (old_to_new_indices_without_i(j+1)-1);
                        else
                            current_new_obs = old_to_new_indices_without_i(j) : sum(all_except_i);
                        end
                        x_times_f = x .* f_without_obs_i;
                        one_minus_v(j, :) = max(x_times_f(current_new_obs));
                    end
                    % assigning to the vector and cell array
                    out_of_sample_one_minus_v(i,:,2) = one_minus_v; % x .* f_without_obs_i;
                    meanssq_var_out_sample(i) = sqrt(mean(one_minus_v.^2)); % sqrt(mean(x .* f_without_obs_i.^2));   % = sqrt(fval / (old_obs_num-1))
                end
            end

            % Max
            if aggregators(3)
                % truncating the matrix and vector for all observations except
                % those that originated in old observation i, and to all cycles
                % except for those that included old observation i
                obs_cycle_mat_without_obs_i = obs_cycle_mat(all_except_i , cycles_that_include_observation_i);
                weights_without_obs_i = weights(all_except_i);

                % Call to Weighted Houtman-Maks recursive algorithm
                [~, w_HM_max_raw, w_HM_max_weight_per_observation] = HPZ_Weighted_Houtman_Maks_Max_Aggregator (obs_cycle_mat_without_obs_i, weights_without_obs_i);
                
                % transform 1-v to the original observations
                one_minus_v = nan(old_obs_num-1, 1);
                for j=1:(old_obs_num-1)
                    if j < old_obs_num-1
                        current_new_obs = old_to_new_indices_without_i(j) : (old_to_new_indices_without_i(j+1)-1);
                    else
                        current_new_obs = old_to_new_indices_without_i(j) : sum(all_except_i);
                    end
                    one_minus_v(j, :) = max(w_HM_max_weight_per_observation(current_new_obs));
                end
                    
                % assigning to the vector and cell array
                out_of_sample_one_minus_v(i,:,3) = one_minus_v;   % w_HM_max_weight_per_observation;
                max_var_out_sample(i) = max(one_minus_v);   % max(w_HM_max_weight_per_observation);   % = w_HM_max_raw
            end
            
        end
        
    end   % end of loop over observations
    
    % finalizing the out-of-sample residuals
    Varian_out_of_sample_residuals = [average_var_out_sample , meanssq_var_out_sample , max_var_out_sample];
    
else
    
    % initializations to avoid errors
    Varian_out_of_sample_residuals = nan(1, obs_num);
    out_of_sample_one_minus_v = [];
    
end






% % if there are no observations, we don't need to perform many calculations
% if obs_num == 0
%     w_HM_raw = [0,0,0];
%     w_HM = w_HM_raw;
%     w_HM_weight_per_observation = [];
%     return
% % if there is only 1 observation, we don't need to perform many calculations
% elseif obs_num == 1
%     if any(obs_cycle_mat(1,:)) % sum(obs_cycle_mat(1,:)) > 0
%         % we need to choose/drop this observation in order to take care of
%         % the remaining cycle/s
%         w_HM_raw = [weights(1) , weights(1) , weights(1)^2];
%         w_HM = w_HM_raw;
%         w_HM_weight_per_observation = [weights(1) , weights(1) , weights(1)];
%     else
%         % we don't need to choose/drop this observation
%         w_HM_raw = [0,0,0];
%         w_HM = w_HM_raw;
%         w_HM_weight_per_observation = w_HM_raw;
%     end
%     return
% end



end