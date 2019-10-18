function [Varian_Bounds, Var_in_sample_component_residuals, Var_out_of_sample_residuals, one_minus_v, Varian_Mat] = HPZ_Varian_Manager (VARIAN_flags, expenditure, identical_choice, index_threshold, SDRP, Varian_algorithm_settings)   % , residuals_waitbar, current_run, total_runs, subject_ID



% number of observations
[obs_num , ~] = size(expenditure);


% whether out-of-sample is required by the user
out_of_sample_required = (VARIAN_flags(2) && VARIAN_flags(4));


% we calculate SDRP with transitivity, instead of DRP with transitivity,
% since in the case of Varian weak relations are of no interest
[~, transitive_SDRP, ~] = GARP_based_on_DRP_and_SDRP(SDRP, SDRP);

% divide the observations to smallest possible subsets such that
% observations from different subsets never belong to one cycle
% (a cycle of strictly revealed relations in this case)
[subsets_of_observations, num_of_subsets] = HPZ_Indices_Problem_Distribute(transitive_SDRP, transitive_SDRP);

% initialization (it must be initialized to 0)
one_minus_v_avg = nan(obs_num, 3);
one_minus_v_meanssq = nan(obs_num, 3);
one_minus_v_max = nan(obs_num, 3);

% this cell array will contain the out-of-sample one-minus-v for each subset   
out_of_sample_one_minus_v_per_subset = cell(1, num_of_subsets);

% now we calculate the Varian separately for each subset
for i=1:num_of_subsets
    
    current_subset = subsets_of_observations{i};
    subset_obs_num = length(current_subset);
    
    if subset_obs_num > 1
        
        % extract the relevant expenditures and identical choices
        subset_expenditure = expenditure(current_subset , current_subset);
        subset_identical_choice = identical_choice(current_subset , current_subset);
        
        % calculate the Varian index for this subset of choices
        [~, subset_Varian_approx_ratio, ~, subset_one_minus_v, ~, subset_out_of_sample_one_minus_v] = HPZ_Varian_index_based_on_Houtman_Maks (subset_expenditure, subset_identical_choice, index_threshold, Varian_algorithm_settings, VARIAN_flags(5:7), out_of_sample_required);
        
        % update lower bound, approximate/exact and upper bound of 1-v
        one_minus_v_avg(current_subset, 1:3)     = [subset_one_minus_v(:,1)/subset_Varian_approx_ratio , subset_one_minus_v(:,1) , subset_one_minus_v(:,1)*subset_Varian_approx_ratio];
        one_minus_v_meanssq(current_subset, 1:3) = [subset_one_minus_v(:,2)/subset_Varian_approx_ratio , subset_one_minus_v(:,2) , subset_one_minus_v(:,2)*subset_Varian_approx_ratio];
        one_minus_v_max(current_subset, 1:3)     = [subset_one_minus_v(:,3)/subset_Varian_approx_ratio , subset_one_minus_v(:,3) , subset_one_minus_v(:,3)*subset_Varian_approx_ratio];
        
        if out_of_sample_required
            % update lower bound, approximate/exact and upper bound of 1-v, for out-of-sample vectors 
            subset_out_of_sample_one_minus_v_approx = nan(subset_obs_num, subset_obs_num-1, 9);
            subset_out_of_sample_one_minus_v_approx(:,:,1) = subset_out_of_sample_one_minus_v(:,:,1)/subset_Varian_approx_ratio;
            subset_out_of_sample_one_minus_v_approx(:,:,2) = subset_out_of_sample_one_minus_v(:,:,1);
            subset_out_of_sample_one_minus_v_approx(:,:,3) = subset_out_of_sample_one_minus_v(:,:,1)*subset_Varian_approx_ratio;
            subset_out_of_sample_one_minus_v_approx(:,:,4) = subset_out_of_sample_one_minus_v(:,:,2)/subset_Varian_approx_ratio;
            subset_out_of_sample_one_minus_v_approx(:,:,5) = subset_out_of_sample_one_minus_v(:,:,2);
            subset_out_of_sample_one_minus_v_approx(:,:,6) = subset_out_of_sample_one_minus_v(:,:,2)*subset_Varian_approx_ratio;
            subset_out_of_sample_one_minus_v_approx(:,:,7) = subset_out_of_sample_one_minus_v(:,:,3)/subset_Varian_approx_ratio;
            subset_out_of_sample_one_minus_v_approx(:,:,8) = subset_out_of_sample_one_minus_v(:,:,3);
            subset_out_of_sample_one_minus_v_approx(:,:,9) = subset_out_of_sample_one_minus_v(:,:,3)*subset_Varian_approx_ratio;
            out_of_sample_one_minus_v_per_subset{i} = subset_out_of_sample_one_minus_v_approx;
        end
    else
        
        % set lower bound, approximate/exact and upper bound of 1-v of this single observation to 0  
        one_minus_v_avg(current_subset, 1:3)     = [0 , 0 , 0];
        one_minus_v_meanssq(current_subset, 1:3) = [0 , 0 , 0];
        one_minus_v_max(current_subset, 1:3)     = [0 , 0 , 0];
        
        if out_of_sample_required
            % define the out-of-sample vector (an empty one though), just in case  
            out_of_sample_one_minus_v_per_subset{i} = nan(subset_obs_num, subset_obs_num-1, 9);
        end
    end
    
end



% the final Varian indices
average_var_lower = mean(one_minus_v_avg(:,1));
average_var_approx = mean(one_minus_v_avg(:,2));
average_var_upper = mean(one_minus_v_avg(:,3));
meanssq_var_lower = sqrt(meansqr(one_minus_v_meanssq(:,1)));
meanssq_var_approx = sqrt(meansqr(one_minus_v_meanssq(:,2)));
meanssq_var_upper = sqrt(meansqr(one_minus_v_meanssq(:,3)));
max_var_lower = max(one_minus_v_max(:,1));
max_var_approx = max(one_minus_v_max(:,2));
max_var_upper = max(one_minus_v_max(:,3));
Varian_Bounds = [average_var_lower , average_var_approx , average_var_upper , ...
                 meanssq_var_lower , meanssq_var_approx , meanssq_var_upper , ...
                 max_var_lower , max_var_approx , max_var_upper];
% we add "eps", since we assume that if we got to Varian calculation, then
% the subject does not satisfy GARP, and we want that in this case the
% Varian won't be 0, but will be a very very low positive number
Varian_Bounds = Varian_Bounds + eps;

% finalizing the in-sample residuals
average_var_in_sample_difference_approx = HPZ_Consistency_Indices_In_Sample_Difference_Residuals_Calc (one_minus_v_avg(:,2), @mean);
meanssq_var_in_sample_difference_approx = HPZ_Consistency_Indices_In_Sample_Difference_Residuals_Calc (one_minus_v_meanssq(:,2), @(x) sqrt(meansqr(x)));
max_var_in_sample_difference_approx     = HPZ_Consistency_Indices_In_Sample_Difference_Residuals_Calc (one_minus_v_max(:,2), @max);
Var_in_sample_component_residuals = [one_minus_v_avg(:,1) , one_minus_v_avg(:,2) , one_minus_v_avg(:,3) , ...
                                     one_minus_v_meanssq(:,1) , one_minus_v_meanssq(:,2) , one_minus_v_meanssq(:,3) , ...
                                     one_minus_v_max(:,1) , one_minus_v_max(:,2) , one_minus_v_max(:,3)];
Var_in_sample_difference_residuals = [average_var_in_sample_difference_approx , ...
                                      meanssq_var_in_sample_difference_approx , ...
                                      max_var_in_sample_difference_approx];

% the best 1-v for each of the aggregators (lower bound, approximate, upper bound)
one_minus_v = [one_minus_v_avg , one_minus_v_meanssq , one_minus_v_max];



% initializations
max_var_out_sample = nan(obs_num, 3);
average_var_out_sample = nan(obs_num, 3);
meanssq_var_out_sample = nan(obs_num, 3);
Varian_Mat = [];

% if out-of-sample residuals are required - combine them from the various subsets 
if out_of_sample_required
    
%     % define the waitbar
%     if (residuals_waitbar)
%         waitbar_name = char(strcat(HPZ_Constants.waitbar_name_estimation, {' '}, '(', HPZ_Constants.current_run_waitbar, {' '}, num2str(current_run), {' '}, HPZ_Constants.total_runs_waitbar, {' '}, num2str(total_runs),')'));
%         waitbar_msg = char(strcat(HPZ_Constants.waitbar_recovery, {' '}, num2str(subject_ID), {' '}, HPZ_Constants.waitbar_residuals_VARIAN));
%         new_bar_val = 0;
%         h_wb = wide_waitbar(new_bar_val, {waitbar_msg, ''}, waitbar_name, HPZ_Constants.waitbar_width_multiplier, [0,0.12]);
%     end
    
    
    for i=1:num_of_subsets
        
        current_subset = subsets_of_observations{i};
        subset_obs_num = length(current_subset);
        
        if subset_obs_num == 1
            
            % it does not involve in any relevant cycle -
            % it does not effect the varian index
            % we only need to multiply by n/(n-1) when it is an average
            max_var_out_sample(current_subset(1), 1:3) = [max_var_lower , max_var_approx , max_var_upper];
            average_var_out_sample(current_subset(1), 1:3) = [average_var_lower , average_var_approx , average_var_upper] * obs_num/(obs_num-1);
            meanssq_var_out_sample(current_subset(1), 1:3) = sqrt( ([meanssq_var_lower , meanssq_var_approx, meanssq_var_upper].^2) * obs_num/(obs_num-1) );
            
%             % updating the waitbar
%             if (residuals_waitbar)
%                 new_bar_val = new_bar_val + 1/obs_num;
%                 waitbar(new_bar_val, h_wb, {waitbar_msg , char(strcat({'Completed '}, num2str(new_bar_val*obs_num), {' observations out of '}, num2str(obs_num)))});
%             end
            
        else
            
            % we need to calculate the out-of-sample for each of these
            % observations *inside* the subset, then add it to what's
            % outside this subset
            % PART A - what's outside
            indices_outside_subset = true(1, obs_num);
            indices_outside_subset(current_subset) = false;
            indices_outside_subset = find(indices_outside_subset);
            % PART B - what's inside (+ combining inside & outside)
            subset_out_of_sample_one_minus_v = out_of_sample_one_minus_v_per_subset{i};
%             subset_expenditure = expenditure(current_subset , current_subset);
%             subset_identical_choice = identical_choice(current_subset , current_subset);
%             subset_SDRP = SDRP(current_subset , current_subset);
            %REDUNDANT: subset_transitive_SDRP = transitive_SDRP(current_subset , current_subset);
            for j=1:subset_obs_num
%                 % creating truncated matrices
%                 remaining_indices = [1:(j-1) , (j+1):length(current_subset)];
%                 truncated_expenditure = subset_expenditure(remaining_indices , remaining_indices);
%                 truncated_identical_choice = subset_identical_choice(remaining_indices , remaining_indices);
%                 truncated_SDRP = subset_SDRP(remaining_indices , remaining_indices);
%                 [~, truncated_transitive_SDRP, ~] = GARP_based_on_DRP_and_SDRP(truncated_SDRP, truncated_SDRP);
%                 %WRONG: truncated_transitive_SDRP = subset_transitive_SDRP(remaining_indices , remaining_indices);
%                 % we are only interested in SDRP relations that are involved in any cycles
%                 [~, truncated_SDRP] = find_relevant_relations(truncated_SDRP, truncated_SDRP);
%                 obs_involved_in_cycle = true(1, length(remaining_indices));
%                 for k=1:length(remaining_indices)
%                     if ~any(truncated_SDRP(k , :))
%                         % then this observation is not involved in any cycle in the truncated data, so we discard it 
%                         obs_involved_in_cycle(k) = false;
%                     end
%                 end
%                 % update remaining indices, after dropping obs that are not part of any cycle 
%                 remaining_indices = find(obs_involved_in_cycle);
%                 num_dropped = sum(~obs_involved_in_cycle);
%                 if isempty(remaining_indices)
%                     % then all out-of-sample values are 0
%                     out_sample_one_minus_v = zeros(num_dropped, 9);
%                     % WRONG CODE:
%                     %max_var_out_sample(current_subset(j),1) = 0;
%                     %max_var_out_sample(current_subset(j),2) = 0;
%                     %average_var_out_sample(current_subset(j),1) = 0;
%                     %average_var_out_sample(current_subset(j),2) = 0;
%                     %meanssq_var_out_sample(current_subset(j),1) = 0;
%                     %meanssq_var_out_sample(current_subset(j),2) = 0;
%                 else
%                     % re-creating truncated matrices
%                     truncated_expenditure = truncated_expenditure(remaining_indices , remaining_indices);
%                     truncated_identical_choice = truncated_identical_choice(remaining_indices , remaining_indices);
%                     truncated_transitive_SDRP = truncated_transitive_SDRP(remaining_indices , remaining_indices);
%                     % calculating varian for the truncated data
%                     [~, ~, ~, out_sample_one_minus_v, ~] = HPZ_Varian_Manager ([1,0,0,0,VARIAN_flags(5:7)], truncated_expenditure, truncated_identical_choice, index_threshold, truncated_transitive_SDRP, Varian_algorithm_settings, false, 0, 0, 0);
%                     % we need to add the dropped observations, with one-minus-v = 0 for each   
%                     out_sample_one_minus_v = [out_sample_one_minus_v ; zeros(num_dropped, 9)]; %#ok<AGROW>
%                     if (obs_num - 1) ~= size(out_sample_one_minus_v,1) + length(indices_outside_subset)
%                         error('There is a crucial bug in Varian out-of-sample residuals');
%                     end
%                 end
                out_sample_one_minus_v = subset_out_of_sample_one_minus_v(j,:,:);
                % now we assign the results
                average_var_out_sample(current_subset(j),1) = mean([out_sample_one_minus_v(:,1)' , one_minus_v_avg(indices_outside_subset,1)']);
                average_var_out_sample(current_subset(j),2) = mean([out_sample_one_minus_v(:,2)' , one_minus_v_avg(indices_outside_subset,2)']);
                average_var_out_sample(current_subset(j),3) = mean([out_sample_one_minus_v(:,3)' , one_minus_v_avg(indices_outside_subset,3)']);
                meanssq_var_out_sample(current_subset(j),1) = sqrt(mean([out_sample_one_minus_v(:,4)' , one_minus_v_meanssq(indices_outside_subset,1)'].^2));
                meanssq_var_out_sample(current_subset(j),2) = sqrt(mean([out_sample_one_minus_v(:,5)' , one_minus_v_meanssq(indices_outside_subset,2)'].^2));
                meanssq_var_out_sample(current_subset(j),3) = sqrt(mean([out_sample_one_minus_v(:,6)' , one_minus_v_meanssq(indices_outside_subset,3)'].^2));
                max_var_out_sample(current_subset(j),1) = max([out_sample_one_minus_v(:,7)' , one_minus_v_max(indices_outside_subset,1)']);
                max_var_out_sample(current_subset(j),2) = max([out_sample_one_minus_v(:,8)' , one_minus_v_max(indices_outside_subset,2)']);
                max_var_out_sample(current_subset(j),3) = max([out_sample_one_minus_v(:,9)' , one_minus_v_max(indices_outside_subset,3)']);
                %max_var_out_sample(current_subset(j)) = max(remaining_max_var, out_sample_Var(1));
                %average_var_out_sample(current_subset(j)) = ( out_sample_Var(2)*(length(current_subset) - 1) + remaining_average_var*(obs_num - length(current_subset)) ) / (obs_num - 1);
                %meanssq_var_out_sample(current_subset(j)) = sqrt( ( (out_sample_Var(3)^2)*(length(current_subset) - 1) + (remaining_meanssq_var^2)*(obs_num - length(current_subset)) ) / (obs_num - 1) );
            
%                 % updating the waitbar
%                 if (residuals_waitbar)
%                     new_bar_val = new_bar_val + 1/obs_num;
%                     waitbar(new_bar_val, h_wb, {waitbar_msg , char(strcat({'Completed '}, num2str(new_bar_val*obs_num), {' observations out of '}, num2str(obs_num)))});
%                 end
                
            end
            
        end
        
    end
    
%     % close the waitbar
%     if (residuals_waitbar)
%         close(h_wb);
%     end
    
end

% finalizing the out-of-sample residuals
Var_out_of_sample_residuals = [average_var_out_sample , meanssq_var_out_sample , max_var_out_sample];



% avg = one_minus_v_avg'
% meanssq = one_minus_v_meanssq'
% maximum = one_minus_v_max'





%% VARIAN residuals - assignment to matrix
if VARIAN_flags(2)

    % initialization of residuals matrix
    % for in-sample we need 18 columns: 
    %   component in-sample = (num aggregates) * [(component in-sample)*(exact , lower , approximate , upper) + (difference in-sample)*(exact , approximate)]    
    % for out-of-sample we need 36 columns:
    %   (num aggregates) * [(exact , lower , approximate, upper)*(full , out-of-sample) + (exact , approximate)*(residual , normalized-residual)]   
    num_of_columns = sum(VARIAN_flags(5:7))*6*VARIAN_flags(3) + sum(VARIAN_flags(5:7))*12*VARIAN_flags(4);
    Varian_Mat = nan(obs_num, num_of_columns);
    
    % initialization of column counter
    col_counter = 1;
    
    % in sample
    if VARIAN_flags(3)

        % we assign each observation its residuals (mean and meanssq) 
        
        if VARIAN_flags(5) % MEAN
            for i=1:obs_num
                if all(Var_in_sample_component_residuals(i,1) == Var_in_sample_component_residuals(i,3))
                    % it is exact
                    Varian_Mat(i, col_counter)     = Var_in_sample_component_residuals(i,2);
                    Varian_Mat(i, col_counter + 1) = NaN;
                    Varian_Mat(i, col_counter + 2) = NaN;
                    Varian_Mat(i, col_counter + 3) = NaN;
                    Varian_Mat(i, col_counter + 4) = Var_in_sample_difference_residuals(i,1);
                    Varian_Mat(i, col_counter + 5) = NaN;
                else
                    % it is not exact
                    Varian_Mat(i, col_counter)     = NaN;
                    Varian_Mat(i, col_counter + 1) = Var_in_sample_component_residuals(i,1);
                    Varian_Mat(i, col_counter + 2) = Var_in_sample_component_residuals(i,2);
                    Varian_Mat(i, col_counter + 3) = Var_in_sample_component_residuals(i,3);
                    Varian_Mat(i, col_counter + 4) = NaN;
                    Varian_Mat(i, col_counter + 5) = Var_in_sample_difference_residuals(i,1);
                end
            end
            % update the counter
            col_counter = col_counter + 6;
        end
        if VARIAN_flags(6) % AVGSSQ
            for i=1:obs_num
                if all(Var_in_sample_component_residuals(i,4) == Var_in_sample_component_residuals(i,6))
                    % it is exact
                    Varian_Mat(i, col_counter)     = Var_in_sample_component_residuals(i,5);
                    Varian_Mat(i, col_counter + 1) = NaN;
                    Varian_Mat(i, col_counter + 2) = NaN;
                    Varian_Mat(i, col_counter + 3) = NaN;
                    Varian_Mat(i, col_counter + 4) = Var_in_sample_difference_residuals(i,2);
                    Varian_Mat(i, col_counter + 5) = NaN;
                else
                    % it is not exact
                    Varian_Mat(i, col_counter)     = NaN;
                    Varian_Mat(i, col_counter + 1) = Var_in_sample_component_residuals(i,4);
                    Varian_Mat(i, col_counter + 2) = Var_in_sample_component_residuals(i,5);
                    Varian_Mat(i, col_counter + 3) = Var_in_sample_component_residuals(i,6);
                    Varian_Mat(i, col_counter + 4) = NaN;
                    Varian_Mat(i, col_counter + 5) = Var_in_sample_difference_residuals(i,2);
                end
            end
            % update the counter
            col_counter = col_counter + 6;
        end
        if VARIAN_flags(7) % MAX        
            for i=1:obs_num
                if all(Var_in_sample_component_residuals(i,7) == Var_in_sample_component_residuals(i,9))
                    % it is exact
                    Varian_Mat(i, col_counter)     = Var_in_sample_component_residuals(i,8);
                    Varian_Mat(i, col_counter + 1) = NaN;
                    Varian_Mat(i, col_counter + 2) = NaN;
                    Varian_Mat(i, col_counter + 3) = NaN;
                    Varian_Mat(i, col_counter + 4) = Var_in_sample_difference_residuals(i,3);
                    Varian_Mat(i, col_counter + 5) = NaN;
                else
                    % it is not exact
                    Varian_Mat(i, col_counter)     = NaN;
                    Varian_Mat(i, col_counter + 1) = Var_in_sample_component_residuals(i,7);
                    Varian_Mat(i, col_counter + 2) = Var_in_sample_component_residuals(i,8);
                    Varian_Mat(i, col_counter + 3) = Var_in_sample_component_residuals(i,9);
                    Varian_Mat(i, col_counter + 4) = NaN;
                    Varian_Mat(i, col_counter + 5) = Var_in_sample_difference_residuals(i,3);
                end
            end
            % update the counter
            col_counter = col_counter + 6;
        end
            
        
    end
    
    % out of sample
    if VARIAN_flags(4)

        % we assign each observation its residuals (mean and meanssq) 
        
        % assigning to the matrix

        if VARIAN_flags(5) % MEAN 

            % normalized indices (relevant for MEAN and AVGSSQ):
            VARIAN_Mean_norm = Var_out_of_sample_residuals(:,2) * (obs_num-1)/obs_num;

            % full index
            if Varian_Bounds(1) == Varian_Bounds(3)
                % it is exact
                Varian_Mat(:, col_counter) = Varian_Bounds(2);
                Varian_Mat(:, col_counter + 1) = NaN;
                Varian_Mat(:, col_counter + 2) = NaN;
                Varian_Mat(:, col_counter + 3) = NaN;
            else
                % it is not exact
                Varian_Mat(:, col_counter) = NaN;
                Varian_Mat(:, col_counter + 1) = Varian_Bounds(1);
                Varian_Mat(:, col_counter + 2) = Varian_Bounds(2);
                Varian_Mat(:, col_counter + 3) = Varian_Bounds(3);
            end
            % update the counter
            col_counter = col_counter + 4;

            % partial index
            for i=1:obs_num
                if Var_out_of_sample_residuals(i,1) == Var_out_of_sample_residuals(i,3)
                    % it is exact
                    Varian_Mat(i, col_counter)     = Var_out_of_sample_residuals(i,2);
                    Varian_Mat(i, col_counter + 1) = NaN;
                    Varian_Mat(i, col_counter + 2) = NaN;
                    Varian_Mat(i, col_counter + 3) = NaN;
                else
                    % it is not exact
                    Varian_Mat(i, col_counter)     = NaN;
                    Varian_Mat(i, col_counter + 1) = Var_out_of_sample_residuals(i,1);
                    Varian_Mat(i, col_counter + 2) = Var_out_of_sample_residuals(i,2);
                    Varian_Mat(i, col_counter + 3) = Var_out_of_sample_residuals(i,3);
                end
            end
            % update the counter
            col_counter = col_counter + 4;

            % difference & normalized difference
            for i=1:obs_num
                % (we added this check, because the whole point of the normalized difference is that it is either positive or zero, but due to calculation issues it sometimes resulted stuff like "2.77555756156289E-16")   
                normalized_difference = Varian_Bounds(2) - VARIAN_Mean_norm(i);
                if abs(normalized_difference) < 10^(-15)
                    normalized_difference = 0;
                end
                if Varian_Bounds(1) == Varian_Bounds(3) && Var_out_of_sample_residuals(i,1) == Var_out_of_sample_residuals(i,3)
                    % it is exact
                    Varian_Mat(i, col_counter)     = Varian_Bounds(2) - Var_out_of_sample_residuals(i,2);
                    Varian_Mat(i, col_counter + 1) = normalized_difference;
                    Varian_Mat(i, col_counter + 2) = NaN;
                    Varian_Mat(i, col_counter + 3) = NaN;
                else
                    % it is not exact
                    Varian_Mat(i, col_counter)     = NaN;
                    Varian_Mat(i, col_counter + 1) = NaN;
                    Varian_Mat(i, col_counter + 2) = Varian_Bounds(2) - Var_out_of_sample_residuals(i,2);
                    Varian_Mat(i, col_counter + 3) = normalized_difference;
                end
            end
            % update the counter
            col_counter = col_counter + 4;
        end

        if VARIAN_flags(6) % AVGSSQ

            % normalized indices (relevant for MEAN and AVGSSQ):
            VARIAN_AVGSSQ_norm = sqrt(Var_out_of_sample_residuals(:,5).^2 * (obs_num-1)/obs_num);

            % full index
            if Varian_Bounds(4) == Varian_Bounds(6)
                % it is exact
                Varian_Mat(:, col_counter)     = Varian_Bounds(5);
                Varian_Mat(:, col_counter + 1) = NaN;
                Varian_Mat(:, col_counter + 2) = NaN;
                Varian_Mat(:, col_counter + 3) = NaN;
            else
                % it is not exact
                Varian_Mat(:, col_counter)     = NaN;
                Varian_Mat(:, col_counter + 1) = Varian_Bounds(4);
                Varian_Mat(:, col_counter + 2) = Varian_Bounds(5);
                Varian_Mat(:, col_counter + 3) = Varian_Bounds(6);
            end
            % update the counter
            col_counter = col_counter + 4;

            % partial index
            for i=1:obs_num
                if Var_out_of_sample_residuals(i,4) == Var_out_of_sample_residuals(i,6)
                    % it is exact
                    Varian_Mat(i, col_counter)     = Var_out_of_sample_residuals(i,5);
                    Varian_Mat(i, col_counter + 1) = NaN;
                    Varian_Mat(i, col_counter + 2) = NaN;
                    Varian_Mat(i, col_counter + 3) = NaN;
                else
                    % it is not exact
                    Varian_Mat(i, col_counter)     = NaN;
                    Varian_Mat(i, col_counter + 1) = Var_out_of_sample_residuals(i,4);
                    Varian_Mat(i, col_counter + 2) = Var_out_of_sample_residuals(i,5);
                    Varian_Mat(i, col_counter + 3) = Var_out_of_sample_residuals(i,6);
                end
            end
            % update the counter
            col_counter = col_counter + 4;

            % difference & normalized difference
            for i=1:obs_num
                % (we added this check, because the whole point of the normalized difference is that it is either positive or zero, but due to calculation issues it sometimes resulted stuff like "2.77555756156289E-16")   
                normalized_difference = Varian_Bounds(5) - VARIAN_AVGSSQ_norm(i);
                if abs(normalized_difference) < 10^(-15)
                    normalized_difference = 0;
                end
                if Varian_Bounds(4) == Varian_Bounds(6) && Var_out_of_sample_residuals(i,4) == Var_out_of_sample_residuals(i,6)
                    % it is exact
                    Varian_Mat(i, col_counter)     = Varian_Bounds(5) - Var_out_of_sample_residuals(i,5);
                    Varian_Mat(i, col_counter + 1) = normalized_difference;
                    Varian_Mat(i, col_counter + 2) = NaN;
                    Varian_Mat(i, col_counter + 3) = NaN;
                else
                    % it is not exact
                    Varian_Mat(i, col_counter)     = NaN;
                    Varian_Mat(i, col_counter + 1) = NaN;
                    Varian_Mat(i, col_counter + 2) = Varian_Bounds(5) - Var_out_of_sample_residuals(i,5);
                    Varian_Mat(i, col_counter + 3) = normalized_difference;
                end
            end
            % update the counter
            col_counter = col_counter + 4;
        end

        if VARIAN_flags(7) % MAX

            % normalized indices (relevant for MEAN and AVGSSQ):
            VARIAN_Max_norm = Var_out_of_sample_residuals(:,8);

            % full index
            if Varian_Bounds(7) == Varian_Bounds(9)
                % it is exact
                Varian_Mat(:, col_counter)     = Varian_Bounds(8);
                Varian_Mat(:, col_counter + 1) = NaN;
                Varian_Mat(:, col_counter + 2) = NaN;
                Varian_Mat(:, col_counter + 3) = NaN;
            else
                % it is not exact
                Varian_Mat(:, col_counter)     = NaN;
                Varian_Mat(:, col_counter + 1) = Varian_Bounds(7);
                Varian_Mat(:, col_counter + 2) = Varian_Bounds(8);
                Varian_Mat(:, col_counter + 3) = Varian_Bounds(9);
            end
            % update the counter
            col_counter = col_counter + 4;

            % partial index
            for i=1:obs_num
                if Var_out_of_sample_residuals(i,7) == Var_out_of_sample_residuals(i,9)
                    % it is exact
                    Varian_Mat(i, col_counter)     = Var_out_of_sample_residuals(i,8);
                    Varian_Mat(i, col_counter + 1) = NaN;
                    Varian_Mat(i, col_counter + 2) = NaN;
                    Varian_Mat(i, col_counter + 3) = NaN;
                else
                    % it is not exact
                    Varian_Mat(i, col_counter)     = NaN;
                    Varian_Mat(i, col_counter + 1) = Var_out_of_sample_residuals(i,7);
                    Varian_Mat(i, col_counter + 2) = Var_out_of_sample_residuals(i,8);
                    Varian_Mat(i, col_counter + 3) = Var_out_of_sample_residuals(i,9);
                end
            end
            % update the counter
            col_counter = col_counter + 4;

            % difference & normalized difference
            for i=1:obs_num
                % (we added this check, because the whole point of the normalized difference is that it is either positive or zero, but due to calculation issues it sometimes resulted stuff like "2.77555756156289E-16")   
                normalized_difference = Varian_Bounds(8) - VARIAN_Max_norm(i);
                if abs(normalized_difference) < 10^(-15)
                    normalized_difference = 0;
                end
                if Varian_Bounds(7) == Varian_Bounds(9) && Var_out_of_sample_residuals(i,7) == Var_out_of_sample_residuals(i,9)
                    % it is exact
                    Varian_Mat(i, col_counter)     = Varian_Bounds(8) - Var_out_of_sample_residuals(i,8);
                    Varian_Mat(i, col_counter + 1) = normalized_difference;
                    Varian_Mat(i, col_counter + 2) = NaN;
                    Varian_Mat(i, col_counter + 3) = NaN;
                else
                    % it is not exact
                    Varian_Mat(i, col_counter)     = NaN;
                    Varian_Mat(i, col_counter + 1) = NaN;
                    Varian_Mat(i, col_counter + 2) = Varian_Bounds(8) - Var_out_of_sample_residuals(i,8);
                    Varian_Mat(i, col_counter + 3) = normalized_difference;
                end
            end
            % update the counter
            col_counter = col_counter + 4; %#ok<NASGU>
        end   
        
    end
    
end



end