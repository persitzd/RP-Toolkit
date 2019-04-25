function [Varian_Bounds, Var_in_sample_residuals, Var_out_sample_residuals, one_minus_v, Varian_Mat] = HPZ_Varian_Manager (VARIAN_flags, expenditure, identical_choice, index_threshold, SDRP, Varian_algorithm_settings, residuals_waitbar, current_run, total_runs, subject_ID)

% expenditure1 = expenditure;
% % for i=1:4
% %     expenditure1(i,i) = expenditure1(i,i) * (1-10^(-13));
% % end
% expenditure1(2,2) = expenditure1(2,2) * (1 - 0.03517); % 0.03517
% expenditure1(3,3) = expenditure1(3,3) * (1 - 0.11703); %0.15958 0.11703
% [GARP, ~, ~, ~, ~] = GARP_based_on_expenditures(expenditure1, identical_choice, index_threshold);
% garp = sum(sum(GARP))
% [frow,fcol] = find(GARP)

% number of observations
[obs_num , ~] = size(expenditure);



% we calculate SDRP with transitivity, instead of DRP with transitivity,
% since in the case of Varian weak relations are of no interest
[~, transitive_SDRP, ~] = GARP_based_on_DRP_and_SDRP(SDRP, SDRP);

% divide the observations to smallest possible subsets such that
% observations from different subsets never belong to one cycle
% (a cycle of strictly revealed relations in this case)
[subsets_of_observations, num_of_subsets] = HPZ_Indices_Problem_Distribute(transitive_SDRP, transitive_SDRP);

% initialization (it must be initialized to 0)
one_minus_v_max = zeros(obs_num, 2);
one_minus_v_avg = zeros(obs_num, 2);
one_minus_v_meanssq = zeros(obs_num, 2);

% now we calculate the Varian separately for each subset
for i=1:num_of_subsets
    
    current_subset = subsets_of_observations{i};
    
    if length(current_subset) > 1
        
        % extract the relevant expenditures and identical choices
        subset_expenditure = expenditure(current_subset , current_subset);
        subset_identical_choice = identical_choice(current_subset , current_subset);
        
        % calculate the Varian index for this subset of choices
        [~, subset_Varian_approx_ratio, ~, subset_one_minus_v] = HPZ_Varian_index_based_on_Houtman_Maks (subset_expenditure, subset_identical_choice, index_threshold, Varian_algorithm_settings);
        % Oriel Notes: for debugging:
        %[~, subset_Varian_exact, ~, subset_one_minus_v] = HPZ_Varian_efficiency_index (subset_expenditure, subset_identical_choice, index_threshold);
        %subset_Varian_approx_ratio = subset_Varian_exact;
        
        % update lower and upper bounds of 1-v
        one_minus_v_max(current_subset, 1:2)     = [subset_one_minus_v(:,1)/subset_Varian_approx_ratio , subset_one_minus_v(:,1)*subset_Varian_approx_ratio];
        one_minus_v_avg(current_subset, 1:2)     = [subset_one_minus_v(:,2)/subset_Varian_approx_ratio , subset_one_minus_v(:,2)*subset_Varian_approx_ratio];
        one_minus_v_meanssq(current_subset, 1:2) = [subset_one_minus_v(:,3)/subset_Varian_approx_ratio , subset_one_minus_v(:,3)*subset_Varian_approx_ratio];
        
    end
    
end



% the final Varian indices
max_var_lower = max(one_minus_v_max(:,1));
max_var_upper = max(one_minus_v_max(:,2));
average_var_lower = mean(one_minus_v_avg(:,1));
average_var_upper = mean(one_minus_v_avg(:,2));
meanssq_var_lower = sqrt(meansqr(one_minus_v_meanssq(:,1)));
meanssq_var_upper = sqrt(meansqr(one_minus_v_meanssq(:,2)));
Varian_Bounds = [max_var_lower, max_var_upper , ...
                 average_var_lower , average_var_upper , ...
                 meanssq_var_lower , meanssq_var_upper];
% we add "eps", since we assume that if we got to Varian calculation, then
% the subject does not satisfy GARP, and we want that in this case the
% Varian won't be 0, but will be a very very low positive number
Varian_Bounds = Varian_Bounds + eps;

% finalizing the in-sample residuals
max_var_in_sample_difference_lower = HPZ_Consistency_Indices_In_Sample_Difference_Residuals_Calc (one_minus_v_max(:,1), @max);
max_var_in_sample_difference_upper = HPZ_Consistency_Indices_In_Sample_Difference_Residuals_Calc (one_minus_v_max(:,2), @max);
average_var_in_sample_difference_lower = HPZ_Consistency_Indices_In_Sample_Difference_Residuals_Calc (one_minus_v_avg(:,1), @mean);
average_var_in_sample_difference_upper = HPZ_Consistency_Indices_In_Sample_Difference_Residuals_Calc (one_minus_v_avg(:,2), @mean);
meanssq_var_in_sample_difference_lower = HPZ_Consistency_Indices_In_Sample_Difference_Residuals_Calc (one_minus_v_meanssq(:,1), @(x) sqrt(meansqr(x)));
meanssq_var_in_sample_difference_upper = HPZ_Consistency_Indices_In_Sample_Difference_Residuals_Calc (one_minus_v_meanssq(:,2), @(x) sqrt(meansqr(x)));
Var_in_sample_residuals = [one_minus_v_max(:,1) , one_minus_v_max(:,2) , ...
                           one_minus_v_avg(:,1) , one_minus_v_avg(:,2) , ...
                           one_minus_v_meanssq(:,1) , one_minus_v_meanssq(:,2) , ...
                           max_var_in_sample_difference_lower , max_var_in_sample_difference_upper , ...
                           average_var_in_sample_difference_lower , average_var_in_sample_difference_upper , ...
                           meanssq_var_in_sample_difference_lower , meanssq_var_in_sample_difference_upper];

% the best 1-v for each of the aggregators (lower and upper bounds)
one_minus_v = [one_minus_v_max, one_minus_v_avg, one_minus_v_meanssq];



% initializations
max_var_out_sample = zeros(obs_num, 2);
average_var_out_sample = zeros(obs_num, 2);
meanssq_var_out_sample = zeros(obs_num, 2);
Varian_Mat = [];

% if out-of-sample residuals are required - calculate them
if VARIAN_flags(2) && VARIAN_flags(4)
    
    % define the waitbar
    if (residuals_waitbar)
        waitbar_name = char(strcat(HPZ_Constants.waitbar_name_estimation, {' '}, '(', HPZ_Constants.current_run_waitbar, {' '}, num2str(current_run), {' '}, HPZ_Constants.total_runs_waitbar, {' '}, num2str(total_runs),')'));
        waitbar_msg = char(strcat(HPZ_Constants.waitbar_recovery, {' '}, num2str(subject_ID), {' '}, HPZ_Constants.waitbar_residuals_VARIAN));
        new_bar_val = 0;
        h_wb = wide_waitbar(new_bar_val, {waitbar_msg, ''}, waitbar_name, HPZ_Constants.waitbar_width_multiplier, [0,0.12]);
    end
    
    
    for i=1:num_of_subsets
        
        current_subset = subsets_of_observations{i};
        
        if length(current_subset) == 1
            
            % it does not involve in any relevant cycle -
            % it does not effect the varian index
            % we only need to multiply by n/(n-1) when it is an average
            max_var_out_sample(current_subset(1), 1:2) = [max_var_lower , max_var_upper];
            average_var_out_sample(current_subset(1), 1:2) = [average_var_lower , average_var_upper] * obs_num/(obs_num-1);
            meanssq_var_out_sample(current_subset(1), 1:2) = sqrt( ([meanssq_var_lower , meanssq_var_upper].^2) * obs_num/(obs_num-1) );
            
            % updating the waitbar
            if (residuals_waitbar)
                new_bar_val = new_bar_val + 1/obs_num;
                waitbar(new_bar_val, h_wb, {waitbar_msg , char(strcat({'Completed '}, num2str(new_bar_val*obs_num), {' observations out of '}, num2str(obs_num)))});
            end
            
        else
            
            % we need to calculate the out-of-sample for each of these
            % observations *inside* the subset, then add it to what's
            % outside this subset
            % PART A - what's outside
            indices_outside_subset = true(1, obs_num);
            indices_outside_subset(current_subset) = false;
            indices_outside_subset = find(indices_outside_subset);
            % PART B - what's inside (+ combining inside & outside)
            subset_expenditure = expenditure(current_subset , current_subset);
            subset_identical_choice = identical_choice(current_subset , current_subset);
            subset_SDRP = SDRP(current_subset , current_subset);
            %REDUNDANT: subset_transitive_SDRP = transitive_SDRP(current_subset , current_subset);
            for j=1:length(current_subset)
                % creating truncated matrices
                remaining_indices = [1:(j-1) , (j+1):length(current_subset)];
                truncated_expenditure = subset_expenditure(remaining_indices , remaining_indices);
                truncated_identical_choice = subset_identical_choice(remaining_indices , remaining_indices);
                truncated_SDRP = subset_SDRP(remaining_indices , remaining_indices);
                [~, truncated_transitive_SDRP, ~] = GARP_based_on_DRP_and_SDRP(truncated_SDRP, truncated_SDRP);
                %WRONG: truncated_transitive_SDRP = subset_transitive_SDRP(remaining_indices , remaining_indices);
                % we are only interested in SDRP relations that are involved in any cycles
                [~, truncated_SDRP] = find_relevant_relations(truncated_SDRP, truncated_SDRP);
                obs_involved_in_cycle = true(1, length(remaining_indices));
                for k=1:length(remaining_indices)
                    if ~any(truncated_SDRP(k , :))
                        % then this observation is not involved in any cycle in the truncated data, so we discard it 
                        obs_involved_in_cycle(k) = false;
                    end
                end
                % update remaining indices, after dropping obs that are not part of any cycle 
                remaining_indices = find(obs_involved_in_cycle);
                num_dropped = sum(~obs_involved_in_cycle);
                if isempty(remaining_indices)
                    % then all out-of-sample values are 0
                    out_sample_one_minus_v = zeros(num_dropped, 6);
                    % WRONG CODE:
                    %max_var_out_sample(current_subset(j),1) = 0;
                    %max_var_out_sample(current_subset(j),2) = 0;
                    %average_var_out_sample(current_subset(j),1) = 0;
                    %average_var_out_sample(current_subset(j),2) = 0;
                    %meanssq_var_out_sample(current_subset(j),1) = 0;
                    %meanssq_var_out_sample(current_subset(j),2) = 0;
                else
                    % re-creating truncated matrices
                    truncated_expenditure = truncated_expenditure(remaining_indices , remaining_indices);
                    truncated_identical_choice = truncated_identical_choice(remaining_indices , remaining_indices);
                    truncated_transitive_SDRP = truncated_transitive_SDRP(remaining_indices , remaining_indices);
                    % calculating varian for the truncated data
                    [~, ~, ~, out_sample_one_minus_v, ~] = HPZ_Varian_Manager ([1,0,0,0], truncated_expenditure, truncated_identical_choice, index_threshold, truncated_transitive_SDRP, Varian_algorithm_settings, false, 0, 0, 0);
                    % we need to add the dropped observations, with one-minus-v = 0 for each   
                    out_sample_one_minus_v = [out_sample_one_minus_v ; zeros(num_dropped, 6)]; %#ok<AGROW>
                    if (obs_num - 1) ~= size(out_sample_one_minus_v,1) + length(indices_outside_subset)
                        error('There is a crucial bug in Varian out-of-sample residuals');
                    end
                end
                % now we assign the results
                max_var_out_sample(current_subset(j),1) = max([out_sample_one_minus_v(:,1)' , one_minus_v_max(indices_outside_subset,1)']);
                max_var_out_sample(current_subset(j),2) = max([out_sample_one_minus_v(:,2)' , one_minus_v_max(indices_outside_subset,2)']);
                average_var_out_sample(current_subset(j),1) = mean([out_sample_one_minus_v(:,3)' , one_minus_v_avg(indices_outside_subset,1)']);
                average_var_out_sample(current_subset(j),2) = mean([out_sample_one_minus_v(:,4)' , one_minus_v_avg(indices_outside_subset,2)']);
                meanssq_var_out_sample(current_subset(j),1) = sqrt(meansqr([out_sample_one_minus_v(:,5)' , one_minus_v_meanssq(indices_outside_subset,1)']));
                meanssq_var_out_sample(current_subset(j),2) = sqrt(meansqr([out_sample_one_minus_v(:,6)' , one_minus_v_meanssq(indices_outside_subset,2)']));
                %max_var_out_sample(current_subset(j)) = max(remaining_max_var, out_sample_Var(1));
                %average_var_out_sample(current_subset(j)) = ( out_sample_Var(2)*(length(current_subset) - 1) + remaining_average_var*(obs_num - length(current_subset)) ) / (obs_num - 1);
                %meanssq_var_out_sample(current_subset(j)) = sqrt( ( (out_sample_Var(3)^2)*(length(current_subset) - 1) + (remaining_meanssq_var^2)*(obs_num - length(current_subset)) ) / (obs_num - 1) );
            
                % updating the waitbar
                if (residuals_waitbar)
                    new_bar_val = new_bar_val + 1/obs_num;
                    waitbar(new_bar_val, h_wb, {waitbar_msg , char(strcat({'Completed '}, num2str(new_bar_val*obs_num), {' observations out of '}, num2str(obs_num)))});
                end
                
            end
            
        end
        
    end
    
    % close the waitbar
    if (residuals_waitbar)
        close(h_wb);
    end
    
end

% finalizing the out-of-sample residuals
Var_out_sample_residuals = [max_var_out_sample , average_var_out_sample , meanssq_var_out_sample];



% avg = one_minus_v_avg'
% meanssq = one_minus_v_meanssq'
% maximum = one_minus_v_max'





%% VARIAN residuals - assignment to matrix
if VARIAN_flags(2)

    % initialization of residuals matrix
    % for in-sample we need 18 columns: 
    %   2 types of in-sample * 3 aggregates * (exact , lower , upper)   
    % for out-of-sample we need 36 columns:
    %   3 aggregates * (exact , lower , upper) * (full , out-of-sample , residual , normalized-residual)   
    num_of_columns = 18*VARIAN_flags(3) + 36*VARIAN_flags(4);
    Varian_Mat = zeros(obs_num, num_of_columns);

    % initialization of column counter
    col_counter = 1;
    
    % in sample
    if VARIAN_flags(3)

        % we assign each observation its residuals (mean and meanssq) 
        for i=1:obs_num
            % MAX
            if all(Var_in_sample_residuals(i,1) == Var_in_sample_residuals(i,2))
                % it is exact
                Varian_Mat(i, col_counter) = Var_in_sample_residuals(i,1);
                Varian_Mat(i, col_counter + 1) = NaN;
                Varian_Mat(i, col_counter + 2) = NaN;
                Varian_Mat(i, col_counter + 3) = Var_in_sample_residuals(i,6+1);
                Varian_Mat(i, col_counter + 4) = NaN;
                Varian_Mat(i, col_counter + 5) = NaN;
            else
                % it is not exact
                Varian_Mat(i, col_counter) = NaN;
                Varian_Mat(i, col_counter + 1) = Var_in_sample_residuals(i,1);
                Varian_Mat(i, col_counter + 2) = Var_in_sample_residuals(i,2);
                Varian_Mat(i, col_counter + 3) = NaN;
                Varian_Mat(i, col_counter + 4) = Var_in_sample_residuals(i,6+1);
                Varian_Mat(i, col_counter + 5) = Var_in_sample_residuals(i,6+2);
            end
            % MEAN 
            if all(Var_in_sample_residuals(i,3) == Var_in_sample_residuals(i,4))
                % it is exact
                Varian_Mat(i, col_counter + 6) = Var_in_sample_residuals(i,3);
                Varian_Mat(i, col_counter + 7) = NaN;
                Varian_Mat(i, col_counter + 8) = NaN;
                Varian_Mat(i, col_counter + 9) = Var_in_sample_residuals(i,6+3);
                Varian_Mat(i, col_counter + 10) = NaN;
                Varian_Mat(i, col_counter + 11) = NaN;
            else
                % it is not exact
                Varian_Mat(i, col_counter + 6) = NaN;
                Varian_Mat(i, col_counter + 7) = Var_in_sample_residuals(i,3);
                Varian_Mat(i, col_counter + 8) = Var_in_sample_residuals(i,4);
                Varian_Mat(i, col_counter + 9) = NaN;
                Varian_Mat(i, col_counter + 10) = Var_in_sample_residuals(i,6+3);
                Varian_Mat(i, col_counter + 11) = Var_in_sample_residuals(i,6+4);
            end
            % AVGSSQ
            if all(Var_in_sample_residuals(i,5) == Var_in_sample_residuals(i,6))
                % it is exact
                Varian_Mat(i, col_counter + 12) = Var_in_sample_residuals(i,5);
                Varian_Mat(i, col_counter + 13) = NaN;
                Varian_Mat(i, col_counter + 14) = NaN;
                Varian_Mat(i, col_counter + 15) = Var_in_sample_residuals(i,6+5);
                Varian_Mat(i, col_counter + 16) = NaN;
                Varian_Mat(i, col_counter + 17) = NaN;
            else
                % it is not exact
                Varian_Mat(i, col_counter + 12) = NaN;
                Varian_Mat(i, col_counter + 13) = Var_in_sample_residuals(i,5);
                Varian_Mat(i, col_counter + 14) = Var_in_sample_residuals(i,6);
                Varian_Mat(i, col_counter + 15) = NaN;
                Varian_Mat(i, col_counter + 16) = Var_in_sample_residuals(i,6+5);
                Varian_Mat(i, col_counter + 17) = Var_in_sample_residuals(i,6+6);
            end     
        end 
        % update the counter
        col_counter = col_counter + 18;
        
    end

    % out of sample
    if VARIAN_flags(4)

        % we assign each observation its residuals (mean and meanssq) 
        for i=1:obs_num
            % normalized indices (relevant for MEAN and AVGSSQ):
            VARIAN_Max_Lower_norm = Var_out_sample_residuals(i,1);
            VARIAN_Max_Upper_norm = Var_out_sample_residuals(i,2);
            VARIAN_Mean_Lower_norm = Var_out_sample_residuals(i,3) * (obs_num-1)/obs_num;
            VARIAN_Mean_Upper_norm = Var_out_sample_residuals(i,4) * (obs_num-1)/obs_num;
            VARIAN_AVGSSQ_Lower_norm = sqrt(Var_out_sample_residuals(i,5)^2 * (obs_num-1)/obs_num);
            VARIAN_AVGSSQ_Upper_norm = sqrt(Var_out_sample_residuals(i,6)^2 * (obs_num-1)/obs_num);

            % assigning to the matrix

            % MAX full index
            if Varian_Bounds(1) == Varian_Bounds(2)
                % it is exact
                Varian_Mat(i, col_counter) = Varian_Bounds(1);
                Varian_Mat(i, col_counter + 1) = NaN;
                Varian_Mat(i, col_counter + 2) = NaN;
            else
                % it is not exact
                Varian_Mat(i, col_counter) = NaN;
                Varian_Mat(i, col_counter + 1) = Varian_Bounds(1);
                Varian_Mat(i, col_counter + 2) = Varian_Bounds(2);
            end
            % MEAN full index
            if Varian_Bounds(3) == Varian_Bounds(4)
                % it is exact
                Varian_Mat(i, col_counter + 3) = Varian_Bounds(3);
                Varian_Mat(i, col_counter + 4) = NaN;
                Varian_Mat(i, col_counter + 5) = NaN;
            else
                % it is not exact
                Varian_Mat(i, col_counter + 3) = NaN;
                Varian_Mat(i, col_counter + 4) = Varian_Bounds(3);
                Varian_Mat(i, col_counter + 5) = Varian_Bounds(4);
            end
            % AVGSSQ full index
            if Varian_Bounds(5) == Varian_Bounds(6)
                % it is exact
                Varian_Mat(i, col_counter + 6) = Varian_Bounds(5);
                Varian_Mat(i, col_counter + 7) = NaN;
                Varian_Mat(i, col_counter + 8) = NaN;
            else
                % it is not exact
                Varian_Mat(i, col_counter + 6) = NaN;
                Varian_Mat(i, col_counter + 7) = Varian_Bounds(5);
                Varian_Mat(i, col_counter + 8) = Varian_Bounds(6);
            end


            % MAX partial index
            if Var_out_sample_residuals(i,1) == Var_out_sample_residuals(i,2)
                % it is exact
                Varian_Mat(i, col_counter + 9) = Var_out_sample_residuals(i,1);
                Varian_Mat(i, col_counter + 10) = NaN;
                Varian_Mat(i, col_counter + 11) = NaN;
            else
                % it is not exact
                Varian_Mat(i, col_counter + 9) = NaN;
                Varian_Mat(i, col_counter + 10) = Var_out_sample_residuals(i,1);
                Varian_Mat(i, col_counter + 11) = Var_out_sample_residuals(i,2);
            end
            % MEAN partial index
            if Var_out_sample_residuals(i,3) == Var_out_sample_residuals(i,4)
                % it is exact
                Varian_Mat(i, col_counter + 12) = Var_out_sample_residuals(i,3);
                Varian_Mat(i, col_counter + 13) = NaN;
                Varian_Mat(i, col_counter + 14) = NaN;
            else
                % it is not exact
                Varian_Mat(i, col_counter + 12) = NaN;
                Varian_Mat(i, col_counter + 13) = Var_out_sample_residuals(i,3);
                Varian_Mat(i, col_counter + 14) = Var_out_sample_residuals(i,4);
            end
            % AVGSSQ partial index
            if Var_out_sample_residuals(i,5) == Var_out_sample_residuals(i,6)
                % it is exact
                Varian_Mat(i, col_counter + 15) = Var_out_sample_residuals(i,5);
                Varian_Mat(i, col_counter + 16) = NaN;
                Varian_Mat(i, col_counter + 17) = NaN;
            else
                % it is not exact
                Varian_Mat(i, col_counter + 15) = NaN;
                Varian_Mat(i, col_counter + 16) = Var_out_sample_residuals(i,5);
                Varian_Mat(i, col_counter + 17) = Var_out_sample_residuals(i,6);
            end


            % MAX difference & MAX normalized difference
            % (we added this check, because the whole point of the normalized difference is that it is either positive or zero, but due to calculation issues it sometimes resulted stuff like "2.77555756156289E-16")   
            normalized_difference_lower = Varian_Bounds(1) - VARIAN_Max_Lower_norm;
            normalized_difference_upper = Varian_Bounds(2) - VARIAN_Max_Upper_norm;
            if abs(normalized_difference_lower) < 10^(-15)
                normalized_difference_lower = 0;
            end
            if abs(normalized_difference_upper) < 10^(-15)
                normalized_difference_upper = 0;
            end
            if Varian_Bounds(1) == Varian_Bounds(2) && Var_out_sample_residuals(i,1) == Var_out_sample_residuals(i,2)
                % it is exact
                Varian_Mat(i, col_counter + 18) = Varian_Bounds(1) - Var_out_sample_residuals(i,1);
                Varian_Mat(i, col_counter + 19) = normalized_difference_lower;
                Varian_Mat(i, col_counter + 20) = NaN;
                Varian_Mat(i, col_counter + 21) = NaN;
                Varian_Mat(i, col_counter + 22) = NaN;
                Varian_Mat(i, col_counter + 23) = NaN;
            else
                % it is not exact
                Varian_Mat(i, col_counter + 18) = NaN;
                Varian_Mat(i, col_counter + 19) = NaN;
                Varian_Mat(i, col_counter + 20) = Varian_Bounds(1) - Var_out_sample_residuals(i,1);
                Varian_Mat(i, col_counter + 21) = Varian_Bounds(2) - Var_out_sample_residuals(i,2);
                Varian_Mat(i, col_counter + 22) = normalized_difference_lower;
                Varian_Mat(i, col_counter + 23) = normalized_difference_upper;
            end

            % MEAN difference & MEAN normalized difference
            % (we added this check, because the whole point of the normalized difference is that it is either positive or zero, but due to calculation issues it sometimes resulted stuff like "2.77555756156289E-16")   
            normalized_difference_lower = Varian_Bounds(3) - VARIAN_Mean_Lower_norm;
            normalized_difference_upper = Varian_Bounds(4) - VARIAN_Mean_Upper_norm;
            if abs(normalized_difference_lower) < 10^(-15)
                normalized_difference_lower = 0;
            end
            if abs(normalized_difference_upper) < 10^(-15)
                normalized_difference_upper = 0;
            end
            if Varian_Bounds(3) == Varian_Bounds(4) && Var_out_sample_residuals(i,3) == Var_out_sample_residuals(i,4)
                % it is exact
                Varian_Mat(i, col_counter + 24) = Varian_Bounds(3) - Var_out_sample_residuals(i,3);
                Varian_Mat(i, col_counter + 25) = normalized_difference_lower;
                Varian_Mat(i, col_counter + 26) = NaN;
                Varian_Mat(i, col_counter + 27) = NaN;
                Varian_Mat(i, col_counter + 28) = NaN;
                Varian_Mat(i, col_counter + 29) = NaN;
            else
                % it is not exact
                Varian_Mat(i, col_counter + 24) = NaN;
                Varian_Mat(i, col_counter + 25) = NaN;
                Varian_Mat(i, col_counter + 26) = Varian_Bounds(3) - Var_out_sample_residuals(i,3);
                Varian_Mat(i, col_counter + 27) = Varian_Bounds(4) - Var_out_sample_residuals(i,4);
                Varian_Mat(i, col_counter + 28) = normalized_difference_lower;
                Varian_Mat(i, col_counter + 29) = normalized_difference_upper;
            end

            % AVGSSQ difference & AVGSSQ normalized difference
            % (we added this check, because the whole point of the normalized difference is that it is either positive  or zero, but due to calculation issues it sometimes resulted stuff like "2.77555756156289E-16")   
            normalized_difference_lower = Varian_Bounds(5) - VARIAN_AVGSSQ_Lower_norm;
            normalized_difference_upper = Varian_Bounds(6) - VARIAN_AVGSSQ_Upper_norm;
            if abs(normalized_difference_lower) < 10^(-15)
                normalized_difference_lower = 0;
            end
            if abs(normalized_difference_upper) < 10^(-15)
                normalized_difference_upper = 0;
            end
            if Varian_Bounds(5) == Varian_Bounds(6) && Var_out_sample_residuals(i,5) == Var_out_sample_residuals(i,6)
                % it is exact
                Varian_Mat(i, col_counter + 30) = Varian_Bounds(5) - Var_out_sample_residuals(i,5);
                Varian_Mat(i, col_counter + 31) = normalized_difference_lower;
                Varian_Mat(i, col_counter + 32) = NaN;
                Varian_Mat(i, col_counter + 33) = NaN;
                Varian_Mat(i, col_counter + 34) = NaN;
                Varian_Mat(i, col_counter + 35) = NaN;
            else
                % it is not exact
                Varian_Mat(i, col_counter + 30) = NaN;
                Varian_Mat(i, col_counter + 31) = NaN;
                Varian_Mat(i, col_counter + 32) = Varian_Bounds(5) - Var_out_sample_residuals(i,5);
                Varian_Mat(i, col_counter + 33) = Varian_Bounds(6) - Var_out_sample_residuals(i,6);
                Varian_Mat(i, col_counter + 34) = normalized_difference_lower;
                Varian_Mat(i, col_counter + 35) = normalized_difference_upper;
            end

        end 
        % update the counter
        col_counter = col_counter + 36; %#ok<NASGU>

    end

end

end