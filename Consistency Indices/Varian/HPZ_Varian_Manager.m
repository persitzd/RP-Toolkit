function [Varian, Varian_exact, Var_in_sample_residuals, Var_out_sample_residuals, one_minus_v] = HPZ_Varian_Manager (expenditure, identical_choice, index_threshold, SRP, out_of_sample)


% number of observations
[obs_num , ~] = size(expenditure);



%[~, SRP] = GARP_based_on_DRP_and_SDRP(SDRP, SDRP);

% calculate GARP with transitivity (in Varian, weak relations are of no
% interest, since a v{i} of epsilon->0 is enough to break them)
relevant_RP = SRP;
relevant_SRP = SRP;

% divide the observations to smallest possible subsets such that
% observations from different subsets never belong to one cycle that
% violates the relevant-GARP
[subsets_of_observations, num_of_subsets] = HPZ_Indices_Problem_Distribute(relevant_RP, relevant_SRP);



% initialization (to: 1=exact)
Varian_exact = 1;

% we need to keep is-exact for each subset separately, for out-of-sample residuals
Varian_exact_per_subset = ones(1, num_of_subsets);

% initialization (it must be initialized to 0)
one_minus_v_avg = zeros(obs_num, 1);
one_minus_v_meanssq = zeros(obs_num, 1);
one_minus_v_min = zeros(obs_num, 1);

% now we calculate the Varian separately for each subset
for i=1:num_of_subsets
    
    current_subset = subsets_of_observations{i};
    
    if length(current_subset) > 1
        
        % extract the relevant expenditures and identical choices
        subset_expenditure = expenditure(current_subset , current_subset);
        subset_identical_choice = identical_choice(current_subset , current_subset);
        
        % calculate the Varian index for this subset of choices
        [~, subset_Varian_exact, ~, subset_one_minus_v] = HPZ_Varian_efficiency_index (subset_expenditure, subset_identical_choice, index_threshold);
        
        % update 1-v
        one_minus_v_min(current_subset)     = subset_one_minus_v(:,1);
        one_minus_v_avg(current_subset)     = subset_one_minus_v(:,2);
        one_minus_v_meanssq(current_subset) = subset_one_minus_v(:,3);
        
        % update the is-exact
        % NOTE! we use max because the bigger the value of this variable,
        % the LESS exact it is. In Houtman-Maks for instance it is the opposite. 
        Varian_exact = max(Varian_exact, subset_Varian_exact);
        
        % these are for out-of-sample
        Varian_exact_per_subset(i) = subset_Varian_exact;
    end
    
end



% the final Varian indices
min_var = max(one_minus_v_min);
average_var = mean(one_minus_v_avg);
meanssq_var = sqrt(meansqr(one_minus_v_meanssq));
Varian = [min_var, average_var , meanssq_var];

% finalizing the in-sample residuals
min_var_in_sample = HPZ_Consistency_Indices_In_Sample_Residuals_Calc (one_minus_v_min, @max);
average_var_in_sample = HPZ_Consistency_Indices_In_Sample_Residuals_Calc (one_minus_v_avg, @mean);
meanssq_var_in_sample = HPZ_Consistency_Indices_In_Sample_Residuals_Calc (one_minus_v_meanssq, @(x) sqrt(meansqr(x)));
Var_in_sample_residuals = [min_var_in_sample , average_var_in_sample , meanssq_var_in_sample];

% the best 1-v for each of the aggregators
one_minus_v = [one_minus_v_min, one_minus_v_avg, one_minus_v_meanssq];



% initializations
min_var_out_sample = zeros(obs_num, 1);
average_var_out_sample = zeros(obs_num, 1);
meanssq_var_out_sample = zeros(obs_num, 1);
is_exact_var_out_sample = zeros(obs_num, 1);

% if out-of-sample residuals are required - calculate them
if out_of_sample
    
    for i=1:num_of_subsets
        
        current_subset = subsets_of_observations{i};
        
        if length(current_subset) == 1
            % it does not involve in any relevant cycle -
            % it does not effect the varian index
            % we only need to multiply by n/(n-1) when it is an average
            is_exact_var_out_sample(current_subset(1)) = Varian_exact;
            min_var_out_sample(current_subset(1)) = min_var;
            average_var_out_sample(current_subset(1)) = average_var*obs_num/(obs_num-1);
            meanssq_var_out_sample(current_subset(1)) = sqrt( (meanssq_var^2)*obs_num/(obs_num-1) );
        else
            % we need to calculate the out-of-sample for each of these
            % observations *inside* the subset, then add it to what's
            % outside this subset
            % PART A - what's outside
            indices_outside_subset = ones(1, obs_num);
            indices_outside_subset(current_subset) = 0;
            indices_outside_subset = find(indices_outside_subset);
            % PART B - what's inside (+ combining inside & outside)
            subset_expenditure = expenditure(current_subset , current_subset);
            subset_identical_choice = identical_choice(current_subset , current_subset);
            subset_SRP = SRP(current_subset , current_subset);
            for j=1:length(current_subset)
                % creating truncated matrices
                remaining_indices = [1:(j-1) , (j+1):length(current_subset)];
                truncated_expenditure = subset_expenditure(remaining_indices , remaining_indices);
                truncated_identical_choice = subset_identical_choice(remaining_indices , remaining_indices);
                truncated_SRP = subset_SRP(remaining_indices , remaining_indices);
                % calculating varian for the truncated data
                [~, out_sample_Var_exact, ~, ~, out_sample_one_minus_v] = HPZ_Varian_Manager (truncated_expenditure, truncated_identical_choice, index_threshold, truncated_SRP, false);
                
                % now we assign the results
                is_exact_all_but_this_subset = max([Varian_exact_per_subset(1:(i-1)), Varian_exact_per_subset((i+1):end)]);
                is_exact_var_out_sample(current_subset(j)) = max(is_exact_all_but_this_subset, out_sample_Var_exact);
                min_var_out_sample(current_subset(j)) = max([out_sample_one_minus_v(:,1)' , one_minus_v_min(indices_outside_subset)']);
                average_var_out_sample(current_subset(j)) = mean([out_sample_one_minus_v(:,2)' , one_minus_v_avg(indices_outside_subset)']);
                meanssq_var_out_sample(current_subset(j)) = sqrt(meansqr([out_sample_one_minus_v(:,3)' , one_minus_v_meanssq(indices_outside_subset)']));
                %min_var_out_sample(current_subset(j)) = max(remaining_min_var, out_sample_Var(1));
                %average_var_out_sample(current_subset(j)) = ( out_sample_Var(2)*(length(current_subset) - 1) + remaining_average_var*(obs_num - length(current_subset)) ) / (obs_num - 1);
                %meanssq_var_out_sample(current_subset(j)) = sqrt( ( (out_sample_Var(3)^2)*(length(current_subset) - 1) + (remaining_meanssq_var^2)*(obs_num - length(current_subset)) ) / (obs_num - 1) );
            end
        end
        
    end
    
end

% finalizing the out-of-sample residuals
Var_out_sample_residuals = [min_var_out_sample , average_var_out_sample , meanssq_var_out_sample , is_exact_var_out_sample];



% avg = one_minus_v_avg'
% meanssq = one_minus_v_meanssq'
% minimum = one_minus_v_min'

end