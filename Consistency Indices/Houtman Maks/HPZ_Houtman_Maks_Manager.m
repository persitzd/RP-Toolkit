function [HM, HM_residuals, HM_raw, HM_raw_residuals, HM_exact] = HPZ_Houtman_Maks_Manager (residuals_flag, DRP, SDRP, RP, SRP, is_2_goods)

% number of observations
[obs_num , ~] = size(DRP);


% thanks to the Rose (1958) Theorem, we can use a very efficient algorithm
% for minimum vertex cover when there are only 2 goods
if is_2_goods
    
    % this algorithm is always exact
    HM_exact = 1;
    
    % a matrix that represents the graph
    GARP_HM = DRP .* SDRP';
    GARP_couples = max(GARP_HM, GARP_HM');
    
    if ~residuals_flag
        % we just calculate the index and that's it
        [HM, HM_raw, ~] = HPZ_Houtman_Maks_Graph_Approach (GARP_couples, 0);
        % these initializations are just to avoid bugs/crashs
        HM_raw_residuals = zeros(obs_num, 1);
        HM_residuals = zeros(obs_num, 1);
    else
        % divide the observations to smallest possible subsets such that
        % observations from different subsets never belong to one 2-length cycle 
        % (this allows the out-of-sample to work faster)
        [subsets_of_observations, num_of_subsets] = HPZ_Indices_Problem_Distribute(DRP, SDRP);
        
        % initializations
        HM_raw_accumulated = 0;
        HM_raw_per_subset = zeros(num_of_subsets, 1);
        HM_raw_pseudo_residuals = zeros(obs_num, 1);
        HM_raw_residuals = zeros(obs_num, 1);
        
        % now we calculate the Houtman-Maks separately for each subset
        for i=1:num_of_subsets

            if length(subsets_of_observations{i}) == 1
                % no need to drop anything
                HM_raw_per_subset(i) = 0;
            else
                % extract the relevant part of the graph
                subset_GARP_couples = GARP_couples(subsets_of_observations{i} , subsets_of_observations{i});

                % calculate the Houtman-Maks index for this subset of choices
                [~, subset_HM_raw, subset_HM_raw_pseudo_residuals] = HPZ_Houtman_Maks_Graph_Approach (subset_GARP_couples, 0);
                HM_raw_per_subset(i) = subset_HM_raw;
                HM_raw_pseudo_residuals(subsets_of_observations{i}) = subset_HM_raw_pseudo_residuals;
                % update the accumulated HM
                HM_raw_accumulated = HM_raw_accumulated + subset_HM_raw;
            end
        end
        
        % the final HM index
        HM_raw = HM_raw_accumulated;
        HM = HM_raw_accumulated / obs_num;

        
        
        % now for the residuals themselves
        for i=1:num_of_subsets

            current_subset = subsets_of_observations{i};
            
            if length(current_subset) == 1
                % no need to drop anything
                HM_raw_residuals(current_subset(1)) = 0;
            else
                % extract the relevant part of the graph
                subset_GARP_couples = GARP_couples(current_subset , current_subset);
                
                for j=1:length(current_subset)
                    
                    if HM_raw_pseudo_residuals(current_subset(j)) == 1
                        % we know it is 1
                        HM_raw_residuals(current_subset(j)) = 1;
                    else
                        % we need to calculate it, to check if it is 0 or 1 
                        % creating truncated matrices
                        remaining_indices = [1:(j-1) , (j+1):length(current_subset)];
                        truncated_GARP_couples = subset_GARP_couples(remaining_indices , remaining_indices);
                        [~, out_sample_HM_raw, ~] = HPZ_Houtman_Maks_Graph_Approach (truncated_GARP_couples, 0);
                        HM_raw_residuals(current_subset(j)) = HM_raw_per_subset(i) - out_sample_HM_raw;
                    end
                end
            end
        end
        
        % finalizing the residuals
        HM_residuals = (HM_raw - HM_raw_residuals) / (obs_num - 1);
    end
    
    return
end





%% all the remaining of the function is for the case of more than 2 goods

% divide the observations to smallest possible subsets such that
% observations from different subsets never belong to one cycle that
% violates the relevant-GARP
[subsets_of_observations, num_of_subsets] = HPZ_Indices_Problem_Distribute(RP, SRP);





% initialization of how many observations should we drop
HM_raw_accumulated = 0;

% initialization (to: 1=exact)
HM_exact = 1;

% initialization (it must be initialized to 0)
% raw_residuals are 0 if deleting this observation will not improve the HM index,  
% and are 1 if deleting this observation will not improve the HM index.
HM_raw_residuals = zeros(obs_num, 1);

% now we calculate the HM separately for each subset
for i=1:num_of_subsets
    
    if length(subsets_of_observations{i}) > 1
        
        %subset = subsets_of_observations{i}
        
        % create new DRP & SDRP matrices with only the observations of this subset 
        subset_DRP = DRP(subsets_of_observations{i} , subsets_of_observations{i});
        subset_SDRP = SDRP(subsets_of_observations{i} , subsets_of_observations{i});
        
        % calculate the HM index for this subset of choices
        [~, ~, subset_HM_raw, subset_HM_raw_residuals, subset_HM_exact] = HPZ_Houtman_Maks_efficiency_index (residuals_flag, subset_DRP, subset_SDRP, is_2_goods);

        % (Next is a work still in progress, that will replace the above function in future versions) 
        %%[subset_GARP, ~] = GARP_based_on_DRP_and_SDRP(subset_DRP, subset_SDRP);
        %[~, subset_HM_raw, subset_HM_raw_pseudo_residuals, subset_HM_exact] = HPZ_Houtman_Maks_Cycles_Approach (residuals_flag, [], 1, subset_DRP, subset_SDRP);
        %subset_HM_raw_residuals = subset_HM_raw_pseudo_residuals;
        
        % update the accumulated HM
        HM_raw_accumulated = HM_raw_accumulated + subset_HM_raw;

        % update the is-exact
        % NOTE! we use min because the smaller the value of this variable,
        % the LESS exact it is. In Varian for instance it is the opposite. 
        HM_exact = min(HM_exact, subset_HM_exact);

        % update the residuals for these observations
        % if the subset residual is bigger than the subset HM, it means
        % that the raw residual was 0, and if it is smaller, it means that
        % the raw residual was 1
        HM_raw_residuals(subsets_of_observations{i}) = subset_HM_raw_residuals;
        
    end
    
end

% the final HM index
HM_raw = HM_raw_accumulated;
HM = HM_raw_accumulated / obs_num;

% finalizing the residuals
HM_residuals = (HM_raw - HM_raw_residuals) / (obs_num - 1);



end