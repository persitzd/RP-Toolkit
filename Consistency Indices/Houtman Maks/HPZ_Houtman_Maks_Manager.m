function [HM, HM_residuals, HM_raw, HM_raw_residuals, HM_Mat] = HPZ_Houtman_Maks_Manager (HOUTMAN_flags, DRP, SDRP, RP, identical_choice, num_of_goods, residuals_waitbar, current_run, total_runs, subject_ID)

% number of observations
[obs_num , ~] = size(DRP);

% for i=1:obs_num
%     DRP(i,i)=0;
% end



    
if num_of_goods == 2
    % if there are only 2 goods
    
    % thanks to the Rose (1958) Theorem, we can use an algorithm
    % for minimum vertex cover when there are only 2 goods
    
    
    
    % a matrix that represents the graph
    GARP_HM = DRP .* SDRP';
    % for SARP (instead of GARP) use this line. if you change to SARP, change the SDRP to DRP in HPZ_Indices_Problem_Distribute ahead as well. 
    % also, the DRP used for SARP must be one that does not include identical choices. 
    %GARP_HM = DRP .* DRP';  
    GARP_couples = max(GARP_HM, GARP_HM');
    
    if ~(HOUTMAN_flags(2) && HOUTMAN_flags(4))
        % we just calculate the index and that's it
        [HM, HM_raw, ~] = HPZ_Houtman_Maks_Two_Goods_Only (GARP_couples, 0);
        % these initializations are just to avoid bugs/crashs
        HM_raw_residuals = zeros(obs_num, 1);
        HM_residuals = zeros(obs_num, 1);
        
    else
        
        % define the waitbar
        if (residuals_waitbar)
            waitbar_name = char(strcat(HPZ_Constants.waitbar_name_estimation, {' '}, '(', HPZ_Constants.current_run_waitbar, {' '}, num2str(current_run), {' '}, HPZ_Constants.total_runs_waitbar, {' '}, num2str(total_runs),')'));
            waitbar_msg = char(strcat(HPZ_Constants.waitbar_recovery, {' '}, num2str(subject_ID), {' '}, HPZ_Constants.waitbar_residuals_HOUTMAN));
            new_bar_val = 0;
            h_wb = wide_waitbar(new_bar_val, {waitbar_msg, ''}, waitbar_name, HPZ_Constants.waitbar_width_multiplier, [0,0.12]);
        end
        
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
                [~, subset_HM_raw, subset_HM_raw_pseudo_residuals] = HPZ_Houtman_Maks_Two_Goods_Only (subset_GARP_couples, 0);
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
                
                % updating the waitbar
                if (residuals_waitbar)
                    new_bar_val = new_bar_val + 1/obs_num;
                    waitbar(new_bar_val, h_wb, {waitbar_msg , char(strcat({'Completed '}, num2str(new_bar_val*obs_num), {' observations out of '}, num2str(obs_num)))});
                end
                
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
                        [~, out_sample_HM_raw, ~] = HPZ_Houtman_Maks_Two_Goods_Only (truncated_GARP_couples, 0);
                        HM_raw_residuals(current_subset(j)) = HM_raw_per_subset(i) - out_sample_HM_raw;
                    end
                    
                    % updating the waitbar
                    if (residuals_waitbar)
                        new_bar_val = new_bar_val + 1/obs_num;
                        waitbar(new_bar_val, h_wb, {waitbar_msg , char(strcat({'Completed '}, num2str(new_bar_val*obs_num), {' observations out of '}, num2str(obs_num)))});
                    end
                    
                end
            end
        end
        
        % finalizing the residuals
        HM_residuals = (HM_raw - HM_raw_residuals) / (obs_num - 1);
        
        % close the waitbar
        if (residuals_waitbar)
            close(h_wb);
        end
        
    end


    
else
    % if there are more than 2 goods
    
    % when there are more than 2 goods, without the Rose (1958) Theorem, 
    % we use an algorithm for minimum vertex cover in hypergraph 
    
    
    
    % divide the observations to smallest possible subsets such that
    % observations from different subsets never belong to one cycle that
    % violates the relevant-GARP
    [subsets_of_observations, num_of_subsets] = HPZ_Indices_Problem_Distribute(RP, RP);
    
    % initialization of how many observations should we drop
    HM_raw_accumulated = 0;
    % initialization of HM raw per subset
    HM_raw_per_subset = zeros(num_of_subsets, 1);
    
    % initialization (it must be initialized to 0)
    % raw_residuals are 0 if deleting this observation will not improve the HM index,  
    % and are 1 if deleting this observation will not improve the HM index.
    HM_raw_pseudo_residuals = zeros(obs_num, 1);

    % now we calculate the HM separately for each subset
    for i=1:num_of_subsets
        
        if length(subsets_of_observations{i}) > 1

            %subset = subsets_of_observations{i}

            % create new DRP & SDRP matrices with only the observations of this subset 
            subset_DRP = DRP(subsets_of_observations{i} , subsets_of_observations{i});
            subset_SDRP = SDRP(subsets_of_observations{i} , subsets_of_observations{i});
            subset_identical_choice = identical_choice(subsets_of_observations{i} , subsets_of_observations{i});
            
            % calculate the HM index for this subset of choices
            %[~, ~, subset_HM_raw, subset_HM_raw_residuals, subset_HM_exact] = HPZ_Houtman_Maks_efficiency_index (residuals_flag, subset_DRP, subset_SDRP, is_2_goods);

            % (Next is a work still in progress, that will replace the above function in future versions) 
            %[subset_GARP, ~, ~] = GARP_based_on_DRP_and_SDRP(subset_DRP, subset_SDRP);
            [~, subset_HM_raw, ~, subset_HM_raw_pseudo_residuals] = HPZ_Houtman_Maks_Any_Number_Of_Goods (subset_DRP, subset_SDRP, subset_identical_choice);
            %subset_weights = ones(1,length(subsets_of_observations{i}));
            %[~, subset_HM_raw, subset_HM_raw_pseudo_residuals, subset_HM_exact] = HPZ_Houtman_Maks_Weighted_Cycles_Approach (subset_SDRP, subset_weights);

            % update the accumulated HM
            HM_raw_accumulated = HM_raw_accumulated + subset_HM_raw;
            % update / assign to the subset HM
            HM_raw_per_subset(i) = subset_HM_raw;
            
            % update the residuals for these observations
            % if the subset residual is bigger than the subset HM, it means
            % that the raw residual was 0, and if it is smaller, it means that
            % the raw residual was 1
            HM_raw_pseudo_residuals(subsets_of_observations{i}) = subset_HM_raw_pseudo_residuals;

        end

    end

    % the final HM index
    HM_raw = HM_raw_accumulated;
    HM = HM_raw_accumulated / obs_num;

    % initializing the residuals (just to prevent program from crashing)
    HM_raw_residuals = 0;
    HM_residuals = 0;

end





% initialization to avoid error when residuals are not required
HM_Mat = [];

%% HOUTMAN-MAKS residuals - assignment to matrix
if HOUTMAN_flags(2)

    % out of sample (the only option actually for residuals)
    if HOUTMAN_flags(4)
        
        % define the waitbar
        if (residuals_waitbar)
            waitbar_name = char(strcat(HPZ_Constants.waitbar_name_estimation, {' '}, '(', HPZ_Constants.current_run_waitbar, {' '}, num2str(current_run), {' '}, HPZ_Constants.total_runs_waitbar, {' '}, num2str(total_runs),')'));
            waitbar_msg = char(strcat(HPZ_Constants.waitbar_recovery, {' '}, num2str(subject_ID), {' '}, HPZ_Constants.waitbar_residuals_HOUTMAN));
            new_bar_val = 0;
            h_wb = wide_waitbar(new_bar_val, {waitbar_msg, ''}, waitbar_name, HPZ_Constants.waitbar_width_multiplier, [0,0.12]);
        end
        
        % initialization
        HM_raw_residuals = zeros(obs_num, 1);
        
        % now for the residuals themselves
        for i=1:num_of_subsets

            current_subset = subsets_of_observations{i};
            
            if length(current_subset) == 1
                
                % no need to drop anything
                HM_raw_residuals(current_subset(1)) = 0;
                
                % updating the waitbar
                if (residuals_waitbar)
                    new_bar_val = new_bar_val + 1/obs_num;
                    waitbar(new_bar_val, h_wb, {waitbar_msg , char(strcat({'Completed '}, num2str(new_bar_val*obs_num), {' observations out of '}, num2str(obs_num)))});
                end
                
            else
                % create new DRP & SDRP matrices with only the observations of this subset 
                subset_DRP = DRP(subsets_of_observations{i} , subsets_of_observations{i});
                subset_SDRP = SDRP(subsets_of_observations{i} , subsets_of_observations{i});
                subset_identical_choice = identical_choice(subsets_of_observations{i} , subsets_of_observations{i});
                
                for j=1:length(current_subset)
                    
                    if HM_raw_pseudo_residuals(current_subset(j)) == 1
                        % we know it is 1
                        HM_raw_residuals(current_subset(j)) = 1;
                    else
                        % we need to calculate it, to check if it is 0 or 1 
                        % creating truncated matrices
                        remaining_indices = [1:(j-1) , (j+1):length(current_subset)];
                        truncated_DRP = subset_DRP(remaining_indices , remaining_indices);
                        truncated_SDRP = subset_SDRP(remaining_indices , remaining_indices);
                        truncated_identical_choice = subset_identical_choice(remaining_indices , remaining_indices);
                        [~, out_sample_HM_raw, ~, ~] = HPZ_Houtman_Maks_Any_Number_Of_Goods (truncated_DRP, truncated_SDRP, truncated_identical_choice);
                        HM_raw_residuals(current_subset(j)) = HM_raw_per_subset(i) - out_sample_HM_raw;
                    end
                    
                    % updating the waitbar
                    if (residuals_waitbar)
                        new_bar_val = new_bar_val + 1/obs_num;
                        waitbar(new_bar_val, h_wb, {waitbar_msg , char(strcat({'Completed '}, num2str(new_bar_val*obs_num), {' observations out of '}, num2str(obs_num)))});
                    end
                end
            end
        end
        
        % finalizing the residuals
        HM_residuals = (HM_raw - HM_raw_residuals) / (obs_num - 1);

        
        % close the waitbar
        if (residuals_waitbar)
            close(h_wb);
        end
        
        
        % initialization of residuals matrix
        num_of_columns = 4;
        HM_Mat = zeros(obs_num, num_of_columns);
        
        % initialization of column counter
        col_counter = 1;
        
        for i=1:obs_num
            HM_Mat(i, col_counter) = HM;                            % full index
            HM_Mat(i, col_counter + 1) = HM_residuals(i);           % partial index
            HM_Mat(i, col_counter + 2) = HM - HM_residuals(i);      % difference
            % (we added this check, because the whole point of the
            % normalized difference is that it is either positive
            % or zero, but due to calculation issues it sometimes
            % resulted stuff like "2.77555756156289E-17")
            normalized_difference = HM - HM_residuals(i)*(obs_num-1)/obs_num;
            if abs(normalized_difference) < 10^(-15)
                normalized_difference = 0;
            end
            HM_Mat(i, col_counter + 3) = normalized_difference;     % normalized difference
        end

        % update the counter
        col_counter = col_counter + 4; %#ok<NASGU>
    end

end


end