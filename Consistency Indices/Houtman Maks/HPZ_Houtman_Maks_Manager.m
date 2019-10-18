function [HM, HM_residuals, HM_raw, HM_raw_residuals, HM_Mat] = HPZ_Houtman_Maks_Manager (HOUTMAN_flags, DRP, SDRP, RP, identical_choice)

% number of observations
[obs_num , ~] = size(DRP);

% for i=1:obs_num
%     DRP(i,i)=0;
% end

% initialization to avoid error when residuals are not required
HM_Mat = [];


% whether out-of-sample is required by the user
out_of_sample_required = (HOUTMAN_flags(2) && HOUTMAN_flags(4));


% divide the observations to smallest possible subsets such that
% observations from different subsets never belong to one cycle that
% violates the relevant-GARP
[subsets_of_observations, num_of_subsets] = HPZ_Indices_Problem_Distribute(RP, RP);

% initialization of how many observations should we drop
HM_raw_accumulated = 0;
% initialization of HM raw per subset
HM_raw_per_subset = nan(num_of_subsets, 1);

% initialization (it must be initialized to 0)
% raw_residuals are 0 if deleting this observation will not improve the HM index,  
% and are 1 if deleting this observation will not improve the HM index.
HM_raw_residuals = nan(obs_num, 1);

% now we calculate the HM separately for each subset
for i=1:num_of_subsets

    if length(subsets_of_observations{i}) > 1

        %subset = subsets_of_observations{i}

        % create new DRP & SDRP matrices with only the observations of this subset 
        subset_DRP = DRP(subsets_of_observations{i} , subsets_of_observations{i});
        subset_SDRP = SDRP(subsets_of_observations{i} , subsets_of_observations{i});
        subset_identical_choice = identical_choice(subsets_of_observations{i} , subsets_of_observations{i});
        
        % calculate the HM index for this subset of choices
        [~, subset_HM_sum, subset_HM_raw_residuals] = HPZ_Houtman_Maks_Index (subset_DRP, subset_SDRP, subset_identical_choice, out_of_sample_required);
        
        % update the accumulated HM
        HM_raw_accumulated = HM_raw_accumulated + subset_HM_sum;
        % update / assign to the subset HM
        HM_raw_per_subset(i) = subset_HM_sum;

        % update the residuals for these observations
        % if the subset residual is bigger than the subset HM, it means
        % that the raw residual was 0, and if it is smaller, it means that
        % the raw residual was 1
        HM_raw_residuals(subsets_of_observations{i}) = subset_HM_raw_residuals; 
        
    else
        
        % we need to assign a raw residual of 0 to this single observation  
        HM_raw_residuals(subsets_of_observations{i}) = 0;
        
    end

end

% the final HM index
HM_raw = HM_raw_accumulated;
HM = HM_raw_accumulated / obs_num;
% finalizing the residuals
HM_residuals = (HM_raw - HM_raw_residuals) / (obs_num - 1);





%% HOUTMAN-MAKS residuals - assignment to matrix
if HOUTMAN_flags(2) && HOUTMAN_flags(4)
    % out of sample (the only option actually for residuals)
        
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