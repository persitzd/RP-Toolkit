function [AFRIAT, Afriat_Mat] = HPZ_Afriat_Manager (AFRIAT_flags, expenditure, index_threshold, SDRP, residuals_waitbar, current_run, total_runs, subject_ID, varargin)

% initialization to avoid errors, just in case
AFRIAT = NaN; %#ok<NASGU>
Afriat_Mat = [];


% number of observations
[obs_num , ~] = size(expenditure);

% we calculate SDRP with transitivity, instead of DRP with transitivity,
% since in the case of Varian weak relations are of no interest
[~, transitive_SDRP, ~] = GARP_based_on_DRP_and_SDRP(SDRP, SDRP);

% divide the observations to smallest possible subsets such that
% observations from different subsets never belong to one cycle
% (a cycle of strictly revealed relations in this case)
[subsets_of_observations, num_of_subsets] = HPZ_Indices_Problem_Distribute(transitive_SDRP, transitive_SDRP);

% initialization
AFRIAT_per_subset = nan(1, num_of_subsets);

for i=1:num_of_subsets

    current_subset = subsets_of_observations{i};

    if length(current_subset) == 1
        
        % insert to the vector of Afriat per subset
        AFRIAT_per_subset(i) = 0;
        
    else   % length(current_subset) > 1

        % extract the relevant expenditures
        subset_expenditure = expenditure(current_subset , current_subset);

        % calculate Afriat index
        AFRIAT_subset = HPZ_Afriat_efficiency_index (subset_expenditure, index_threshold);

        % insert to the vector of Afriat per subset
        AFRIAT_per_subset(i) = AFRIAT_subset;
    end

end

AFRIAT = max(AFRIAT_per_subset);
% we add "eps", since we assume that if we got to Varian calculation, then
% the subject does not satisfy GARP, and we want that in this case the
% Varian won't be 0, but will be a very very low positive number
AFRIAT = AFRIAT + eps;


% AFRIAT residuals
if AFRIAT_flags(2)

    % out of sample (the only option actually for Afriat residuals)
    if AFRIAT_flags(4)
        
        % define the waitbar
        if (residuals_waitbar)
            waitbar_name = char(strcat(HPZ_Constants.waitbar_name_estimation, {' '}, '(', HPZ_Constants.current_run_waitbar, {' '}, num2str(current_run), {' '}, HPZ_Constants.total_runs_waitbar, {' '}, num2str(total_runs),')'));
            waitbar_msg = char(strcat(HPZ_Constants.waitbar_recovery, {' '}, num2str(subject_ID), {' '}, HPZ_Constants.waitbar_residuals_AFRIAT));
            new_bar_val = 0;
            h_wb = wide_waitbar(new_bar_val, {waitbar_msg, ''}, waitbar_name, HPZ_Constants.waitbar_width_multiplier, [0,0.12]);
        end
        
        % initialization of residuals matrix
        num_of_columns = 3;
        Afriat_Mat = nan(obs_num, num_of_columns);

        % GARP_vector is a vector that has 0 in the i'th place, iff
        % after dropping the i'th observation, the subject satisfies GARP  
        if ~isempty(varargin)
            GARP_vector = varargin{1};
        else
            GARP_vector = ones(obs_num, 1);
        end
        
        % initialization of column counter
        col_counter = 1;
        
        for i=1:num_of_subsets

            current_subset = subsets_of_observations{i};

            if length(current_subset) == 1
                
                Afriat_Mat(current_subset(1), col_counter) = AFRIAT;               % full index
                Afriat_Mat(current_subset(1), col_counter + 1) = AFRIAT;           % partial index
                Afriat_Mat(current_subset(1), col_counter + 2) = AFRIAT - AFRIAT;  % difference = 0
                
                % updating the waitbar
                if (residuals_waitbar)
                    new_bar_val = new_bar_val + 1/obs_num;
                    waitbar(new_bar_val, h_wb, {waitbar_msg , char(strcat({'Completed '}, num2str(new_bar_val*obs_num), {' observations out of '}, num2str(obs_num)))});
                end
                
            else   % length(current_subset) > 1

                % extract the relevant expenditures
                subset_expenditure = expenditure(current_subset , current_subset);
                subset_GARP_vector = GARP_vector(current_subset);
                
                % Afriat of the rest of the observations without this subset  
                AFRIAT_other_subsets = max(AFRIAT_per_subset([1:(i-1),(i+1):end]));
                
                % we calculate for each observation its residual
                for j=1:length(current_subset)
                    
                    if subset_GARP_vector(j) == 0
                        % perfectly consistent
                        AFRIAT_partial_full = 0;
                    else
                        % not perfectly consistent
                        truncated_expenditure = subset_expenditure([1:(j-1) , (j+1):end] , [1:(j-1) , (j+1):end]);
                        AFRIAT_partial = HPZ_Afriat_efficiency_index (truncated_expenditure, index_threshold);
                        AFRIAT_partial_full = max(AFRIAT_partial , AFRIAT_other_subsets);
                        % we add "eps", since we assume that if we got to Varian calculation, then
                        % the subject does not satisfy GARP, and we want that in this case the
                        % Varian won't be 0, but will be a very very low positive number
                        AFRIAT_partial_full = AFRIAT_partial_full + eps;
                    end
                    Afriat_Mat(current_subset(j), col_counter) = AFRIAT;                            % full index
                    Afriat_Mat(current_subset(j), col_counter + 1) = AFRIAT_partial_full;           % partial index
                    Afriat_Mat(current_subset(j), col_counter + 2) = AFRIAT - AFRIAT_partial_full;  % difference
                
                    % updating the waitbar
                    if (residuals_waitbar)
                        new_bar_val = new_bar_val + 1/obs_num;
                        waitbar(new_bar_val, h_wb, {waitbar_msg , char(strcat({'Completed '}, num2str(new_bar_val*obs_num), {' observations out of '}, num2str(obs_num)))});
                    end
                end 
                
            end

        end

        % update the counter (not needed, was needed in older code where other indexes residuals were calculated in the same function)  
        col_counter = col_counter + 3; %#ok<NASGU>
        
        % close the waitbar
        if (residuals_waitbar)
            close(h_wb);
        end
    end

end


end