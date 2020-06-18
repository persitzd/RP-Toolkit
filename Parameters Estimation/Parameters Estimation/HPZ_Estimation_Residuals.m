function Mat = HPZ_Estimation_Residuals (action_flag, data_matrix, obs_num, choice_set_type, treatment, function_flag, param1_restrictions, param2_restrictions, fix_corners, metric_flag, aggregation_flag, asymmetric_flag, in_sample_flag, out_sample_flag, param_1, param_2, main_criterion, pref_class, numeric_flag, max_time_estimation, min_counter, max_starting_points, BI_threshold, debugger_mode, active_waitbar, current_run, total_runs)

% this function performs residual calculations for the estimation of a
% single subject. it may calculate in-sample residuals, out-of-sample
% residuals, or both, depending on what it is set to do.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 

% the function returns a matrix of results, as follow:
%   column 1: subject number
%   column 2: 1st parameter value 
%   column 3: 2nd parameter value
%   column 4: original criterion value 
%   column 5: observation number 
% if in-sample is required, the next column is the in-sample residual for
% this observation.
% if out-of-sample is required, the next 3 columns will be:
%   1: alternative 1st parameter value (when estimating without ths observation) 
%   2: alternative 2nd parameter value (when estimating without ths observation) 
%   3: alternative criterion (when estimating without ths observation) 
%   4: out-of-sample residual for this observation 



if (fix_corners)
    % Choi et al. (2007) correction is applied for corner choices
    subject_data = HPZ_No_Corners (data_matrix, obs_num, 1);
else
    % no correction should be applied for corner choices
    subject_data = data_matrix;
end


% determining how many columns are required in the result matrix
number_of_columns = 5;
if in_sample_flag == true
    number_of_columns = number_of_columns + 2;
end
if out_sample_flag == true
    number_of_columns = number_of_columns + 4;
end
% creating the result matrix
Mat = nan (obs_num, number_of_columns);

% the given parameters
param = [param_1 , param_2];

% basic data - subject number, parameters values, 
% criterion value and observation number 
Mat(:,1) = repelem(subject_data(1,1), obs_num);
Mat(:,2) = repelem(param(1), obs_num);
Mat(:,3) = repelem(param(2), obs_num);
Mat(:,4) = repelem(main_criterion, obs_num);
Mat(:,5) = 1:obs_num;

% current (next) column to be printed
current_col = 6;

if choice_set_type ~= HPZ_Constants.choice_set_finite_set
    Choices(:, 1:4) = subject_data(1:obs_num, 3:6);
    expenditure = (Choices(:,1)*Choices(:,3)' + Choices(:,2)*Choices(:,4)')';
    endowments = diag(expenditure);
else % i.e. choice_set_type == HPZ_Constants.choice_set_finite_set
    Choices = subject_data(1:obs_num, 3:end);
    expenditure = nan; %#ok<NASGU>
    endowments = nan;
end


% define the waitbar
if (active_waitbar)
    waitbar_name = char(strcat(HPZ_Constants.waitbar_name_estimation, {' '}, '(', HPZ_Constants.current_run_waitbar, {' '}, num2str(current_run), {' '}, HPZ_Constants.total_runs_waitbar, {' '}, num2str(total_runs),')'));
    waitbar_msg = char(strcat(HPZ_Constants.waitbar_recovery, {' '}, num2str(data_matrix(1,1)), {' '}, HPZ_Constants.waitbar_residuals));
    new_bar_val = 0;
    h_wb = wide_waitbar(new_bar_val, {waitbar_msg, ''}, waitbar_name, HPZ_Constants.waitbar_width_multiplier, [0,0.12]);
end


%% In Sample
if (in_sample_flag)
    
    % (in-sample is very short, so we don't bother updating the waitbar here) 
    
    % NLLS
    if action_flag == HPZ_Constants.NLLS_action
        
        % the optimal cohices for the subject given these parameters
        if choice_set_type ~= HPZ_Constants.choice_set_finite_set
            if numeric_flag == HPZ_Constants.numeric
                [predicted_choices, ~] = HPZ_NLLS_Choices_Numeric (param, Choices(:,3:4), endowments, treatment, function_flag, asymmetric_flag, pref_class, debugger_mode);
            elseif numeric_flag == HPZ_Constants.analytic
                [predicted_choices, ~] = HPZ_NLLS_Choices_Analytic(param, Choices(:,1:4), function_flag, pref_class, debugger_mode);
            end
        else
            [predicted_choices, ~] = HPZ_NLLS_Choices_Finite_Set(param, Choices(:,1:end), function_flag, pref_class, debugger_mode);
        end
        % (there is no: "numeric_flag == HPZ_Constants.semi_numeric"
        % since semi-numeric is unique for MMI and BI)
        
        %observed_choices = Choices(:,1:2);
        
        % calculating the NLLS in-sample component residuals
        in_sample_component = HPZ_NLLS_Criterion_Per_Observation(Choices, predicted_choices, metric_flag);
        Mat (:, current_col) = in_sample_component;
        
        % taking the desired metric   ( (x,y) = (observations, optimal_bundles) )
        if metric_flag == HPZ_Constants.euclidean_metric
            metric_function = @(x,y) sum(HPZ_NLLS_Criterion_Euclid(x,y));
        elseif metric_flag == HPZ_Constants.CFGK_metric
            metric_function = @(x,y) sum(HPZ_NLLS_Criterion_Ldr(x,y));
        elseif metric_flag == HPZ_Constants.normalized_euclidean_metric
            metric_function = @(x,y) mean(HPZ_NLLS_Criterion_Euclid_normalized(x,y));
        end
        
        % calculate in-sample difference residuals
        for i=1:obs_num
            Mat (i, current_col+1) = main_criterion - metric_function(Choices([1:(i-1),(i+1):end], 1:2) , predicted_choices([1:(i-1),(i+1):end], 1:2));
        end
        
        
    % MMI
    elseif action_flag == HPZ_Constants.MMI_action
        
        % calculating the MMI in-sample component residuals
        in_sample_component = HPZ_MMI_Criterion_Per_Observation (param, endowments, Choices, treatment, function_flag, pref_class, numeric_flag, debugger_mode);
        Mat (:, current_col) = in_sample_component;
        
        % taking the desired aggregator
        if aggregation_flag == HPZ_Constants.MMI_Max
            aggregator = @max;
        elseif aggregation_flag == HPZ_Constants.MMI_Mean
            aggregator = @mean;
        elseif aggregation_flag == HPZ_Constants.MMI_AVGSSQ
            aggregator = @(x) sqrt(meansqr(x));
        end
        
        % calculate in-sample difference residuals
        for i=1:obs_num
            Mat (i, current_col+1) = main_criterion - aggregator(in_sample_component([1:(i-1),(i+1):end]));
        end
        
        
    % BI
    elseif action_flag == HPZ_Constants.BI_action
        
        % calculating the MMI in-sample component residuals
        MMI_residual_i = HPZ_MMI_Criterion_Per_Observation (param, endowments, Choices, treatment, function_flag, pref_class, numeric_flag, debugger_mode);
        % if MMI criterion < BI_threshold then we will get 0,
        % otherwise we will get 1
        in_sample_component = ceil(MMI_residual_i - BI_threshold);
        Mat (:, current_col) = in_sample_component;
        
        % calculate in-sample difference residuals
        for i=1:obs_num
            Mat (i, current_col+1) = main_criterion - mean(in_sample_component([1:(i-1),(i+1):end]));
        end
    end
    
end

% updating the index of the next (current) column
current_col = current_col + 2;



%% Out of Sample
if (out_sample_flag) 
    
    for i=1:obs_num
        
        if i < obs_num
            adjusted_data = [subject_data(1:(i-1),:) ; subject_data((i+1):obs_num,:)];
        elseif i == obs_num
            adjusted_data = subject_data(1:(obs_num-1),:);
        end
        
        % we first perform an estimation on the truncated data
        [param_i , criterion_i , ~] = HPZ_Estimation(adjusted_data, (obs_num-1), choice_set_type, action_flag, treatment, function_flag, param1_restrictions, param2_restrictions, ...
                                    fix_corners, metric_flag, asymmetric_flag, aggregation_flag, pref_class, numeric_flag, false, ...
                                    max_time_estimation, min_counter, max_starting_points, BI_threshold, debugger_mode, false, current_run, total_runs);
        
        % we save the resulting parameters and the resulting criterion 
        % to the result matrix 
        Mat(i, current_col)     = param_i(1);
        Mat(i, current_col+1)   = param_i(2);
        Mat(i, current_col+2)   = criterion_i;
        Mat(i, current_col+3)   = main_criterion - criterion_i;
        
        % updating the waitbar
        if (active_waitbar)
            new_bar_val = i / obs_num;
            waitbar(new_bar_val, h_wb, {waitbar_msg , char(strcat({'Completed '}, num2str(i), {' observations out of '}, num2str(obs_num)))});
        end
        
    end   % end of loop
    
end


% close the waitbar
if (active_waitbar)
    close(h_wb);
end


end

