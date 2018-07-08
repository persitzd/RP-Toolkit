function [param, main_criterion, final_output] = HPZ_Estimation (data, obs_num, action_flag, treatment, function_flag, param1_restrictions, param2_restrictions, fix_corners, metric_flag, asymmetric_flag, aggregation_flag, pref_class, numeric_flag, write_all_flag, active_waitbar, current_run, total_runs, max_time_estimation, min_counter, max_starting_points)

% this function estimates the parameters of the utility function
% according to the choices made by the subject.
% it does all its calculations for a single subject.
% in order to deal with all the subjects in the database, 
% a loop is in place in HPZ_Estimation_Manager, that passes each 
% time the data of a single subject to this function.

% the best "param" is returned seperately for easier use in bootstrap.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% the entire time that user wants to invest on estimation
% (the user specified minutes, we need to translate to seconds)
% (if it is Inf, multiplying by 60 does not affect it)
allowed_time = max_time_estimation * 60;

% choose proper algorithm for the fminsearchbnd optimization procedure
% the algorithm is for interior points
% not displaying error messages
options = optimset('Algorithm','interior-point', 'Display','off');



% subject_data is a matrix with obs_num rows and 4 columns.
% each row is one choice of the subject.
% The first column is the quantity of good 1 chosen by the subject.
% The second column is the quantity of good 2 chosen by the subject.
% The third column is the price of good 1. 
% The fourth column is the price of good 2. 
subject_data = data(1:obs_num,3:6);

% a matrix of the size: obs_num X obs_num
% in the i,j cell we have the expenditure on the bundle
% that was chosen in the i'th observation given the
% prices of the j'th observation.
expenditure = (subject_data(:,1)*subject_data(:,3)' + subject_data(:,2)*subject_data(:,4)')';

% the vector of length obs_num which indicates the 
% endowment in each observation
endowments = diag(expenditure);



% Choi et al. (2007) correction
if fix_corners == true
    % Choi et al. (2007) correction is applied for corner choices
    observations = HPZ_No_Corners (subject_data, obs_num ,1);
else
    observations = subject_data;
end



%% Initializing the structures that hold the estimation results
 
% (this matrix keeps track of all estimations of parameters)
results = zeros(max_starting_points, 2);

% (keeps track of all function values)
criterion = zeros(max_starting_points, 1);

% a 2-D matrix with a size of: 3 x min_fvals 
% (min_fvals is the max number of convergence points,
% it is basically equal to min_counter, unless min_counter == inf, then we 
% take the max number of initial points)
% the 3 columns represent, respectively: 
%               1) 1st Parameter (beta/alpha)
%               2) 2nd parameter (rho/A)
%               3) Function Value (Fval)
optimal_parameter_matrix = zeros(min(min_counter, max_starting_points), 3);



% in risk preferences, when the functional form is 
% DA with CRRA\CARA we wish to have [0,0] as one of 
% the starting points for the optimization routine.
% zero_rho_initial:
% true - [0,0] is one of the initial points.
% false - [0,0] is not one of the initial points.
zero_rho_initial = false;
if ((pref_class == HPZ_Constants.risk_pref) && (function_flag == HPZ_Constants.CRRA_func || function_flag == HPZ_Constants.CARA_func))
    zero_rho_initial = true;
end



% Handle the corner solutions when CRRA is chosen as the utility function: 
% a flag that restricts rho to be in the range of [0,1] for subjects with
% corner solutions where DA2-CRRA is chosen as the utility function
restricted_rho = false;





% Setting initial points for the fminsearchbnd optimization procedure
% (The first two initial points are the point chosen by CFGK (2007) and the
% rest are points chosen randomly.)
% and setting the first (constant) point to be evaluated ((0,0) or (0.5,0)) 
if pref_class == HPZ_Constants.risk_pref   % risk preferences
    
    % initial points for risk preferences 
    initial_points = HPZ_Initial_Points_Risk (param1_restrictions, param2_restrictions, restricted_rho, zero_rho_initial, max_starting_points);

    % base point to have its criterion calculated always, prior to the loop 
    base_point = [0 , 0];
    
elseif pref_class == HPZ_Constants.OR_pref   % other regarding preferences
    
    % initial points for other regarding (OR) preferences
    initial_points = HPZ_Initial_Points_OR (max_starting_points);
    
    % base point to have its criterion calculated always, prior to the loop 
    base_point = [0.5 , 0];
end



% calculating the function value for the initial point (0,0) or (0.5,0)  
if (action_flag == HPZ_Constants.NLLS_action)   % NLLS
    [criterion_base_point] = HPZ_NLLS_Criterion(base_point, endowments, observations, treatment, function_flag, fix_corners, metric_flag, asymmetric_flag, pref_class, numeric_flag);
elseif (action_flag == HPZ_Constants.MMI_action)   % MMI
    [criterion_base_point] = HPZ_MMI_Criterion(base_point, endowments, observations, treatment, function_flag, aggregation_flag, pref_class, numeric_flag);
elseif (action_flag == HPZ_Constants.BI_action)   % BI
    [criterion_base_point] = HPZ_BI_Criterion(base_point, endowments, observations, treatment, function_flag, pref_class, numeric_flag);
end
% entering the result to all the relevant variables and matrices
optimal_parameter_matrix(1,1:2) = base_point;
optimal_parameter_matrix(1,3) = criterion_base_point;
fval_temp_min = criterion_base_point;



% Now we want to check the criterion value for every pair of (beta,rho)
%   such that beta=p-1 and rho=0, when p is an intermediate of two of the 
%   prices ratios (p>=1) that the DM was exposed to.
%   If any of them has a lower (hence better) criterion - we replace what
%   we got earlier with it. In the end, the first row of the matrice and
%   fval_temp_min will contain the best (beta,rho) pair out of the default 
%   initial bundle (0,0) and all the (p-1,0) bundles.
% we check this only for DA-2 with CRRA, and only if beta can be greater than 0 and rho can be equal to 0. 
if (pref_class == HPZ_Constants.risk_pref) && (function_flag == HPZ_Constants.CRRA_func) && (param1_restrictions(2) > 0) && (param2_restrictions(1) <= 0) && (param2_restrictions(2) >= 0)

    % finding the (beta, 0) parameter combination with the best criterion
    [param_temp, criterion_temp] = HPZ_Check_Rho_Zero_Cases(obs_num, endowments, observations, treatment, function_flag, param1_restrictions, fix_corners, metric_flag, aggregation_flag, asymmetric_flag, pref_class, action_flag, numeric_flag);

    if  (criterion_temp < fval_temp_min)
        % if we got a better estimation, we enter the new result to all 
        % the relevant variables and matrices
        optimal_parameter_matrix(1,1:2) = param_temp(1:2);
        optimal_parameter_matrix(1,3) = criterion_temp;
        fval_temp_min = criterion_temp;
    end
end





%% preparations before the loop

% initialization of 2 arbitrary thresholds on the aggregated criterion value.  
% minimal values of function value will be considered equal if they are 
% less than (one of) these thresholds from one another.
% the first one (fval_threshold) is a threshold in distance,
% while the second one(log_fval_threshold) is a threshold in distance of
% logarithmic (with base 2) values.
% the second one is mainly needed for NLLS with euclidean metric, since the
% criterion's size depends on the number of observations and on the prices.
if action_flag == HPZ_Constants.NLLS_action   % NLLS
    fval_threshold = eps;   % meaningless with NLLS euclidean
    log_fval_threshold = log2(1.001);
elseif (action_flag == HPZ_Constants.MMI_action || action_flag == HPZ_Constants.BI_action)   % MMI / BI
    fval_threshold = 10^(-6);
    log_fval_threshold = log2(1.0001);
end

% number of consecutive equal (+-fval_threshold) 
% function values.
% process will stop once equal_fval_counter 
% will be equal to min_counter.
equal_fval_counter = 1;

% define the waitbar
if (active_waitbar)
    waitbar_name = char(strcat(HPZ_Constants.waitbar_name_estimation, {' '}, '(', HPZ_Constants.current_run_waitbar, {' '}, num2str(current_run), {' '}, HPZ_Constants.total_runs_waitbar, {' '}, num2str(total_runs), ')'));
    waitbar_msg = char(strcat(HPZ_Constants.waitbar_recovery, {' '}, num2str(data(1,1)), {' '}, HPZ_Constants.waitbar_preferences));
    new_bar_val = 0;
    h_wb = wide_waitbar(new_bar_val, waitbar_msg, waitbar_name, HPZ_Constants.waitbar_width_multiplier);
end

% keeps track of accumulative time for all initial points
time_init_accum = 0;
% time for a certain initial point
time_init_pt = 0; %#ok<NASGU>

% main_multiplier is the time designated (out of 1.0) for the estimation itself,
% whereas remainder_multiplier is the time designated for calculating the
% corresponding criterion values after determining the parameters.
main_multiplier = 0.80;
remainder_multplier = 1 - main_multiplier;


%% loop over all starting points 
% (will be halted midway if enough convergence points 
%  were found or if time limitation was reached)
for j=1:max_starting_points
        
    %% updating the waitbar
    % update the value of the progress bar
    % the value is the maximum between the estimated progress according
    % to counting convergence points, the estimated progress according to
    % counting initial points, and the estimated progress according to 
    % time limit
    if (active_waitbar)
        new_bar_val = max([equal_fval_counter/min_counter , j/max_starting_points , time_init_accum/allowed_time]);
        waitbar(new_bar_val, h_wb);
    end

    % start clock ticking, in order to find the total estimation time 
    % for this initial point
    t_init = tic;

    % check if estimation should continue - that is, if we should try more
    % initial points, based on 2 factors:
    %   1) required number of convergence points have not yet been reached, AND
    %   2) the amount of time passed so far plus the estimated amount
    %      of time for the next initial point (1/j from the accumulated)
    %      is less than the required time limitation.
    %   (also, if j==1 we always continue, cause there must be at least 1
    %    estimation attempt...)
    % if TRUE, continue processing.
    % else, do nothing (the loop continues till
    % j==max_starting_points, but does nothing
    if (j == 1) || ( (equal_fval_counter < min_counter) && (time_init_accum*(j+1)/j < allowed_time) )
        
        % the range of the parameters - the minimal and maximal values of each of them:
        min_values = [param1_restrictions(1), param2_restrictions(1)];
        max_values = [param1_restrictions(2), param2_restrictions(2)];
        
        % computations take place here
        if (action_flag == HPZ_Constants.NLLS_action)   % NLLS
            [results(j,:), criterion(j)] = fminsearchbnd(@(param) HPZ_NLLS_Criterion(param, endowments, observations, treatment, function_flag, fix_corners, metric_flag, asymmetric_flag, pref_class, numeric_flag), initial_points(j,:), min_values, max_values, options); 
        elseif (action_flag == HPZ_Constants.MMI_action)   % MMI
            [results(j,:), criterion(j)] = fminsearchbnd(@(param) HPZ_MMI_Criterion(param, endowments, observations, treatment, function_flag, aggregation_flag, pref_class, numeric_flag), initial_points(j,:), min_values, max_values, options); 
        elseif (action_flag == HPZ_Constants.BI_action)   % BI
            [results(j,:),criterion(j)] = patternsearch(@(param) HPZ_BI_Criterion(param, endowments, subject_data, treatment, function_flag, pref_class, numeric_flag), initial_points(j,:), [], [], [], [], min_values, max_values, [], options);
        end
        
        % here we check if the resulting criterion is equal 
        % (+- fval_threshold) to computed minimum criterion 
        if ( fval_temp_min < criterion(j) - fval_threshold ) && ( log2(fval_temp_min) < log2(criterion(j)) - log_fval_threshold ) 
            % do nothing - bigger than current minimum and outside the threshold
        elseif ( fval_temp_min > criterion(j) + fval_threshold ) && ( log2(fval_temp_min) > log2(criterion(j)) + log_fval_threshold )
            % smaller than current minimum and the current minimum is outside its
            % threshold
            equal_fval_counter = 1;
            optimal_parameter_matrix(1 , 1:2) = results(j,:);
            optimal_parameter_matrix(1 , 3  ) = criterion(j);
            fval_temp_min = criterion(j);
            % note: we do not bother to erase the unneeded data from the matrix,
            % because we can keep track using equal_fval_counter
        else
            % we need to add it to the current optimal matrix, and also to sort it,
            % and if it is smaller than the current minimum - possibly get rid of
            % some of the old optimal points that are now outside the new threshold
            equal_fval_counter = equal_fval_counter + 1;
            optimal_parameter_matrix(equal_fval_counter , 1:2) = results(j , 1:2);
            optimal_parameter_matrix(equal_fval_counter , 3  ) = criterion(j);
            % sort the optimal parameter matrix by fval
            optimal_parameter_matrix(1:equal_fval_counter , 1:3) = sortrows(optimal_parameter_matrix(1:equal_fval_counter , 1:3) , 3);

            if ( fval_temp_min <= criterion(j) )
                % we just need to increment the number of optimal points,
                % but we already did that
            else
                % there's a new minimum and a new threshold, therefore we need to
                % count from zero
                fval_temp_min = criterion(j);
                k = 1;
                while (k < equal_fval_counter)
                    if  ( criterion(k+1) > fval_temp_min + fval_threshold ) && ( log2(criterion(k+1)) > log2(fval_temp_min) + log_fval_threshold )
                        break;
                    else
                        k = k + 1;
                    end
                end
                equal_fval_counter = k;
                % note: we do not bother to erase the unneeded data from the matrix,
                % because we can keep track using equal_fval_counter
            end
        end
        
        
        
        % measure the time for each initial point
        time_init_pt = toc(t_init);
        % measure the accumulative time for all initial points so far
        time_init_accum = time_init_accum + time_init_pt;
        
    end
 
end % end of loop



% update the waitbar
if (active_waitbar)
    new_bar_val = main_multiplier;
    waitbar(new_bar_val, h_wb);
end


    
% sorting the matrix of the best points discovered
optimal_parameter_matrix = sortrows(optimal_parameter_matrix(1:equal_fval_counter,:), 3);

% whether the user wants all of the results with approximately the same
% function value printed, or only the result with the best (lowest)
% function value amongst them
if write_all_flag == false
    % show only the best estimation
    max_index = 1;
elseif write_all_flag == true
    % show all resulting estimations
    max_index = equal_fval_counter; 
end

% initializing the output matrix
final_output = zeros(max_index, 8);    

% for each point we store in final_output the values 
% that we already calculated (point coordinates), 
% as well as the NLLS criterions (using Euclidean metric and using CFGK metric), 
% three MMI criterions (rows 5-7, for MAX, MEAN and AVGSSQ aggregators), 
% and the BI criterion.  
for i = 1 : max_index

%     % this code was replaced by a better code about 10 lines ahead
%     final_output(i,1) = optimal_parameter_matrix(i, 1);
%     final_output(i,2) = optimal_parameter_matrix(i, 2);
    
    % NLLS Criterion (Euclidean metric, Choi et al. (2007) metric)
    [final_output(i,3), final_output(i,4), param_NLLS] = HPZ_NLLS_Metrics (optimal_parameter_matrix(i,1:2), endowments, observations, treatment, function_flag, fix_corners, asymmetric_flag, pref_class, numeric_flag);
    
    % MMI Criterion (Max Waste, Mean Waste, Sum of Squares Wastes)
    [final_output(i,5), final_output(i,6), final_output(i,7), param_MMI] = HPZ_MMI_Aggregates(optimal_parameter_matrix(i,1:2), endowments, observations, treatment, function_flag, pref_class, numeric_flag);
    
    % BI Criterion
    [final_output(i,8), param_BI] = HPZ_BI_Criterion(optimal_parameter_matrix(i,1:2), endowments, observations, treatment, function_flag, pref_class, numeric_flag);
    
    % When we perform analytic estimation, we somtimes round number,
    % e.g. in Risk (DA) we round beta to -1 when it is close to 1, and in
    % CRRA we round rho to 0 when it is too close to 0, and etc.
    % Say that the minimization algorithm stopped in beta = -0.99999999,
    % and say that we round beta < -0.999999 to -1, then we want the
    % results files to have -1 printed in them, and not -0.99999999.
    if (action_flag == HPZ_Constants.NLLS_action)
        final_output(i,1:2) = param_NLLS;
    elseif (action_flag == HPZ_Constants.MMI_action)
        final_output(i,1:2) = param_MMI;
    elseif (action_flag == HPZ_Constants.BI_action)
        final_output(i,1:2) = param_BI;
    end
    
    % update the waitbar
    if (active_waitbar)
        waitbar(new_bar_val + ((i/max_index)*remainder_multplier));
    end
end



% we need these values to make bootstrap code easier, for example
param = final_output(1,1:2);

if (action_flag == HPZ_Constants.NLLS_action)
    if (metric_flag == HPZ_Constants.euclidean_metric)
        criterion_index = 3;
    elseif (metric_flag == HPZ_Constants.CFGK_metric)
        criterion_index = 4;
    end
elseif (action_flag == HPZ_Constants.MMI_action)
    if (aggregation_flag == HPZ_Constants.MMI_Max)
        criterion_index = 5;
        elseif (aggregation_flag == HPZ_Constants.MMI_Mean)
        criterion_index = 6;
    elseif (aggregation_flag == HPZ_Constants.MMI_AVGSSQ)
        criterion_index = 7;
    end
elseif (action_flag == HPZ_Constants.BI_action)
    criterion_index = 8;
end
main_criterion = final_output(1,criterion_index);


% close the waitbar
if (active_waitbar)
    close(h_wb);
end

end