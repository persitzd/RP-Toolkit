function HPZ_Estimation_Manager(data_matrix, subjects_index, choice_set_type, action_flag, pref_class, function_flag, numeric_flag, param1_restrictions, param2_restrictions, fix_corners, metric_flag, aggregation_flag, max_time_estimation, min_counter, parallel_flag, output_file_config, write_all_flag, bootstrap_flag, file_val_str, residual_flag, in_sample_flag, out_sample_flag, print_precision, bootstrap_sample_sizes, bootstrap_significance_level, BI_threshold, max_starting_points, one_residuals_file, debugger_mode, waitbar_settings, main_folder_for_results, current_run, total_runs)

% this function is the container for the whole parameter estimation module - 
% all parameter estimation procedures (estimation, printing results to
% files, bootstrap, residuals) are handled from here.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 

% explanations about these can be found in: HPZ_Global_Variables 
global warnings_nan
global warnings_minus_inf
global warnings_plus_inf
global warnings_bigger_than_1
global warnings_smaller_than_0
global current_subject



% precision of numbers when printing
precision_string = strcat('%10.', num2str(print_precision), 'g');



% 'stable' is essential to keep the original order of subjects as it is
% in the file data
subjects = unique(data_matrix(:,1), 'stable');  % array of subjects' IDs
num_subjects = length(subjects);
% chosen subjects
chosen_subjects = subjects(subjects_index);
chosen_subjects_num = length(subjects_index);
% an array of location indices - the i'th element 		
% is the index of the first observation of the subject i 
[rows,~] = size(data_matrix);
first_obs_row = zeros(num_subjects, 1);
counter_subjects = 1;
for i=1:rows
    if data_matrix(i,2) == 1
        first_obs_row (counter_subjects) = i;
        counter_subjects = counter_subjects + 1;
    end
end



% determining whether to use a single waitbar or a separate waitbar per subject 
if chosen_subjects_num < waitbar_settings(1) 
    single_waitbar = false;
    active_waitbar = true;
    % per subject waitbar for residuals and for bootstrap
    residuals_waitbar = true;
    bootstrap_waitbar = true;
else
    single_waitbar = true;
    active_waitbar = false;
    % per subject waitbar for residuals and for bootstrap
    residuals_waitbar = waitbar_settings(3);
    bootstrap_waitbar = waitbar_settings(4);
end




% initialize (or re-initialize) global variables for this estimation
HPZ_Global_Variables(chosen_subjects_num);



if (bootstrap_flag)
    % analytic approach is computed faster than numeric, 
    % therefore we allow a bigger sample.
    if (numeric_flag == HPZ_Constants.analytic)
        number_of_samples = bootstrap_sample_sizes(1);
    elseif (numeric_flag == HPZ_Constants.numeric)
        number_of_samples = bootstrap_sample_sizes(2);
    elseif (numeric_flag == HPZ_Constants.semi_numeric)
        number_of_samples = bootstrap_sample_sizes(3);
    end
end





%% prepare results file/s, including column headers
[file_handle, file_name_str, file_residuals_path_str, file_warnings_path_str] = HPZ_Write_Result_File_Initializer(output_file_config, action_flag, pref_class, function_flag, file_val_str, bootstrap_flag, bootstrap_significance_level, main_folder_for_results); %#ok<ASGLU>

% make sure the file closes when execution stops for any reason 
cleanup_close_file = onCleanup(@() fclose(file_handle));

% getting the list of column headers for the residual file/files
if (residual_flag)
    residuals_col_headers = HPZ_Estimation_Residuals_Headers (in_sample_flag, out_sample_flag, pref_class, action_flag, function_flag, metric_flag, aggregation_flag);
end


%% Set matlabpool (Parallel Computing)
if parallel_flag
    % close all other session          
    poolobj = gcp('nocreate');
    delete(poolobj);

    % use all resources
    parpool;
end





% code for a single waitbar instead of multiple waitbars (part 1 out of 3) 
if (single_waitbar)
    % disable the per-subject waitbar
    active_waitbar = false;
    % define the waitbar
    waitbar_name = char(strcat(HPZ_Constants.waitbar_name_estimation, {' - '}, HPZ_Constants.all_actions_short_names{action_flag}, {' , '}, HPZ_Constants.all_numeric_options_names{numeric_flag+1}, {' '}, '(', HPZ_Constants.current_run_waitbar, {' '}, num2str(current_run), {' '}, HPZ_Constants.total_runs_waitbar, {' '}, num2str(total_runs), ')'));
    waitbar_msg = char(strcat(HPZ_Constants.waitbar_recovery, 's', {' '}, HPZ_Constants.waitbar_preferences));
    new_bar_val = 0;
    h_wb = wide_waitbar(new_bar_val, {waitbar_msg, '', ''}, waitbar_name, HPZ_Constants.waitbar_width_multiplier);
end



%% This loop goes through the subjects chosen one by one and 		
% estimates the parameters or the indices for that subject, according 		
% to the method chosen by the user, then prints the results for that
% subject to the results file/s
for i=1:chosen_subjects_num
    
    
    % code for a single waitbar instead of multiple waitbars (part 2 out of 3) 
    if (single_waitbar)
        % update the waitbar
        new_bar_val = (i-1) / chosen_subjects_num;
        current_subject_str = char(strcat({'Estimating Subject '}, num2str(chosen_subjects(i))));
        waitbar(new_bar_val, h_wb, {waitbar_msg, char(strcat({'Completed '}, num2str(i-1), {' Subjects out of '}, num2str(chosen_subjects_num))), current_subject_str});
    end
    
    
    % this index is for the warnings module; we keep a separate count of
    % warnings for each subject
    current_subject = i;
    
    
    % extracting the data that is necessary for all methods and calculations
    if subjects_index(i) < num_subjects
        data = data_matrix(first_obs_row(subjects_index(i)):((first_obs_row(subjects_index(i)+1))-1),:);
        obs_num = data_matrix((first_obs_row(subjects_index(i)+1))-1,2);
    else
        data = data_matrix(first_obs_row(subjects_index(i)):rows,:);
        obs_num = data_matrix(rows,2);
    end
    
    asymmetric_flag = 1;    % As is
    treatment = 1;          % As is
    
    % Perform the estimation (NLLS/MMI/BI)
    [param, main_criterion, final_output] = HPZ_Estimation (data, obs_num, choice_set_type, action_flag, treatment, function_flag, param1_restrictions, param2_restrictions, ...
                                                fix_corners, metric_flag, asymmetric_flag, aggregation_flag, pref_class, ...
                                                numeric_flag, write_all_flag, max_time_estimation, min_counter, max_starting_points, BI_threshold, debugger_mode, active_waitbar, current_run, total_runs);

    % number of optimal results to be shown
    [rows_out,~] = size(final_output);

    if bootstrap_flag == true
        % if the user desired and set so, print the 95% upper	        
        % and lower estimates of the samples.
        % (not necessarily 95%, but in general it is (1-significance_level)) 
        % the samples take randomly obs_num observations from the
        % subject's observations, with repeat, and find the estimation 
        % for this new set of observations.
        [me_1, me_2, se_1, se_2, le_1, le_2, ue_1, ue_2] = HPZ_Bootstrap_SE_Helper(data, choice_set_type, treatment, number_of_samples, action_flag, ...
                                                        function_flag, param1_restrictions, param2_restrictions, fix_corners, metric_flag, ...
                                                        aggregation_flag, asymmetric_flag, pref_class, numeric_flag, bootstrap_significance_level, rows_out, ...
                                                        max_time_estimation, min_counter, max_starting_points, BI_threshold, debugger_mode, bootstrap_waitbar, current_run, total_runs);

        final_output_new = [final_output, me_1, me_2, se_1, se_2, le_1, le_2, ue_1, ue_2];

    else

        final_output_new = final_output;

    end

    for j=1:rows_out

        % writing the results to the results file
        HPZ_Write_Result_File_Finalizer(file_handle, output_file_config, final_output_new, ...
                                        j, num2str(data(1,1)), obs_num, bootstrap_flag, print_precision);

    end

    if (residual_flag)
        % if the user desired and set so, print the residuals to a separate file
        Mat = HPZ_Estimation_Residuals (action_flag, data, obs_num, choice_set_type, treatment, function_flag, param1_restrictions, param2_restrictions, fix_corners, metric_flag, ...
                                    aggregation_flag, asymmetric_flag, in_sample_flag, out_sample_flag, param(1), param(2), main_criterion, pref_class, numeric_flag, ...
                                    max_time_estimation, min_counter, max_starting_points, BI_threshold, debugger_mode, residuals_waitbar, current_run, total_runs);

        % make sure the Results directory exists, if not - create it
        dir_exists = exist(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir) , 'dir');
        if ~dir_exists
            mkdir(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir));
        end
        
        if (one_residuals_file)   % print to one file, different sheets
            
            % printing column headers to the file
            xlrange = char(strcat('A1:', HPZ_Constants.excel_alphabet(max(size(residuals_col_headers))), '1'));
            xlswrite(file_residuals_path_str, residuals_col_headers, i, xlrange);
            xlswrite(file_residuals_path_str, residuals_col_headers, i, xlrange);
            
            % printing the results to the file
            [mat_rows, mat_cols] = size(Mat);
            xlrange = char(strcat('A2:', HPZ_Constants.excel_alphabet(mat_cols), num2str(mat_rows+1)));
            xlswrite(file_residuals_path_str, Mat, i, xlrange);
            
        else   % print to each subject to a separate file
            
            % name of residuals file of this subject - we make sure that there isn't an existing file by the same name    
            file_residuals_str_i = char(strcat(file_residuals_path_str, {' - Subject '}, num2str(data(1,1))));
            while (true)
                if ~exist(strcat(file_residuals_str_i, '.csv'), 'file')
                    break
                else
                    file_residuals_str_i = char(strcat(file_residuals_str_i, {' - '}, generate_random_str));
                end
            end
            % printing the results to the file
            print_matrix_to_file (file_residuals_str_i, Mat, precision_string, col_headers, 0);
        end                  

    end   % end of residuals
    
    
%     % if there were any unreasonable criterion values - print a summarizing warning 
%     if sum([warnings_nan(i), warnings_minus_inf(i) ,warnings_plus_inf(i), warnings_bigger_than_1(i) , warnings_smaller_than_0(i)])
%         if (action_flag == HPZ_Constants.NLLS_action)
%             warning('In Subject %s , there were (accumulated through all the iterations. criterions per observation) %d NaN criterions, %d -Inf criterions, %d +Inf criterions.', num2str(subjects(i)), warnings_nan(i), warnings_minus_inf(i) ,warnings_plus_inf(i));
%         elseif (action_flag == HPZ_Constants.MMI_action) || (action_flag == HPZ_Constants.BI_action)
%             warning('In Subject %s , there were (accumulated through all the iterations. criterions per observation) %d NaN criterions, %d -Inf criterions, %d +Inf criterions, %d criterions > 1, %d criterions < 0.', num2str(subjects(i)), warnings_nan(i), warnings_minus_inf(i) ,warnings_plus_inf(i), warnings_bigger_than_1(i) , warnings_smaller_than_0(i));
%         end
%     end

    
end   % end of loop over chosen subjects



% code for a single waitbar instead of multiple waitbars (part 3 out of 3) 
if (single_waitbar)
    % close the waitbar
    close(h_wb);
end



% if there were any unreasonable criterion values 
if sum(sum([warnings_nan, warnings_minus_inf ,warnings_plus_inf, warnings_bigger_than_1 , warnings_smaller_than_0]))
    % print the warnings for this estimation to a file
    if (action_flag == HPZ_Constants.NLLS_action)
        warning_matrix = [subjects(subjects_index) , warnings_nan , warnings_minus_inf , warnings_plus_inf];
        warning_headers = {'Subject' , 'NaN criterions' , 'Minus -Inf criterions' , 'Plus +Inf criterions'};
    elseif (action_flag == HPZ_Constants.MMI_action) || (action_flag == HPZ_Constants.BI_action)
        warning_matrix = [subjects(subjects_index) , warnings_nan , warnings_minus_inf , warnings_plus_inf , warnings_bigger_than_1 , warnings_smaller_than_0];
        warning_headers = {'Subject' , 'NaN criterions' , 'Minus -Inf criterions' , 'Plus +Inf criterions', 'criterions > 1' , 'criterions < 0'};
    end
    print_matrix_to_file (file_warnings_path_str, warning_matrix, precision_string, warning_headers, 0);
end



end