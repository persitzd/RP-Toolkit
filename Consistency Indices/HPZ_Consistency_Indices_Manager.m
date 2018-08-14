function HPZ_Consistency_Indices_Manager(data_matrix, subjects_index, GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags, one_residuals_file, max_time_estimation, print_precision, waitbar_options, main_folder_for_results, current_run, total_runs) %#ok<INUSL>

% this function is the container for the whole inconsistency indices module - 
% all inconsistency indices procedures (consistency check, printing results to
% files, calculating inconsistency indices, residuals) are handled from here.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% determining whether to use a single waitbar or a seperate waitbar per subject 
if (waitbar_options(1) == HPZ_Constants.waitbar_none)
    single_waitbar = false;
    active_waitbar = false;
elseif (waitbar_options(1) == HPZ_Constants.waitbar_per_subject)
    single_waitbar = false;
    active_waitbar = true;
elseif (waitbar_options(1) == HPZ_Constants.waitbar_single)
    single_waitbar = true;
    active_waitbar = false;
end



% these are just for convenience
consistency_basic_flags = [GARP_flags(1) , AFRIAT_flags(1) , VARIAN_flags(1) , HOUTMAN_flags(1) , MPI_flags(1)];
consistency_residuals_flags_temp = [GARP_flags(2) , AFRIAT_flags(2) , VARIAN_flags(2) , HOUTMAN_flags(2) , MPI_flags(2)];
consistency_residuals_flags = min(consistency_basic_flags, consistency_residuals_flags_temp);
consistency_in_sample_flags = [GARP_flags(3) , AFRIAT_flags(3) , VARIAN_flags(3) , HOUTMAN_flags(3) , MPI_flags(3)]; %#ok<NASGU>
consistency_out_sample_flags = [GARP_flags(4) , AFRIAT_flags(4) , VARIAN_flags(4) , HOUTMAN_flags(4) , MPI_flags(4)]; %#ok<NASGU>



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





%% Setting the main result file

% make sure the Results directory exists, if not - create it
dir_exists = exist(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir) , 'dir');
if ~dir_exists
    mkdir(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir));
end

% naming and creating the main result file - we make sure that there isn't an existing file by the same name    
while (true)
    file_name_str = char(strcat(HPZ_Constants.Consistency_action_file_name, {' - '}, 'Results - Date', {' '}, date, {' - '}, generate_random_str));
    file_path_str = strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir, '/', file_name_str, '.csv');
    if ~exist(file_path_str, 'file')
        break
    end
end
file_handle = fopen(file_path_str, 'w');

% make sure the file closes when execution stops for any reason 
cleanup_close_file = onCleanup(@() fclose(file_handle));

% name and headers for residuals file
if sum(consistency_residuals_flags) ~= 0
    % getting the list of column headers for the residual file/files
    residuals_col_headers = HPZ_Consistency_Indices_Residuals_Headers (GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags);
    % name of residuals file - we make sure that there isn't an existing file by the same name    
    while (true)
        consistency_residuals_str = char(strcat(HPZ_Constants.Consistency_action_file_name, {' - '}, 'Results Residuals - Date', {' '}, date, {' - '}, generate_random_str));
        consistency_residuals_path_str = strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir, '/', consistency_residuals_str);
        if ~exist(strcat(consistency_residuals_path_str, '.csv'), 'file')
            break
        end
    end
end



%% finalize the output main file headers for distinct results
fprintf(file_handle, '%s,', 'Subject');
if (GARP_flags(1) == 1)
    fprintf(file_handle, '%s,%s,%s,%s,%s,%s,', ...
        'WARP Violations', '#Pairs', ...
        'GARP Violations', '#Pairs', ...
        'SARP Violations', '#Pairs');
end
if (AFRIAT_flags(1) == 1)
    fprintf(file_handle, '%s,', 'AFRIAT_Index');
end
if (VARIAN_flags(1) == 1)
    fprintf(file_handle, '%s,%s,%s,%s,', ...
        'VARIAN Index Min', 'VARIAN Index Mean', 'VARIAN Index AVG(SSQ)', 'Is exact');
end
if (HOUTMAN_flags(1) == 1)
    fprintf(file_handle, '%s,%s,', ...
        'Houtman Maks Index', 'Is exact');
end
if (MPI_flags(1) == 1)
    fprintf(file_handle, '%s,%s,', ...
        'MPI Index Mean', 'MPI Index Median');
end
fprintf(file_handle, '\n');





% code for a single waitbar instead of multiple waitbars (part 1 out of 3) 
if (single_waitbar)
    % disable the per-subject waitbar
    active_waitbar = false;
    % define the waitbar
    basic_flags = [GARP_flags(1) , AFRIAT_flags(1) , VARIAN_flags(1) , HOUTMAN_flags(1) , MPI_flags(1)];
    indices_strings = {' GARP' , ' Afriat', ' Varian', ' Houtman-Maks', ' Money-Pump'};
    for i=1:length(basic_flags)
        if basic_flags(i) == 0
            indices_strings{i} = '';
        end
    end
    waitbar_name = char(strcat(HPZ_Constants.waitbar_name_calculation, {' -'}, indices_strings{:}, {' '}, '(', HPZ_Constants.current_run_waitbar, {' '}, num2str(current_run), {' '}, HPZ_Constants.total_runs_waitbar, {' '}, num2str(total_runs), ')'));
    waitbar_msg = char(strcat(HPZ_Constants.waitbar_calculation, 's', {' '}, HPZ_Constants.waitbar_indices));
    new_bar_val = 0;
    h_wb = wide_waitbar(new_bar_val, {waitbar_msg, '', ''}, waitbar_name, HPZ_Constants.waitbar_width_multiplier);
end



%% This loop goes through the subjects chosen one by one and 		
% calculates the indices for that subject, then prints the results for that
% subject to the results file/s
for i=1:chosen_subjects_num
    
    
    % code for a single waitbar instead of multiple waitbars (part 2 out of 3) 
    if (single_waitbar)
        % update the waitbar
        new_bar_val = (i-1) / chosen_subjects_num;
        current_subject_str = char(strcat({'Calculating Subject '}, num2str(chosen_subjects(i))));
        waitbar(new_bar_val, h_wb, {waitbar_msg, char(strcat({'Completed '}, num2str(i-1), {' Subjects out of '}, num2str(chosen_subjects_num))), current_subject_str});
    end
    
    
    % extracting the data that is necessary for all methods and calculations
    if subjects_index(i) < num_subjects
        data = data_matrix(first_obs_row(subjects_index(i)):((first_obs_row(subjects_index(i)+1))-1),:);
        obs_num = data_matrix((first_obs_row(subjects_index(i)+1))-1,2); %#ok<NASGU>
    else
        data = data_matrix(first_obs_row(subjects_index(i)):rows,:);
        obs_num = data_matrix(rows,2); %#ok<NASGU>
    end
    
    
    
    % Calculating consistency tests and inconsistency indices

    % the function can return the disaggregated matrices of violations, we
    % don't report it.
    [~,~,~,~,~,~,~,~,~, VIO_PAIRS, VIOLATIONS, AFRIAT, VARIAN, Var_exact, HM, HM_exact, MPI, Mat] = ...
                            HPZ_Subject_Consistency (data, GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags, active_waitbar, current_run, total_runs);

    % precision of numbers when printing
    precision_string = strcat('%10.', num2str(print_precision), 'g');

    %% printing the results
    fprintf(file_handle, '%s,', num2str(data(1,1)));
    if (GARP_flags(1) == 1)
        fprintf(file_handle, '%s,%s,%s,%s,%s,%s,', ...
            num2str(VIOLATIONS(1), precision_string), ...
            num2str(VIO_PAIRS(1), precision_string), ...
            num2str(VIOLATIONS(2), precision_string), ...
            num2str(VIO_PAIRS(2), precision_string), ....
            num2str(VIOLATIONS(3), precision_string), ...
            num2str(VIO_PAIRS(3), precision_string));
    end
    if (AFRIAT_flags(1) == 1)
        fprintf(file_handle, '%s,', num2str(AFRIAT, precision_string));
    end
    if (VARIAN_flags(1) == 1)
        fprintf(file_handle, '%s,%s,%s,%s,', ...
            num2str(VARIAN(1), precision_string), ...
            num2str(VARIAN(2), precision_string), ...
            num2str(VARIAN(3), precision_string), ...
            num2str(Var_exact, precision_string));
    end
    if (HOUTMAN_flags(1) == 1)
        fprintf(file_handle, '%s,%s,', ...
            num2str(HM, precision_string), ...
            num2str(HM_exact, precision_string));
    end
    if (MPI_flags(1) == 1)
        fprintf(file_handle, '%s,%s', ...
            num2str(MPI(1), precision_string), ...
            num2str(MPI(2), precision_string));
    end
    fprintf(file_handle, '\n');

    %% printing residuals
    if sum(consistency_residuals_flags) ~= 0
        % if the user desired and set so, print the residuals to a separate file

        % make sure the Results directory exists, if not - create it
        dir_exists = exist(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir) , 'dir');
        if ~dir_exists
            mkdir(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir));
        end
        
        if (one_residuals_file)   % print to one file, different sheets
            
            % printing column headers to the file
            xlrange = char(strcat('A1:', HPZ_Constants.excel_alphabet(max(size(residuals_col_headers))), '1'));
            xlswrite(consistency_residuals_path_str, residuals_col_headers, i, xlrange);

            % printing the results to the file
            [mat_rows, mat_cols] = size(Mat);
            xlrange = char(strcat('A2:', HPZ_Constants.excel_alphabet(mat_cols), num2str(mat_rows+1)));
            xlswrite(consistency_residuals_path_str, Mat, i, xlrange);
            
        else   % print to each subject to a separate file
            
            % name of residuals file of this subject - we make sure that there isn't an existing file by the same name    
            consistency_residuals_path_str_i = char(strcat(consistency_residuals_path_str, {' - Subject '}, num2str(data(1,1))));
            while (true)
                if ~exist(strcat(consistency_residuals_path_str_i, '.csv'), 'file')
                    break
                else
                    consistency_residuals_path_str_i = char(strcat(consistency_residuals_path_str_i, {' - '}, generate_random_str));
                end
            end
            % printing the results to the file
            print_matrix_to_file (consistency_residuals_path_str_i, Mat, precision_string, residuals_col_headers, 0);
        
        end
        
    end   % end of residuals printing
    

end   % end of loop over chosen subjects



% code for a single waitbar instead of multiple waitbars (part 3 out of 3) 
if (single_waitbar)
    % close the waitbar
    close(h_wb);
end



end