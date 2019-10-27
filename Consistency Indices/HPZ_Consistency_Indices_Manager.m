function HPZ_Consistency_Indices_Manager(data_matrix, subjects_index, Graph_flags, power_test_settings, GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags, max_time_estimation, print_precision, one_residuals_file, waitbar_settings, Varian_algorithm_settings, main_folder_for_results, current_run, total_runs) %#ok<INUSL>

% this function is the container for the whole inconsistency indices module - 
% all inconsistency indices procedures (consistency check, printing results to
% files, calculating inconsistency indices, residuals) are handled from here.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% these are just for convenience
consistency_basic_flags = [GARP_flags(1) , AFRIAT_flags(1) , VARIAN_flags(1) , HOUTMAN_flags(1) , MPI_flags(1)];
consistency_residuals_flags_temp = [GARP_flags(2) , AFRIAT_flags(2) , VARIAN_flags(2) , HOUTMAN_flags(2) , MPI_flags(2)];
consistency_residuals_flags = min(consistency_basic_flags, consistency_residuals_flags_temp);
consistency_in_sample_flags = [GARP_flags(3) , AFRIAT_flags(3) , VARIAN_flags(3) , HOUTMAN_flags(3) , MPI_flags(3)]; %#ok<NASGU>
consistency_out_sample_flags = [GARP_flags(4) , AFRIAT_flags(4) , VARIAN_flags(4) , HOUTMAN_flags(4) , MPI_flags(4)]; %#ok<NASGU>
% initialization to avoid errors
Mat_num_of_columns = 0;

% whether to perform power test
perform_power_test_flag = power_test_settings(1);



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
if chosen_subjects_num < waitbar_settings(1) && ~perform_power_test_flag
    single_waitbar = false;
    active_waitbar = true;
    % per subject waitbar for residuals
    residuals_waitbar = true;
else
    single_waitbar = true;
    active_waitbar = false;
    % per subject waitbar for residuals and for bootstrap
    residuals_waitbar = waitbar_settings(2);
end
% % DUE TO A TEMPORARY BUG I HAD IN MY COMPUTER THAT DID NOT ALLOW WAITBARS...  
% % IT WAS EVENTUALLY FIXED BY REMOVING SOME WINDOWS UPDATES >:(
% single_waitbar = false;
% active_waitbar = false;
% residuals_waitbar = false;




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
    [residuals_col_headers , Mat_num_of_columns] = HPZ_Consistency_Indices_Residuals_Headers (GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags);
    % name of residuals file - we make sure that there isn't an existing file by the same name    
    while (true)
        consistency_residuals_str = char(strcat(HPZ_Constants.Consistency_action_file_name, {' - '}, 'Results Residuals - Date', {' '}, date, {' - '}, generate_random_str));
        consistency_residuals_path_str = strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir, '/', consistency_residuals_str);
        if ~exist(strcat(consistency_residuals_path_str, '.csv'), 'file')
            break
        end
    end
end



% we will add this number in the end of each .png / .fig file,
% in order not to run-over possible existing images
Graph_flags(end+1) = str2num(generate_random_str); %#ok<ST2NM>



%% finalize the output main file headers for distinct results
fprintf(file_handle, '%s,', 'Subject');
fprintf(file_handle, '%s,', 'Num of Observations');
if (GARP_flags(1) == 1)
    fprintf(file_handle, '%s,%s,%s,%s,%s,%s,', ...
        'WARP Violations', ...
        'GARP Violations', 'GARP Pairs', ...
        'SARP Violations');
end
if (AFRIAT_flags(1) == 1)
    fprintf(file_handle, '%s,', 'AFRIAT Index');
end
if (VARIAN_flags(1) == 1)
    if (VARIAN_flags(5) == 1)
        fprintf(file_handle, '%s,%s,%s,%s,', ...
            'VARIAN Index Mean', 'VARIAN Index Mean Lower Bound', 'VARIAN Index Mean Approximate', 'VARIAN Index Mean Upper Bound');
    end
    if (VARIAN_flags(6) == 1)
        fprintf(file_handle, '%s,%s,%s,%s,', ...
            'VARIAN Index AVG(SSQ)', 'VARIAN Index AVG(SSQ) Lower Bound', 'VARIAN Index AVG(SSQ) Approximate', 'VARIAN Index AVG(SSQ) Upper Bound');
    end
    if (VARIAN_flags(7) == 1)
        fprintf(file_handle, '%s,%s,%s,%s,', ...
            'VARIAN Index Max', 'VARIAN Index Max Lower Bound', 'VARIAN Index Max Approximate', 'VARIAN Index Max Upper Bound');
    end
end
if (HOUTMAN_flags(1) == 1)
    fprintf(file_handle, '%s,', ...
        'Houtman Maks Index');
end
if (MPI_flags(1) == 1)
    fprintf(file_handle, '%s,%s,', ...
        'MPI Index Mean', 'MPI Index Median');
end
% Power Test
if perform_power_test_flag
    
    fprintf(file_handle, '%s,%s,', ...
        'Number of Simulations', 'Number of Simulations Satisfying GARP');
    if (GARP_flags(1) == 1)
        fprintf(file_handle, '%s,%s,%s, %s,%s,%s, %s,%s,%s, %s,%s,%s,', ...
            'WARP Violations - simulation median', 'WARP Violations - simulation mean', 'WARP Violations - simulation std', ...
            'GARP Violations - simulation median', 'GARP Violations - simulation mean', 'GARP Violations - simulation std', ...
            'GARP Pairs - simulation median', 'GARP Pairs - simulation mean', 'GARP Pairs - simulation std', ...
            'SARP Violations - simulation median', 'SARP Violations - simulation mean', 'SARP Violations - simulation std');
    end
    if (AFRIAT_flags(1) == 1)
        fprintf(file_handle, '%s,%s,%s,', ...
            'AFRIAT Index - simulation median', 'AFRIAT Index - simulation mean', 'AFRIAT Index - simulation std');
    end
    if (VARIAN_flags(1) == 1)
        if (VARIAN_flags(5) == 1) % Mean Aggregator
            fprintf(file_handle, '%s,%s,%s,', ...
                'VARIAN Index Mean - simulation median', ...
                'VARIAN Index Mean - simulation mean', ...
                'VARIAN Index Mean - simulation std');
        end
        if (VARIAN_flags(6) == 1) % AVG(SSQ) Aggregator
            fprintf(file_handle, '%s,%s,%s,', ...
                'VARIAN Index AVG(SSQ) - simulation median', ...
                'VARIAN Index AVG(SSQ) - simulation mean', ...
                'VARIAN Index AVG(SSQ) - simulation std');
        end
        if (VARIAN_flags(7) == 1) % Maximum Aggregator
            fprintf(file_handle, '%s,%s,%s,', ...
                'VARIAN Index Max - simulation median', ...
                'VARIAN Index Max - simulation mean', ...
                'VARIAN Index Max - simulation std');
        end
    end
    if (HOUTMAN_flags(1) == 1)
        fprintf(file_handle, '%s,%s,%s,', ...
            'Houtman Maks Index median', 'Houtman Maks Index mean', 'Houtman Maks Index std');
    end
%     if (MPI_flags(1) == 1)
%         fprintf(file_handle, '%s,%s,%s, %s,%s,%s,', ...
%             num2str(MPI_stats(1, 1), precision_string), ...
%             num2str(MPI_stats(1, 2), precision_string), ...
%             num2str(MPI_stats(1, 3), precision_string), ...
%             num2str(MPI_stats(2, 1), precision_string), ...
%             num2str(MPI_stats(2, 2), precision_string), ...
%             num2str(MPI_stats(2, 3), precision_string));
%     end
end
% finally, move to the second row
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
    %i
    %subject = chosen_subjects(i)
    
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
        obs_num = data_matrix((first_obs_row(subjects_index(i)+1))-1,2);
    else
        data = data_matrix(first_obs_row(subjects_index(i)):rows,:);
        obs_num = data_matrix(rows,2);
    end
    
%     start_time = now;
%     subject = num2str(data(1,1))
    
    
    % Calculating consistency tests and inconsistency indices

    % the function can return the disaggregated matrices of violations, we
    % don't report it.
    [~,~,~,~,~,~,~,~,~, VIO_PAIRS, VIOLATIONS, AFRIAT, VARIAN_Bounds, HoutmanMaks, MPI, Mat] = ...
                            HPZ_Subject_Consistency (data, Graph_flags, GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags, Mat_num_of_columns, Varian_algorithm_settings, main_folder_for_results, active_waitbar, residuals_waitbar, current_run, total_runs);

%     end_time = now;
%     running_time = datevec(end_time - start_time);
%     months = running_time(2);
%     days = running_time(3);
%     hours = running_time(4);
%     minutes = running_time(5);
%     seconds = running_time(6);
%     fprintf('\ntotal running time was:\n%d months, %d days, %d hours, %d minutes and %.3f seconds.\n',...
%                                                                     months, days, hours, minutes, seconds);                    

    % precision of numbers when printing
    precision_string = strcat('%10.', num2str(print_precision), 'g');

    %% printing the results
    fprintf(file_handle, '%s,', num2str(data(1,1)));
    fprintf(file_handle, '%s,', num2str(obs_num));
    if (GARP_flags(1) == 1)
        fprintf(file_handle, '%s,%s,%s,%s,', ...
            num2str(VIOLATIONS(1), precision_string), ...
            num2str(VIOLATIONS(2), precision_string), ...
            num2str(VIO_PAIRS(2), precision_string), ....
            num2str(VIOLATIONS(3), precision_string));
    end
    if (AFRIAT_flags(1) == 1)
        fprintf(file_handle, '%s,', num2str(AFRIAT, precision_string));
    end
    if (VARIAN_flags(1) == 1)
        if (VARIAN_flags(5) == 1) % Mean Aggregator
            if VARIAN_Bounds(1) == VARIAN_Bounds(3)
                % it is exact
                fprintf(file_handle, '%s,%s,%s,%s,', ...
                    num2str(VARIAN_Bounds(2), precision_string), '-', '-', '-');
            else
                % it is not exact
                fprintf(file_handle, '%s,%s,%s,%s,', ...
                    '-', num2str(VARIAN_Bounds(1), precision_string), ...
                         num2str(VARIAN_Bounds(2), precision_string), ...
                         num2str(VARIAN_Bounds(3), precision_string));
            end
        end
        if (VARIAN_flags(6) == 1) % AVG(SSQ) Aggregator
            if VARIAN_Bounds(4) == VARIAN_Bounds(6)
                % it is exact
                fprintf(file_handle, '%s,%s,%s,%s,', ...
                    num2str(VARIAN_Bounds(5), precision_string), '-', '-', '-');
            else
                % it is not exact
                fprintf(file_handle, '%s,%s,%s,%s,', ...
                    '-', num2str(VARIAN_Bounds(4), precision_string), ...
                         num2str(VARIAN_Bounds(5), precision_string), ...
                         num2str(VARIAN_Bounds(6), precision_string));
            end
        end
        if (VARIAN_flags(7) == 1) % Maximum Aggregator
            if VARIAN_Bounds(7) == VARIAN_Bounds(9)
                % it is exact
                fprintf(file_handle, '%s,%s,%s,%s,', ...
                    num2str(VARIAN_Bounds(8), precision_string), '-', '-', '-');
            else
                % it is not exact
                fprintf(file_handle, '%s,%s,%s,%s,', ...
                    '-', num2str(VARIAN_Bounds(7), precision_string), ...
                         num2str(VARIAN_Bounds(8), precision_string), ...
                         num2str(VARIAN_Bounds(9), precision_string));
            end
        end
    end
    if (HOUTMAN_flags(1) == 1)
        fprintf(file_handle, '%s,', ...
            num2str(HoutmanMaks, precision_string));
    end
    if (MPI_flags(1) == 1)
        fprintf(file_handle, '%s,%s', ...
            num2str(MPI(1), precision_string), ...
            num2str(MPI(2), precision_string));
    end
    %fprintf(file_handle, '\n');

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
    
    
    
    
    
    %% Power Test
    
    % if the user desired to perform Power Test per subject
    if perform_power_test_flag
        
        % call to subject consistency with same calculate-or-not flags,
        % but without creating images and without any residuals
        Graph_flags_for_power_test = Graph_flags; 
        GARP_flags_for_power_test = GARP_flags; 
        AFRIAT_flags_for_power_test = AFRIAT_flags;  
        VARIAN_flags_for_power_test = VARIAN_flags; 
        HOUTMAN_flags_for_power_test = HOUTMAN_flags;  
        MPI_flags_for_power_test = MPI_flags; 
        Graph_flags_for_power_test(:) = 0;   % no Graphs required
        GARP_flags_for_power_test(2:end) = 0; 
        AFRIAT_flags_for_power_test(2:end) = 0;  
        VARIAN_flags_for_power_test(2:4) = 0;   % 5,6,7 are for the three aggregators, which of them to calculate 
        HOUTMAN_flags_for_power_test(2:end) = 0;  
        MPI_flags_for_power_test(2:end) = 0; 
        
        % number of observations and calculation of number of goods -
        % we currently only allow power test for the case of 2 goods
        [obs_num, num_of_columns] = size(data);
        num_of_goods = (num_of_columns - 2) / 2;
        if num_of_goods == 2
        
            base_data_for_power_test = data;
            base_data_for_power_test(:, 3:4) = nan;   % remove the observed choices
            max_quantities = 1 ./ data(:, 5:6);
            safe_bundles = 1 ./ sum(data(:, 5:6), 2);   % safe bundles (only in risk)
            
            num_of_simulations = power_test_settings(2);
            bundle_randomization_type = power_test_settings(3);   % 0 = standard , 1 = given FOSD
            
            % initilizations of list of indexes for each simulation
            VIO_PAIRS_list = nan(3, num_of_simulations);
            VIOLATIONS_list = nan(3, num_of_simulations);
            AFRIAT_list = nan(1, num_of_simulations);
            VARIAN_Bounds_list = nan(9, num_of_simulations);
            HoutmanMaks_list = nan(1, num_of_simulations);
            MPI_list = nan(2, num_of_simulations);
            
            for simulation = 1:num_of_simulations
                
                if bundle_randomization_type == 0
                    % create a data with same budget lines with uniformly random choices   
                    data_for_power_test = base_data_for_power_test;
                    random_uniform_numbers = rand(obs_num, 1);
                    random_uniform_choice_bundle_1 = max_quantities(:,1) .* random_uniform_numbers;
                    random_uniform_choice_bundle_2 = max_quantities(:,2) .* (1-random_uniform_numbers);
                    data_for_power_test(:, 3:4) = [random_uniform_choice_bundle_1 , random_uniform_choice_bundle_2];
                elseif bundle_randomization_type == 1
                    % create a data with same budget lines, with choices randomly 
                    % chosen from points on the budget line that adhere to FOSD-satisfying preferences,
                    % i.e. the more expensive good will never have bigger chosen quantity 
                    data_for_power_test = base_data_for_power_test;
                    random_uniform_numbers = rand(obs_num, 1);
                    for obs = 1:obs_num
                        if max_quantities(obs, 1) > max_quantities(obs, 2)   % the 1st good is cheaper
                            single_random_uniform_choice_1 = safe_bundles(obs) + (1-random_uniform_numbers(obs)) * (max_quantities(obs, 1) - safe_bundles(obs));
                            single_random_uniform_choice_2 = random_uniform_numbers(obs) * safe_bundles(obs);
                        elseif max_quantities(obs, 1) < max_quantities(obs, 2)   % the 2nd good is cheaper
                            single_random_uniform_choice_1 = random_uniform_numbers(obs) * safe_bundles(obs);
                            single_random_uniform_choice_2 = safe_bundles(obs) + (1-random_uniform_numbers(obs)) * (max_quantities(obs, 2) - safe_bundles(obs));
                        else   % equal prices to both goods
                            single_random_uniform_choice_1 = max_quantities(obs,1) * random_uniform_numbers(obs);
                            single_random_uniform_choice_2 = max_quantities(obs,2) * (1-random_uniform_numbers(obs));
                        end
                        data_for_power_test(obs, 3:4) = [single_random_uniform_choice_1 , single_random_uniform_choice_2];
                    end
                end
                
                [~,~,~,~,~,~,~,~,~, VIO_PAIRS, VIOLATIONS, AFRIAT, VARIAN_Bounds, HoutmanMaks, MPI, ~] = ...
                                HPZ_Subject_Consistency (data_for_power_test, Graph_flags_for_power_test, GARP_flags_for_power_test, AFRIAT_flags_for_power_test, VARIAN_flags_for_power_test, HOUTMAN_flags_for_power_test, MPI_flags_for_power_test, 0, Varian_algorithm_settings, '', false, false, 0, 0);
                
                % assigning the results to the simulations list
                VIO_PAIRS_list(:,simulation) = VIO_PAIRS;
                VIOLATIONS_list(:,simulation) = VIOLATIONS;
                AFRIAT_list(:,simulation) = AFRIAT;
                VARIAN_Bounds_list(:,simulation) = VARIAN_Bounds;
                HoutmanMaks_list(:,simulation) = HoutmanMaks;
                MPI_list(:,simulation) = MPI;
                
            end
            
            % we ignore lower/upper bound, we only use the value that is the exact/approximation for Varian  
            VARIAN_list = VARIAN_Bounds_list([2,5,8] , :);
            
            % currently, "stats" contain median, mean, std (3 stats)
            num_stats = 3;
            VIO_PAIRS_stats = nan(3, num_stats);
            VIO_PAIRS_stats(:,1) = median(VIO_PAIRS_list , 2);
            VIO_PAIRS_stats(:,2) = mean(VIO_PAIRS_list , 2);
            VIO_PAIRS_stats(:,3) = std(VIO_PAIRS_list , [] , 2);
            VIOLATIONS_stats = nan(3, num_stats);
            VIOLATIONS_stats(:,1) = median(VIOLATIONS_list , 2);
            VIOLATIONS_stats(:,2) = mean(VIOLATIONS_list , 2);
            VIOLATIONS_stats(:,3) = std(VIOLATIONS_list , [] , 2);
            AFRIAT_stats = nan(1, num_stats);
            AFRIAT_stats(:,1) = median(AFRIAT_list , 2);
            AFRIAT_stats(:,2) = mean(AFRIAT_list , 2);
            AFRIAT_stats(:,3) = std(AFRIAT_list , [] , 2);
            VARIAN_stats = nan(3, num_stats);
            VARIAN_stats(:,1) = median(VARIAN_list , 2);
            VARIAN_stats(:,2) = mean(VARIAN_list , 2);
            VARIAN_stats(:,3) = std(VARIAN_list , [] , 2);
            HoutmanMaks_stats = nan(1, num_stats);
            HoutmanMaks_stats(:,1) = median(HoutmanMaks_list , 2);
            HoutmanMaks_stats(:,2) = mean(HoutmanMaks_list , 2);
            HoutmanMaks_stats(:,3) = std(HoutmanMaks_list , [] , 2);
%             MPI_stats = nan(2, num_stats);
%             MPI_stats(:,1) = median(MPI_list , 2);
%             MPI_stats(:,2) = mean(MPI_list , 2);
%             MPI_stats(:,3) = std(MPI_list , [] , 2);
            % we also check the number of simulations in which the randomized data was consistent with GARP   
            num_of_GARP_consistent_simulations = sum(~VIO_PAIRS_list(2,:));
            
            
            % print the results of the power test to the file
            fprintf(file_handle, '%s,%s,', ...
                num2str(num_of_simulations, precision_string), ...
                num2str(num_of_GARP_consistent_simulations, precision_string));
            if (GARP_flags(1) == 1)
                fprintf(file_handle, '%s,%s,%s, %s,%s,%s, %s,%s,%s, %s,%s,%s,', ...
                    num2str(VIOLATIONS_stats(1, 1), precision_string), ...
                    num2str(VIOLATIONS_stats(1, 2), precision_string), ...
                    num2str(VIOLATIONS_stats(1, 3), precision_string), ...
                    num2str(VIOLATIONS_stats(2, 1), precision_string), ...
                    num2str(VIOLATIONS_stats(2, 2), precision_string), ...
                    num2str(VIOLATIONS_stats(2, 3), precision_string), ...
                    num2str(VIO_PAIRS_stats(2, 1), precision_string), ...
                    num2str(VIO_PAIRS_stats(2, 2), precision_string), ...
                    num2str(VIO_PAIRS_stats(2, 3), precision_string), ...
                    num2str(VIOLATIONS_stats(3, 1), precision_string), ...
                    num2str(VIOLATIONS_stats(3, 2), precision_string), ...
                    num2str(VIOLATIONS_stats(3, 3), precision_string));
            end
            if (AFRIAT_flags(1) == 1)
                fprintf(file_handle, '%s,%s,%s,', ...
                    num2str(AFRIAT_stats(1, 1), precision_string), ...
                    num2str(AFRIAT_stats(1, 2), precision_string), ...
                    num2str(AFRIAT_stats(1, 3), precision_string));
            end
            if (VARIAN_flags(1) == 1)
                if (VARIAN_flags(5) == 1) % Mean Aggregator
                    fprintf(file_handle, '%s,%s,%s,', ...
                        num2str(VARIAN_stats(1, 1), precision_string), ...
                        num2str(VARIAN_stats(1, 2), precision_string), ...
                        num2str(VARIAN_stats(1, 3), precision_string));
                end
                if (VARIAN_flags(6) == 1) % AVG(SSQ) Aggregator
                    fprintf(file_handle, '%s,%s,%s,', ...
                        num2str(VARIAN_stats(2, 1), precision_string), ...
                        num2str(VARIAN_stats(2, 2), precision_string), ...
                        num2str(VARIAN_stats(2, 3), precision_string));
                end
                if (VARIAN_flags(7) == 1) % Maximum Aggregator
                    fprintf(file_handle, '%s,%s,%s,', ...
                        num2str(VARIAN_stats(3, 1), precision_string), ...
                        num2str(VARIAN_stats(3, 2), precision_string), ...
                        num2str(VARIAN_stats(3, 3), precision_string));
                end
            end
            if (HOUTMAN_flags(1) == 1)
                fprintf(file_handle, '%s,%s,%s,', ...
                    num2str(HoutmanMaks_stats(1, 1), precision_string), ...
                    num2str(HoutmanMaks_stats(1, 2), precision_string), ...
                    num2str(HoutmanMaks_stats(1, 3), precision_string));
            end
%             if (MPI_flags(1) == 1)
%                 fprintf(file_handle, '%s,%s,%s, %s,%s,%s,', ...
%                     num2str(MPI_stats(1, 1), precision_string), ...
%                     num2str(MPI_stats(1, 2), precision_string), ...
%                     num2str(MPI_stats(1, 3), precision_string), ...
%                     num2str(MPI_stats(2, 1), precision_string), ...
%                     num2str(MPI_stats(2, 2), precision_string), ...
%                     num2str(MPI_stats(2, 3), precision_string));
%             end
            
        
        else
            
            % we still don't implement power test for more than 2 goods
        
        end
        
    end
    
    
    
    % move to the next line in the main results file, for the next subject  
    fprintf(file_handle, '\n');
    
    
    
end   % end of loop over chosen subjects



% code for a single waitbar instead of multiple waitbars (part 3 out of 3) 
if (single_waitbar)
    % close the waitbar
    close(h_wb);
end



end