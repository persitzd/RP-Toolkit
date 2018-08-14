function [file_handle, file_name_str, file_residuals_path_str, file_warnings_path_str] = HPZ_Write_Result_File_Initializer(output_file_config, action_flag, pref_class, function_flag, criterion_value_str, bootstrap_flag, significance_level, main_folder_for_results)

% this function creates the main results file for parameters estimation
% results,
% and prints the column headers to the file.

% more details about the formatting of this file can be found in the
% documentation of HPZ_Write_Result_File_Finalizer.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



%% the names of the 2 parameters that are being estimated

if pref_class == HPZ_Constants.risk_pref  && function_flag == HPZ_Constants.CRRA_func  % CRRA 
    param_1 = 'Beta';
    param_2 = 'Rho';
elseif pref_class == HPZ_Constants.risk_pref  && function_flag == HPZ_Constants.CARA_func   % CARA
    param_1 = 'Beta';
    param_2 = 'A';
elseif pref_class == HPZ_Constants.OR_pref && function_flag == HPZ_Constants.CES_func  % CES
    param_1 = 'Alpha';
    param_2 = 'Rho';
end



%% determining the names of the files, and opening (creating) the files

% NLLS
if action_flag == HPZ_Constants.NLLS_action
    action_str = HPZ_Constants.NLLS_action_file_name;
% MMI
elseif action_flag == HPZ_Constants.MMI_action
    action_str = HPZ_Constants.MMI_action_file_name;
% BI
elseif action_flag == HPZ_Constants.BI_action
    action_str = HPZ_Constants.BI_action_file_name;     
end

% make sure the Results directory exists, if not - create it
dir_exists = exist(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir) , 'dir');
if ~dir_exists
    mkdir(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir));
end

% name of main results file - we make sure that there isn't an existing file by the same name    
while (true)
    file_name_str = char(strcat(action_str, {' - '}, 'Results - Date', {' '}, date, {' - '}, generate_random_str));
    file_path_str = strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir, '/', file_name_str, '.csv');
    if ~exist(file_path_str, 'file')
        break
    end
end
% opening the main results file
file_handle = fopen(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir, '/', file_name_str, '.csv'), 'w'); 

% name of residuals file (if needed) - we make sure that there isn't an existing file by the same name    
while (true)
    file_residuals_str = char(strcat(action_str, {' - '}, 'Results Residuals - Date', {' '}, date, {' - '}, generate_random_str));
    file_residuals_path_str = strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir, '/', file_residuals_str);
    if ~exist(strcat(file_residuals_path_str, '.csv'), 'file')
        break
    end
end

% name of warnings file (if needed) - we make sure that there isn't an existing file by the same name    
while (true)
    file_warnings_str = char(strcat(action_str, {' - '}, 'Warnings - Date', {' '}, date, {' - '}, generate_random_str));
    file_warnings_path_str = strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir, '/', file_warnings_str);
    if ~exist(strcat(file_warnings_path_str, '.csv'), 'file')
        break
    end
end





%% printing the headers of each column in the results file

% there are at least three strings - the first is the subject's index,
% and the second and third are the estimated parameters
fprintf(file_handle, '%s,%s,%s', 'Subject', param_1, param_2);

% if bootstrap was asked for, we need to add 8 more columns
if (bootstrap_flag == true)
    
    confidence_level_percent = num2str((1-significance_level)*100);
    
    fprintf(file_handle, ',%s,%s,%s,%s,%s,%s,%s,%s', ...
            char(strcat({'Simulated Mean '}, param_1)), ...
            char(strcat({'Simulated Mean '}, param_2)), ...
            char(strcat({'Simulated SD '}, param_1)), ...
            char(strcat({'Simulated SD '}, param_2)), ...
            char(strcat({'Simulated Low '}, confidence_level_percent, {'% '}, param_1)), ...
            char(strcat({'Simulated Low '}, confidence_level_percent, {'% '}, param_2)), ...
            char(strcat({'Simulated High '}, confidence_level_percent, {'% '}, param_1)), ...
            char(strcat({'Simulated High '}, confidence_level_percent, {'% '}, param_2)));
end

% for each of the criterions, we check whether we were asked to print
% that criterion or not
% reminder: the criterions are:
%   1. NLLS - Euclidean metric
%   2. NLLS - CFGK metric
%   3. NLLS - normalized-Euclidean metric
%   4. MMI - Max aggregator
%   5. MMI - Mean aggregator
%   6. MMI - Average Sum of Squares aggregator
%   7. BI

full_str = criterion_value_str;
full_str{1} = char(strcat({'NLLS '}, criterion_value_str{1}, {' Euclidean Metric'}));
full_str{2} = char(strcat({'NLLS '}, criterion_value_str{2}, {' CFGK Metric'}));
full_str{3} = char(strcat({'NLLS '}, criterion_value_str{3}, {' normalized-Euclidean Metric'}));
full_str{4} = char(strcat({'MMI '}, criterion_value_str{4}, {' Max(Wastes)'}));
full_str{5} = char(strcat({'MMI '}, criterion_value_str{5}, {' Mean(Wastes)'}));
full_str{6} = char(strcat({'MMI '}, criterion_value_str{6}, {' Avg(SumOfSquares(Wastes))'}));
full_str{7} = char(strcat({'BI '}, criterion_value_str{7}));

for i=1:length(output_file_config)
    if (output_file_config(i) == 1)
        fprintf(file_handle, ',%s', full_str{i});
    end
end

% after printing the headers line, we go down to the next line when
% we start to print the results for the first subject
fprintf(file_handle, '\n');

end