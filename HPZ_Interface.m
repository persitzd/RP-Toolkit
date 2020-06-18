function HPZ_Interface

% This file collects the information on the required estimation and reports
% the results to the user.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% before we start meddling with the warnings ('on' and 'off'), we first
% save the warning settings, and we call onCleanup in order to restore
% these settings when the function finishes or stops for any other reason
original_state = warning;
cleanup_warnings = onCleanup(@() warning(original_state));

% adding all the subfolders to the code path:
% supress warning. the warning is due to function 'assert' in the 'test'
% folder of 'matlab_bgl', that has the same name as MATLAB builtin
warning('off','MATLAB:dispatcher:nameConflict');
% determine where your m-file's folder is (we need it in other functions as well) 
main_folder = fileparts(which(mfilename));
% add that folder plus all subfolders to the path
addpath(genpath(main_folder));
% as default, the results files will be in a "Results" folder under the
% main folder which contains the code. but if you wish to have a different
% folder for the results, you can change the following line as needed
% ( e.g. - main_folder_for_results = pwd; )
main_folder_for_results = main_folder;



% set off warning when adding a sheet using xlswrite
warning('off','MATLAB:xlswrite:AddSheet');





% initialization of cell arrays (the definition as cells is required since 
% we have an option of multiple runs) 
data_matrix = cell(1,1);
subjects_index = cell(1,1);
action_flag = cell(1,1);
Graph_flags = cell(1,1);
power_test_settings = cell(1,1);
GARP_flags = cell(1,1);
AFRIAT_flags = cell(1,1); 
VARIAN_flags = cell(1,1);
HOUTMAN_flags = cell(1,1); 
MPI_flags = cell(1,1);
pref_class = cell(1,1);
choice_set_type = cell(1,1);
function_flag = cell(1,1); 
numeric_flag = cell(1,1); 
param1_restrictions = cell(1,1);
param2_restrictions = cell(1,1);
fix_corners = cell(1,1);
metric_flag = cell(1,1);
aggregation_flag = cell(1,1);
max_time_estimation = cell(1,1);
min_counter = cell(1,1); 
parallel_flag = cell(1,1); 
output_file_config = cell(1,1);
write_all_flag = cell(1,1); 
bootstrap_flag = cell(1,1);
file_val_str = cell(1,1);
residual_flag = cell(1,1);
in_sample_flag = cell(1,1);
out_sample_flag = cell(1,1);
bootstrap_sample_sizes = cell(1,1);
bootstrap_significance_level = cell(1,1);
BI_threshold = cell(1,1);
max_starting_points = cell(1,1);
one_residuals_file = cell(1,1);
debugger_mode = cell(1,1);
waitbar_settings = cell(1,1);
Varian_algorithm_settings = cell(1,1);
print_precision = cell(1,1);

% this loop is needed for multiple runs (each round of the loop = one run) 
another_run = true;
runs_counter = 0;
while (another_run)
    
    % updating number of runs
    runs_counter = runs_counter + 1;
    % i is just for easier use
    i = runs_counter;
    
    % (1) Showing the user the user-interface screens one after another,
    %     and saving all the user's decisions and preferences in these screens
    % (2) In the last screen, the user is asked whether he desires to make
    %     more runs (with a different data set / method / etc.), or whther
    %     he is finished and the program should start running (another_run==false) 
    [ok, data_matrix{i}, subjects_index{i}, action_flag{i}, Graph_flags{i}, power_test_settings{i}, GARP_flags{i}, AFRIAT_flags{i}, VARIAN_flags{i}, HOUTMAN_flags{i}, MPI_flags{i}, pref_class{i}, choice_set_type{i}, function_flag{i}, numeric_flag{i}, param1_restrictions{i}, param2_restrictions{i}, fix_corners{i}, metric_flag{i}, aggregation_flag{i}, max_time_estimation{i}, min_counter{i}, parallel_flag{i}, output_file_config{i}, write_all_flag{i}, bootstrap_flag{i}, file_val_str{i}, residual_flag{i}, in_sample_flag{i}, out_sample_flag{i}, bootstrap_sample_sizes{i}, bootstrap_significance_level{i}, BI_threshold{i}, max_starting_points{i}, one_residuals_file{i}, debugger_mode{i}, waitbar_settings{i}, Varian_algorithm_settings{i}] = HPZ_All_Screens_Manager(main_folder, runs_counter);
    % if the user pressed "cancel" in one of the screens, then this run was cancelled 
    if (ok == 0)
        runs_counter = runs_counter - 1;
    end
    
    if (runs_counter == 0)
        % then the user pressed "cancel" while in the first run - so we just end the program  
        return
    end
    
    [ok, another_run] = HPZ_Screen_Another_Run(runs_counter);
    % if the user pressed "cancel All" - we just end the program
    if (ok == 0)
        return
    end
    
end



% this loop performs all the runs one after another
for i = 1:runs_counter
    
    start_time = now;   % we want to tell the user how long the whole run took
    
    print_precision{i} = HPZ_Constants.print_precision;
    
    try
        % All calculations and printing to results file occur here
        if action_flag{i} == HPZ_Constants.Consistency_action
            % Consistency Indices
            HPZ_Consistency_Indices_Manager(data_matrix{i}, subjects_index{i}, choice_set_type{i}, Graph_flags{i}, power_test_settings{i}, GARP_flags{i}, AFRIAT_flags{i}, VARIAN_flags{i}, HOUTMAN_flags{i}, MPI_flags{i}, max_time_estimation{i}, print_precision{i}, one_residuals_file{i}, waitbar_settings{i}, Varian_algorithm_settings{i}, main_folder_for_results, i, runs_counter);  

        elseif (action_flag{i} == HPZ_Constants.NLLS_action) || (action_flag{i} == HPZ_Constants.MMI_action) || (action_flag{i} == HPZ_Constants.BI_action)
            % Parameters Estimation (NLLS / MMI / BI)
            HPZ_Estimation_Manager(data_matrix{i}, subjects_index{i}, choice_set_type{i}, action_flag{i}, pref_class{i}, function_flag{i}, numeric_flag{i}, param1_restrictions{i}, param2_restrictions{i}, fix_corners{i}, metric_flag{i}, aggregation_flag{i}, max_time_estimation{i}, min_counter{i}, parallel_flag{i}, output_file_config{i}, write_all_flag{i}, bootstrap_flag{i}, file_val_str{i}, residual_flag{i}, in_sample_flag{i}, out_sample_flag{i}, print_precision{i}, bootstrap_sample_sizes{i}, bootstrap_significance_level{i}, BI_threshold{i}, max_starting_points{i}, one_residuals_file{i}, debugger_mode{i}, waitbar_settings{i}, main_folder_for_results, i, runs_counter);

        end
    catch error_var
        HPZ_Print_Error_Message(error_var, main_folder_for_results, i, runs_counter);
    end
    
    % tell the user how long the whole run took
    end_time = now;
    running_time = datevec(end_time - start_time);
    days = running_time(3);
    hours = running_time(4);
    minutes = running_time(5);
    seconds = running_time(6);
    if days
        fprintf('\nRun number %d running time was:\n%d days, %d hours, %d minutes and %.0f seconds.\n\n',...
                                                                        i, days, hours, minutes, seconds);
    else
        fprintf('\nRun number %d running time was:\n%d hours, %d minutes and %.3f seconds.\n\n',...
                                                                        i, hours, minutes, seconds);
    end
end



% Address the folder where the result files were saved
msgbox(char(strcat('All files are saved under the following path:', {' '}, main_folder_for_results, '\', HPZ_Constants.results_files_dir)), 'Output File Path','help');
% Opening the folder (works only in Windows)
if strcmp(computer('arch'),'win64') || strcmp(computer('arch'),'win32')
    winopen(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir));
end



end