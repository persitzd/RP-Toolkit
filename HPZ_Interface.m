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



% what type of waitbar to show (none / single per run / seperate per subject) 
waitbar_options = HPZ_Constants.waitbar_per_subject;    % show a seperate waitbar for each subject 
%waitbar_options = HPZ_Constants.waitbar_single;        % show a single waitbar for all subjects 
%waitbar_options = HPZ_Constants.waitbar_none;          % do not show any waitbar  



% If 1, the function doesn't print the results, if 0 it prints the results
% normally.
% It is normally set to 0, 
% but while residuals are being calculated, it is set to 1
%global do_not_write

% set off warning when adding a sheet using xlswrite
warning('off','MATLAB:xlswrite:AddSheet');





% initialization of cell arrays (the definition as cells is required since 
% we have an option of multiple runs) 
data_matrix = cell(1,1);
subjects_index = cell(1,1);
action_flag = cell(1,1);
GARP_flags = cell(1,1);
AFRIAT_flags = cell(1,1); 
VARIAN_flags = cell(1,1);
HOUTMAN_flags = cell(1,1); 
pref_class = cell(1,1);
function_flag = cell(1,1); 
numeric_flag = cell(1,1); 
param1_restrictions = cell(1,1);
param2_restrictions = cell(1,1);
fix_corners = cell(1,1);
metric_flag = cell(1,1);
aggregation_flag = cell(1,1);
max_time_estimation = cell(1,1);
max_starting_points = cell(1,1);
min_counter = cell(1,1); 
parallel_flag = cell(1,1); 
output_file_config = cell(1,1);
write_all_flag = cell(1,1); 
bootstrap_flag = cell(1,1);
file_val_str = cell(1,1);
residual_flag = cell(1,1);
in_sample_flag = cell(1,1);
out_sample_flag = cell(1,1);
one_residuals_file = cell(1,1);
significance_level = cell(1,1);
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
    [ok, data_matrix{i}, subjects_index{i}, action_flag{i}, GARP_flags{i}, AFRIAT_flags{i}, VARIAN_flags{i}, HOUTMAN_flags{i}, pref_class{i}, function_flag{i}, numeric_flag{i}, param1_restrictions{i}, param2_restrictions{i}, fix_corners{i}, metric_flag{i}, aggregation_flag{i}, max_time_estimation{i}, max_starting_points{i}, min_counter{i}, parallel_flag{i}, output_file_config{i}, write_all_flag{i}, bootstrap_flag{i}, file_val_str{i}, residual_flag{i}, in_sample_flag{i}, out_sample_flag{i}] = HPZ_All_Screens_Manager(main_folder, runs_counter);
    %HPZ_All_Screens_Manager(runs_counter);
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
    
    one_residuals_file{i} = true;   % change this to false in order to print the residuals as a separate file for each subject, instead of one file with separate sheets
    significance_level{i} = HPZ_Constants.significance_level;
    print_precision{i} = HPZ_Constants.print_precision;
    
    try
        % All calculations and printing to results file occur here
        if action_flag{i} == HPZ_Constants.Consistency_action
            % Consistency Indexes
            HPZ_Consistency_Indices_Manager(data_matrix{i}, subjects_index{i}, GARP_flags{i}, AFRIAT_flags{i}, VARIAN_flags{i}, HOUTMAN_flags{i}, one_residuals_file{i}, max_time_estimation{i}, print_precision{i}, waitbar_options, main_folder_for_results, i, runs_counter);  

        elseif (action_flag{i} == HPZ_Constants.NLLS_action) || (action_flag{i} == HPZ_Constants.MMI_action) || (action_flag{i} == HPZ_Constants.BI_action)
            % Parameters Estimation (NLLS / MMI / BI)
            HPZ_Estimation_Manager(data_matrix{i}, subjects_index{i}, action_flag{i}, pref_class{i}, function_flag{i}, numeric_flag{i}, param1_restrictions{i}, param2_restrictions{i}, fix_corners{i}, metric_flag{i}, aggregation_flag{i}, max_time_estimation{i}, min_counter{i}, max_starting_points{i}, parallel_flag{i}, output_file_config{i}, write_all_flag{i}, bootstrap_flag{i}, file_val_str{i}, residual_flag{i}, in_sample_flag{i}, out_sample_flag{i}, one_residuals_file{i}, significance_level{i}, print_precision{i}, waitbar_options, main_folder_for_results, i, runs_counter);

        end
    catch error_var
        HPZ_Print_Error_Message(error_var, main_folder_for_results, i, runs_counter);
    end
end



% Address the folder where the result files were saved
msgbox(char(strcat('All files are saved under the following path:', {' '}, main_folder_for_results, '\', HPZ_Constants.results_files_dir)), 'Output File Path','help');
% Opening the folder (works only in Windows)
if strcmp(computer('arch'),'win64') || strcmp(computer('arch'),'win32')
    winopen(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir));
end



end