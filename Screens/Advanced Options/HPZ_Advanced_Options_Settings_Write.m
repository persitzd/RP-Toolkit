function HPZ_Advanced_Options_Settings_Write(bootstrap_sample_sizes, bootstrap_significance_level, BI_threshold, max_starting_points, possible_num_convergence_points, one_residuals_file, debugger_mode, waitbar_settings, Varian_algorithm_settings, main_folder)

% This function is called when the user wishes (and it is also the
% default...) to save the choices he made now for the next time.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% make sure the Settings Files directory exists, if not - create it
dir_exists = exist(strcat(main_folder, '/', HPZ_Constants.settings_files_dir) , 'dir');
if ~dir_exists
    mkdir(HPZ_Constants.settings_files_dir);
end

% creating/opening the settings file
settings_file = fopen(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.advanced_options_settings_file_name, '.csv'), 'w');

% close file when this function end
cleanup_close_file = onCleanup(@() fclose(settings_file));



%% printing the new settings and choices

% bootstrap_sample_sizes (for analytic, numeric, semi-numeric)
fprintf(settings_file, '%s,%s,%s\n', num2str(bootstrap_sample_sizes(1)), ...
                    num2str(bootstrap_sample_sizes(2)), num2str(bootstrap_sample_sizes(3)));   % line (1)
% bootstrap_significance_level
fprintf(settings_file, '%s\n', num2str(bootstrap_significance_level));   % line (2)

% BI_threshold
fprintf(settings_file, '%s\n', num2str(BI_threshold));   % line (3)

% max_starting_points
fprintf(settings_file, '%s\n', num2str(max_starting_points));   % line (4)
% possible_num_convergence_points
for i=1:length(possible_num_convergence_points)
    fprintf(settings_file, '%s,', num2str(possible_num_convergence_points(i)));   % line (5)
end
fprintf(settings_file, '\n');

% one_residuals_file (1 = one file for all subjects , 0 = one file per subject) 
fprintf(settings_file, '%s\n', num2str(one_residuals_file));   % line (6)

% debugger_mode (1 = activated , 0 = not activated)
fprintf(settings_file, '%s\n', num2str(debugger_mode));   % line (7)

% waitbar_settings (minimum num of subjects for single waitbar, show
% residuals waitbar, show bootstrap waitbar)
fprintf(settings_file, '%s,%s,%s,%s\n', num2str(waitbar_settings(1)), num2str(waitbar_settings(2)), ...
                                        num2str(waitbar_settings(3)), num2str(waitbar_settings(4)));   % line (8)

% Varian_algorithm_settings (max number of observations in new problem,
% minimum multiply of num of  new observations compared to original num of observations 
fprintf(settings_file, '%s,%s\n', num2str(Varian_algorithm_settings(1)), ...
                                num2str(Varian_algorithm_settings(2)));   % line (9)



end