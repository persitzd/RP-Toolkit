function HPZ_Optimization_Settings_Write(aggregation_flag, metric_flag, max_time_estimation, min_counter, parallel_flag, main_folder)

% This function is called when the user wishes (and it is also the
% default...) to save the choices he made now for the next time.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% preparations for saving the current decisions and settings of the user 
if max_time_estimation == inf
    max_time_estimation = 0;
end
if min_counter == inf
    min_counter = HPZ_Constants.infinity;
end



% make sure the Settings Files directory exists, if not - create it
dir_exists = exist(strcat(main_folder, '/', HPZ_Constants.settings_files_dir) , 'dir');
if ~dir_exists
    mkdir(HPZ_Constants.settings_files_dir);
end

% creating/opening the settings file
settings_file = fopen(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.optimization_settings_file_name, '.csv'), 'w');

% close file when this function end
cleanup_close_file = onCleanup(@() fclose(settings_file));



%% printing the new settings and choices

% aggregation_flag for MMI (Max / Mean / AVGSSQ)
fprintf(settings_file, '%s\n', num2str(aggregation_flag));   % line (1)
% metric_flag for NLLS (euclidean or CFGK)
fprintf(settings_file, '%s\n', num2str(metric_flag));   % line (2)
% max_time_estimation in minutes ('0' if there is no limit)
if isempty(max_time_estimation)
    max_time_estimation = 0;
end
fprintf(settings_file, '%s\n', num2str(max_time_estimation));   % line (3)
% min_counter
fprintf(settings_file, '%s\n', num2str(min_counter));   % line (4)
% parallel_flag
fprintf(settings_file, '%s\n', num2str(parallel_flag));   % line (5)



end