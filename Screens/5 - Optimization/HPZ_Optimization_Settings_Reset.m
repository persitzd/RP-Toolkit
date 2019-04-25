function HPZ_Optimization_Settings_Reset(main_folder)

% This function is used to restore the Settings file, in case it was
% lost/deleted for some reason, or if it was corrupted.
% This function is called when an error occurs while attempting to read Settings.csv, 
% but it can also run directly if you desire to reset all the settings.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% make sure the Settings Files directory exists, if not - create it
dir_exists = exist(strcat(main_folder, '/', HPZ_Constants.settings_files_dir) , 'dir');
if ~dir_exists
    mkdir(strcat(main_folder, '/', HPZ_Constants.settings_files_dir));
end

% creating/opening the settings file
settings_file = fopen(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.optimization_settings_file_name, '.csv'), 'w');

% close file when this function end
cleanup_close_file = onCleanup(@() fclose(settings_file));



%% printing the default settings, which are:

% aggregation_flag for MMI (Max / Mean / AVGSSQ) default is AVGSSQ
fprintf(settings_file, '%s\n', '3');   % line (1)
% metric_flag for NLLS (euclidean or CFGK) default is euclidean
fprintf(settings_file, '%s\n', '1');   % line (2)
% max_time_estimation in minutes ('0' if there is no limit)
fprintf(settings_file, '%s\n', '0');   % line (3)
% min_counter default is the highest possible except for inf
fprintf(settings_file, '%s\n', '100');   % line (4)
% parallel_flag default is false
fprintf(settings_file, '%s\n', '0');   % line (5)



end

