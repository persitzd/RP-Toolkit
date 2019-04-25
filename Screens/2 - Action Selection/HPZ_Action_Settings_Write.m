function HPZ_Action_Settings_Write(action_choice, main_folder)

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
settings_file = fopen(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.action_settings_file_name, '.csv'), 'w');

% close file when this function end
cleanup_close_file = onCleanup(@() fclose(settings_file));



%% printing the new settings and choices

% action
fprintf(settings_file, '%s\n', num2str(action_choice));   % line (1)



end