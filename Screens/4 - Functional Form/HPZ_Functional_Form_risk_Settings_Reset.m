function HPZ_Functional_Form_risk_Settings_Reset(main_folder)

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
settings_file = fopen(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.functional_form_risk_settings_file_name, '.csv'), 'w');

% close file when this function end
cleanup_close_file = onCleanup(@() fclose(settings_file));



%% printing the default settings, which are:

% not numeric (analytic is default for risk)
fprintf(settings_file, '%s,%s\n', '0','0');   % line (1)
% functional form (set to first for risk)
fprintf(settings_file, '%s,%s\n', '1','1');   % line (2)
% restrictions on the 1st parameter (beta) in risk preferences
fprintf(settings_file, '%s,%s\n', '-1', num2str(HPZ_Constants.infinity));   % line (3)
% restrictions on the 2nd parameter (rho/A) in risk preferences
fprintf(settings_file, '%s,%s\n', '0', num2str(HPZ_Constants.infinity));   % line (4)
% fix corners (CFGK) or not - default is not
fprintf(settings_file, '%s\n', '0');   % line (5)



end

