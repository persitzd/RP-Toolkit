function HPZ_Functional_Form_risk_Settings_Write(numeric_flag, function_flag, param1_restrictions, param2_restrictions, fix_corners, main_folder)

% This function is called when the user wishes (and it is also the
% default...) to save the choices he made now for the next time.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% preparations for saving the current decisions and settings of the user 
if (param1_restrictions(1) == -inf) 
    param1_restrictions(1) = -HPZ_Constants.infinity; end
if (param1_restrictions(2) == inf) 
    param1_restrictions(2) = HPZ_Constants.infinity; end
if (param2_restrictions(1) == -inf) 
    param2_restrictions(1) = -HPZ_Constants.infinity; end
if (param2_restrictions(2) == inf) 
    param2_restrictions(2) = HPZ_Constants.infinity; end



% make sure the Settings Files directory exists, if not - create it
dir_exists = exist(strcat(main_folder, '/', HPZ_Constants.settings_files_dir) , 'dir');
if ~dir_exists
    mkdir(HPZ_Constants.settings_files_dir);
end

% creating/opening the settings file
settings_file = fopen(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.functional_form_risk_settings_file_name, '.csv'), 'w');

% close file when this function end
cleanup_close_file = onCleanup(@() fclose(settings_file));



%% printing the new settings and choices
% numeric (risk)
fprintf(settings_file, '%s\n', num2str(numeric_flag));   % line (1)
% functional form (risk)
fprintf(settings_file, '%s\n', num2str(function_flag));   % line (2)
% restrictions on the 1st parameter (beta) in risk preferences
fprintf(settings_file, '%s,%s\n', num2str(param1_restrictions(1)), num2str(param1_restrictions(2)));   % line (3)
% restrictions on the 2nd parameter (rho/A) in risk preferences
fprintf(settings_file, '%s,%s\n', num2str(param2_restrictions(1)), num2str(param2_restrictions(2)));   % line (4)
% fix corners (CFGK) or not
fprintf(settings_file, '%s\n', num2str(fix_corners));   % line (5)



end