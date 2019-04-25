function [action_choice] = HPZ_Action_Settings_Read(main_folder)

% this function reads the user-interface screens saved settings

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% read settings from file 
try
    % reading from the file
    [settings] = csvread(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.action_settings_file_name, '.csv'));
    
    % in the next line, we call the 1st row and 1st column in the settings
    % matrix, because we have 24 rows and the max number of columns is 1.
    % if it doesn't exists, it means that the settings file is corrupted,
    % cause it doesn't have all required data.
    % if in the future you will add more rows (>1) or use longer rows (num
    % of columns > 1), you should increase these numbers respectively.
    check_settings = settings(1 , 1); %#ok<NASGU>
catch
    % A) if we failed to read the file,
    % we assume it might be because it was accidently deleted, or because
    % it is not formatted in a way that allows it to be read.
    % B) if the file's data was corrupted, we assume someone might've
    % played with it and made changes to it by mistake.
    % A+B) either way, we reset the file, then read it again.
    HPZ_Action_Settings_Reset(main_folder);
    [settings] = csvread(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.action_settings_file_name, '.csv'));
end

action_choice = settings(1,1);



%% for (almost) every variable, we check if it has a valid value, otherwise we reset it

% action flag
if (action_choice ~= HPZ_Constants.Consistency_action && action_choice ~= HPZ_Constants.NLLS_action && action_choice ~= HPZ_Constants.MMI_action && action_choice ~= HPZ_Constants.BI_action)
    action_choice = HPZ_Constants.Consistency_action;
end



end