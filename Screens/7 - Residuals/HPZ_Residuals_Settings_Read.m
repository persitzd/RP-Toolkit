function [residual_flag, in_sample_flag, out_sample_flag] = HPZ_Residuals_Settings_Read(main_folder)

% this function reads the user-interface screens saved settings

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% read settings from file 
try
    % reading from the file
    [settings] = csvread(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.residuals_settings_file_name, '.csv'));
    
    % in the next line, we call the 1st row and 3rd column in the settings
    % matrix, because we have 1 rows and the max number of columns is 3.
    % if it doesn't exists, it means that the settings file is corrupted,
    % cause it doesn't have all required data.
    % if in the future you will add more rows (>1) or use longer rows (num
    % of columns > 3), you should increase these numbers respectively.
    check_settings = settings(1 , 3); %#ok<NASGU>
catch
    % A) if we failed to read the file,
    % we assume it might be because it was accidently deleted, or because
    % it is not formatted in a way that allows it to be read.
    % B) if the file's data was corrupted, we assume someone might've
    % played with it and made changes to it by mistake.
    % A+B) either way, we reset the file, then read it again.
    HPZ_Residuals_Settings_Reset(main_folder);
    [settings] = csvread(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.residuals_settings_file_name, '.csv'));
end

% assign initial values to variables in accordance with the settings
residual_flag = settings(1,1);
in_sample_flag = settings(1,2);
out_sample_flag = settings(1,3);



%% for (almost) every variable, we check if it has a valid value, otherwise we reset it

% residual flags
if (residual_flag ~= 0 && residual_flag ~= 1)
    residual_flag = 1;
end
if (in_sample_flag ~= 0 && in_sample_flag ~= 1)
    in_sample_flag = 1;
end
if (out_sample_flag ~= 0 && out_sample_flag ~= 1)
    out_sample_flag = 1;
end



end