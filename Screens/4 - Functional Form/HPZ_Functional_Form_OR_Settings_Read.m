function [numeric_flag, function_flag, param1_restrictions, param2_restrictions] = HPZ_Functional_Form_OR_Settings_Read(main_folder)

% this function reads the user-interface screens saved settings

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% read settings from file 
try
    % reading from the file
    [settings] = csvread(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.functional_form_OR_settings_file_name, '.csv'));
    
    % in the next line, we call the 4th row and 2nd column in the settings
    % matrix, because we have 4 rows and the max number of columns is 7
    % (it is 2 for param restrictions).
    % if it doesn't exists, it means that the settings file is corrupted,
    % cause it doesn't have all required data.
    % if in the future you will add more rows (>4) or use longer rows (num
    % of columns > 2), you should increase these numbers respectively.
    check_settings = settings(4 , 2); %#ok<NASGU>
catch
    % A) if we failed to read the file,
    % we assume it might be because it was accidently deleted, or because
    % it is not formatted in a way that allows it to be read.
    % B) if the file's data was corrupted, we assume someone might've
    % played with it and made changes to it by mistake.
    % A+B) either way, we reset the file, then read it again.
    HPZ_Functional_Form_OR_Settings_Reset(main_folder);
    [settings] = csvread(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.functional_form_OR_settings_file_name, '.csv'));
end

% assign initial values to variables in accordance with the settings
numeric_flag = settings(1,1);
function_flag = settings(2,1); %
param1_restrictions = settings(3,1:2);
param2_restrictions = settings(4,1:2);
if (param1_restrictions(1) == -HPZ_Constants.infinity) 
    param1_restrictions(1) = -inf; end
if (param1_restrictions(2) == HPZ_Constants.infinity) 
    param1_restrictions(2) = inf; end
if (param2_restrictions(1) == -HPZ_Constants.infinity) 
    param2_restrictions(1) = -inf; end
if (param2_restrictions(2) == HPZ_Constants.infinity) 
    param2_restrictions(2) = inf; end



%% for (almost) every variable, we check if it has a valid value, otherwise we reset it

% numeric flag
if (numeric_flag ~= HPZ_Constants.analytic && numeric_flag ~= HPZ_Constants.numeric && numeric_flag ~= HPZ_Constants.semi_numeric)
    numeric_flag = HPZ_Constants.analytic;
end
% function flag
if (function_flag ~= HPZ_Constants.CES_func)
   function_flag = HPZ_Constants.CES_func;
end



end