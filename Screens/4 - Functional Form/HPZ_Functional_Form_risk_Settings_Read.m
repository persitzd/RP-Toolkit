function [numeric_flag, function_flag, param1_restrictions, param2_restrictions, fix_corners] = HPZ_Functional_Form_risk_Settings_Read(main_folder)

% this function reads the user-interface screens saved settings

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% read settings from file 
try
    % reading from the file
    [settings] = csvread(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.functional_form_risk_settings_file_name, '.csv'));
    
    % in the next line, we call the 5th row and 2nd column in the settings
    % matrix, because we have 5 rows and the max number of columns is 2
    % (it is 2 for param restrictions).
    % if it doesn't exists, it means that the settings file is corrupted,
    % cause it doesn't have all required data.
    % if in the future you will add more rows (>5) or use longer rows (num
    % of columns > 2), you should increase these numbers respectively.
    check_settings = settings(5 , 2); %#ok<NASGU>
catch
    % A) if we failed to read the file,
    % we assume it might be because it was accidently deleted, or because
    % it is not formatted in a way that allows it to be read.
    % B) if the file's data was corrupted, we assume someone might've
    % played with it and made changes to it by mistake.
    % A+B) either way, we reset the file, then read it again.
    HPZ_Functional_Form_risk_Settings_Reset(main_folder);
    [settings] = csvread(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.functional_form_risk_settings_file_name, '.csv'));
end

% assign initial values to variables in accordance with the settings
numeric_flag = settings(1,1);
function_flag = settings(2,1);
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
fix_corners = settings(5,1);




%% for (almost) every variable, we check if it has a valid value, otherwise we reset it

% numeric flag
if (numeric_flag ~= HPZ_Constants.analytic && numeric_flag ~= HPZ_Constants.numeric && numeric_flag ~= HPZ_Constants.semi_numeric)
    numeric_flag = HPZ_Constants.analytic;
end
% function flag
if (function_flag ~= HPZ_Constants.CRRA_func && function_flag ~= HPZ_Constants.CARA_func)
   function_flag = HPZ_Constants.CRRA_func;
end
% fix corners
if (fix_corners ~= 0 && fix_corners ~= 1)
    fix_corners = 1;
end



end