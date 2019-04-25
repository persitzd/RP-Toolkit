function [output_file_config, write_all_flag, bootstrap_flag] = HPZ_Output_Format_Settings_Read(main_folder)

% this function reads the user-interface screens saved settings

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% read settings from file 
try
    % reading from the file
    [settings] = csvread(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.output_format_settings_file_name, '.csv'));
    
    % in the next line, we call the 3th row and 8th column in the settings
    % matrix, because we have 3 rows and the max number of columns is 7
    % (it is 7 for output file config).
    % if it doesn't exists, it means that the settings file is corrupted,
    % cause it doesn't have all required data.
    % if in the future you will add more rows (>3) or use longer rows (num
    % of columns > 7), you should increase these numbers respectively.
    check_settings = settings(3 , 7); %#ok<NASGU>
catch
    % A) if we failed to read the file,
    % we assume it might be because it was accidently deleted, or because
    % it is not formatted in a way that allows it to be read.
    % B) if the file's data was corrupted, we assume someone might've
    % played with it and made changes to it by mistake.
    % A+B) either way, we reset the file, then read it again.
    HPZ_Output_Format_Settings_Reset(main_folder);
    [settings] = csvread(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.output_format_settings_file_name, '.csv'));
end

% assign initial values to variables in accordance with the settings
output_file_config = settings(1,1:7);
write_all_flag = settings(2,1);
bootstrap_flag = settings(3,1);



%% for (almost) every variable, we check if it has a valid value, otherwise we reset it
% output-related flags
for i=1:length(output_file_config)
    if (output_file_config(i) ~= 0 && output_file_config(i) ~= 1)
        output_file_config(i) = 1;
    end
end
if (write_all_flag ~= 0 && write_all_flag ~= 1)
    write_all_flag = 1;
end
% bootstrap flag
if (bootstrap_flag ~= 0 && bootstrap_flag ~= 1)
    bootstrap_flag = 1;
end



end