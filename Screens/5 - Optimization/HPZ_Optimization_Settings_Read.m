function [aggregation_flag, metric_flag, max_time_estimation, min_counter, parallel_flag] = HPZ_Optimization_Settings_Read(main_folder)

% this function reads the user-interface screens saved settings

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% read settings from file 
try
    % reading from the file
    [settings] = csvread(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.optimization_settings_file_name, '.csv'));
    
    % in the next line, we call the 5th row and 1st column in the settings
    % matrix, because we have 5 rows and the max number of columns is 1.
    % if it doesn't exists, it means that the settings file is corrupted,
    % cause it doesn't have all required data.
    % if in the future you will add more rows (>5) or use longer rows (num
    % of columns > 1), you should increase these numbers respectively.
    check_settings = settings(5 , 1); %#ok<NASGU>
catch
    % A) if we failed to read the file,
    % we assume it might be because it was accidently deleted, or because
    % it is not formatted in a way that allows it to be read.
    % B) if the file's data was corrupted, we assume someone might've
    % played with it and made changes to it by mistake.
    % A+B) either way, we reset the file, then read it again.
    HPZ_Optimization_Settings_Reset(main_folder);
    [settings] = csvread(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.optimization_settings_file_name, '.csv'));
end

% assign initial values to variables in accordance with the settings
aggregation_flag = settings(1,1); %
metric_flag = settings(2,1); %
max_time_estimation = settings(3,1);
if max_time_estimation <= 0
    max_time_estimation = inf;
end
min_counter = settings(4,1); %
if min_counter == HPZ_Constants.infinity
    min_counter = inf;
end
parallel_flag = settings(5,1);



%% for (almost) every variable, we check if it has a valid value, otherwise we reset it

% aggregation and metric flags
if (aggregation_flag ~= HPZ_Constants.MMI_Max && aggregation_flag ~= HPZ_Constants.MMI_Mean && aggregation_flag ~= HPZ_Constants.MMI_AVGSSQ)
   aggregation_flag = HPZ_Constants.MMI_AVGSSQ;
end
if (metric_flag ~= HPZ_Constants.euclidean_metric && metric_flag ~= HPZ_Constants.CFGK_metric && metric_flag ~= HPZ_Constants.normalized_euclidean_metric)
   metric_flag = HPZ_Constants.euclidean_metric;
end
% parallel flag
if (parallel_flag ~= 0 && parallel_flag ~= 1)
    parallel_flag = 1;
end



end