function [save_user_choices, fix_endowments, action_choice, GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, risk_numeric_flag, OR_numeric_flag, risk_function_flag, OR_function_flag, risk_param1_restrictions, risk_param2_restrictions, OR_param1_restrictions, OR_param2_restrictions, fix_corners, aggregation_flag, metric_flag, max_time_estimation, min_counter, parallel_flag, output_file_config, write_all_flag, bootstrap_flag, residual_flag, in_sample_flag, out_sample_flag] = HPZ_Settings_Read(main_folder)

% this function reads the user-interface screens saved settings

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% read settings from file 
try
    % reading from the file
    [settings] = csvread(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.settings_file_name, '.csv'));
    
    % in the next line, we call the 23th row and 8th column in the settings
    % matrix, because we have 23 rows and the max number of columns is 7
    % (it is 7 for GARP flags).
    % if it doesn't exists, it means that the settings file is corrupted,
    % cause it doesn't have all required data.
    % if in the future you will add more rows (>23) or use longer rows (num
    % of columns > 7), you should increase these numbers respectively.
    check_settings = settings(23 , 7); %#ok<NASGU>
catch
    % A) if we failed to read the file,
    % we assume it might be because it was accidently deleted, or because
    % it is not formatted in a way that allows it to be read.
    % B) if the file's data was corrupted, we assume someone might've
    % played with it and made changes to it by mistake.
    % A+B) either way, we reset the file, then read it again.
    HPZ_Settings_Reset(main_folder);
    [settings] = csvread(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.settings_file_name, '.csv'));
end

% assign initial values to variables in accordance with the settings
save_user_choices = settings(1,1);
fix_endowments = settings(2,1);
action_choice = settings(3,1); %
GARP_flags = settings(4,1:7);
AFRIAT_flags = settings(5,1:4);
VARIAN_flags = settings(6,1:4);
HOUTMAN_flags = settings(7,1:4);
risk_numeric_flag = settings(8,1);
OR_numeric_flag = settings(8,2);
risk_function_flag = settings(9,1); %
OR_function_flag = settings(9,2); %
risk_param1_restrictions = settings(10,1:2);
risk_param2_restrictions = settings(11,1:2);
if (risk_param1_restrictions(1) == -HPZ_Constants.infinity) 
    risk_param1_restrictions(1) = -inf; end
if (risk_param1_restrictions(2) == HPZ_Constants.infinity) 
    risk_param1_restrictions(2) = inf; end
if (risk_param2_restrictions(1) == -HPZ_Constants.infinity) 
    risk_param2_restrictions(1) = -inf; end
if (risk_param2_restrictions(2) == HPZ_Constants.infinity) 
    risk_param2_restrictions(2) = inf; end
OR_param1_restrictions = settings(12,1:2);
OR_param2_restrictions = settings(13,1:2);
if (OR_param1_restrictions(1) == -HPZ_Constants.infinity) 
    OR_param1_restrictions(1) = -inf; end
if (OR_param1_restrictions(2) == HPZ_Constants.infinity) 
    OR_param1_restrictions(2) = inf; end
if (OR_param2_restrictions(1) == -HPZ_Constants.infinity) 
    OR_param2_restrictions(1) = -inf; end
if (OR_param2_restrictions(2) == HPZ_Constants.infinity) 
    OR_param2_restrictions(2) = inf; end
fix_corners = settings(14,1);
aggregation_flag = settings(15,1); %
metric_flag = settings(16,1); %
max_time_estimation = settings(17,1);
if max_time_estimation <= 0
    max_time_estimation = inf;
end
min_counter = settings(18,1); %
if min_counter == HPZ_Constants.infinity
    min_counter = inf;
end
parallel_flag = settings(19,1);
output_file_config = settings(20,1:6);
write_all_flag = settings(21,1);
bootstrap_flag = settings(22,1);
residual_flag = settings(23,1);
in_sample_flag = settings(23,2);
out_sample_flag = settings(23,3);



%% for (almost) every variable, we check if it has a valid value, otherwise we reset it
% whether to save the choices (whether to use this function)
if (save_user_choices ~= 0 && save_user_choices ~= 1)
    save_user_choices = 1;
end
% fix endowments
if (fix_endowments ~= 0 && fix_endowments ~= 1)
    fix_endowments = 1;
end
% action flag
if (action_choice ~= HPZ_Constants.Consistency_action && action_choice ~= HPZ_Constants.NLLS_action && action_choice ~= HPZ_Constants.MMI_action && action_choice ~= HPZ_Constants.BI_action)
    action_choice = HPZ_Constants.Consistency_action;
end
% consistency flags
for i=1:length(GARP_flags)
    if (GARP_flags(i) ~= 0 && GARP_flags(i) ~= 1)
        GARP_flags(i) = 1;
    end
end
for i=1:length(AFRIAT_flags)
    if (AFRIAT_flags(i) ~= 0 && AFRIAT_flags(i) ~= 1)
        AFRIAT_flags(i) = 1;
    end
end
for i=1:length(VARIAN_flags)
    if (VARIAN_flags(i) ~= 0 && VARIAN_flags(i) ~= 1)
        VARIAN_flags(i) = 1;
    end
end
for i=1:length(HOUTMAN_flags)
    if (HOUTMAN_flags(i) ~= 0 && HOUTMAN_flags(i) ~= 1)
        HOUTMAN_flags(i) = 1;
    end
end
% numeric flags
if (risk_numeric_flag ~= 0 && risk_numeric_flag ~= 1)
    risk_numeric_flag = 1;
end
if (OR_numeric_flag ~= 0 && OR_numeric_flag ~= 1)
    OR_numeric_flag = 1;
end
% function flags
if (risk_function_flag ~= HPZ_Constants.CRRA_func && risk_function_flag ~= HPZ_Constants.CARA_func)
   risk_function_flag = HPZ_Constants.CRRA_func;
end
if (OR_function_flag ~= HPZ_Constants.CES_func)
   OR_function_flag = HPZ_Constants.CES_func;
end
% fix corners
if (fix_corners ~= 0 && fix_corners ~= 1)
    fix_corners = 1;
end
% aggregation and metric flags
if (aggregation_flag ~= HPZ_Constants.MMI_Max && aggregation_flag ~= HPZ_Constants.MMI_Mean && aggregation_flag ~= HPZ_Constants.MMI_AVGSSQ)
   aggregation_flag = HPZ_Constants.MMI_AVGSSQ;
end
if (metric_flag ~= HPZ_Constants.euclidean_metric && metric_flag ~= HPZ_Constants.CFGK_metric)
   metric_flag = HPZ_Constants.euclidean_metric;
end
% we check if min_counter is one of those that appear in the list in
% HPZ_Constants, otherwise we reset it
min_counter_values = HPZ_Constants.min_counter_values;
num_of_values = length(min_counter_values);
is_valid = false;
for i=1:num_of_values
    if min_counter == str2num(min_counter_values{i}) %#ok<*ST2NM>
        is_valid = true;
    end
end
if ~is_valid
    min_counter = str2num(min_counter_values{num_of_values});
end
% parallel flag
if (parallel_flag ~= 0 && parallel_flag ~= 1)
    parallel_flag = 1;
end
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