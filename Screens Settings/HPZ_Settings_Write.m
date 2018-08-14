function HPZ_Settings_Write(save_user_choices, fix_endowments, action_choice, GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags, risk_numeric_flag, OR_numeric_flag, risk_function_flag, OR_function_flag, risk_param1_restrictions, risk_param2_restrictions, OR_param1_restrictions, OR_param2_restrictions, fix_corners, aggregation_flag, metric_flag, max_time_estimation, min_counter, parallel_flag, output_file_config, write_all_flag, bootstrap_flag, residual_flag, in_sample_flag, out_sample_flag, main_folder)

% This function is called when the user wishes (and it is also the
% default...) to save the choices he made now for the next time.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% preparations for saving the current decisions and settings of the user 
if max_time_estimation == inf
    max_time_estimation = 0;
end
if min_counter == inf
    min_counter = HPZ_Constants.infinity;
end
if (risk_param1_restrictions(1) == -inf) 
    risk_param1_restrictions(1) = -HPZ_Constants.infinity; end
if (risk_param1_restrictions(2) == inf) 
    risk_param1_restrictions(2) = HPZ_Constants.infinity; end
if (risk_param2_restrictions(1) == -inf) 
    risk_param2_restrictions(1) = -HPZ_Constants.infinity; end
if (risk_param2_restrictions(2) == inf) 
    risk_param2_restrictions(2) = HPZ_Constants.infinity; end
if (OR_param1_restrictions(1) == -inf) 
    OR_param1_restrictions(1) = -HPZ_Constants.infinity; end
if (OR_param1_restrictions(2) == inf) 
    OR_param1_restrictions(2) = HPZ_Constants.infinity; end
if (OR_param2_restrictions(1) == -inf) 
    OR_param2_restrictions(1) = -HPZ_Constants.infinity; end
if (OR_param2_restrictions(2) == inf) 
    OR_param2_restrictions(2) = HPZ_Constants.infinity; end



% make sure the Settings Files directory exists, if not - create it
dir_exists = exist(strcat(main_folder, '/', HPZ_Constants.settings_files_dir) , 'dir');
if ~dir_exists
    mkdir(HPZ_Constants.settings_files_dir);
end

% creating/opening the settings file
settings_file = fopen(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.settings_file_name, '.csv'), 'w');

% close file when this function end
cleanup_close_file = onCleanup(@() fclose(settings_file));



%% printing the new settings and choices

% save user choices for next time
fprintf(settings_file, '%s\n', num2str(save_user_choices));   % line (1)
% fix_endowments
fprintf(settings_file, '%s\n', num2str(fix_endowments));   % line (2)
% action
fprintf(settings_file, '%s\n', num2str(action_choice));   % line (3)
% GARP - perform GARP WARP and SARP, with/without residuals
for i=1:max(size(GARP_flags))
    fprintf(settings_file, '%s,', num2str(GARP_flags(i)));   % line (4)
end
fprintf(settings_file, '\n');
% AFRIAT - perform AFRIAT, with/without residuals
for i=1:max(size(AFRIAT_flags))
    fprintf(settings_file, '%s,', num2str(AFRIAT_flags(i)));   % line (5)
end
fprintf(settings_file, '\n');
% VARIAN - perform VARIAN, with/without residuals
for i=1:max(size(VARIAN_flags))
    fprintf(settings_file, '%s,', num2str(VARIAN_flags(i)));   % line (6)
end
fprintf(settings_file, '\n');
% HOUTMAN-MAKS - perform HOUTMAN-MAKS, with/without residuals
for i=1:max(size(HOUTMAN_flags))
    fprintf(settings_file, '%s,', num2str(HOUTMAN_flags(i)));   % line (7)
end
fprintf(settings_file, '\n');
% MPI - perform MPI, with/without residuals
for i=1:max(size(MPI_flags))
    fprintf(settings_file, '%s,', num2str(MPI_flags(i)));   % line (8)
end
fprintf(settings_file, '\n');
% numeric (risk and OR)
fprintf(settings_file, '%s,%s\n', num2str(risk_numeric_flag), num2str(OR_numeric_flag));   % line (8)
% functional form (risk and OR)
fprintf(settings_file, '%s,%s\n', num2str(risk_function_flag), num2str(OR_function_flag));   % line (9)
% restrictions on the 1st parameter (beta) in risk preferences
fprintf(settings_file, '%s,%s\n', num2str(risk_param1_restrictions(1)), num2str(risk_param1_restrictions(2)));   % line (10)
% restrictions on the 2nd parameter (rho/A) in risk preferences
fprintf(settings_file, '%s,%s\n', num2str(risk_param2_restrictions(1)), num2str(risk_param2_restrictions(2)));   % line (11)
% restrictions on the 1st parameter (alpha) in OR preferences
fprintf(settings_file, '%s,%s\n', num2str(OR_param1_restrictions(1)), num2str(OR_param1_restrictions(2)));   % line (12)
% restrictions on the 2nd parameter (rho) in OR preferences
fprintf(settings_file, '%s,%s\n', num2str(OR_param2_restrictions(1)), num2str(OR_param2_restrictions(2)));   % line (13)
% fix corners (CFGK) or not
fprintf(settings_file, '%s\n', num2str(fix_corners));   % line (14)

% aggregation_flag for MMI (Max / Mean / AVGSSQ)
fprintf(settings_file, '%s\n', num2str(aggregation_flag));   % line (15)
% metric_flag for NLLS (euclidean or CFGK)
fprintf(settings_file, '%s\n', num2str(metric_flag));   % line (16)
% max_time_estimation in minutes ('0' if there is no limit)
if isempty(max_time_estimation)
    max_time_estimation = 0;
end
fprintf(settings_file, '%s\n', num2str(max_time_estimation));   % line (17)
% min_counter
fprintf(settings_file, '%s\n', num2str(min_counter));   % line (18)
% parallel_flag
fprintf(settings_file, '%s\n', num2str(parallel_flag));   % line (19)

% output_file_config (NLLS-euclidean , NLLS-CFGK , NLLS-normalized-euclidean , MMI-Max , MMI-Mean , MMI-AVGSSQ, BI)
for i=1:max(size(output_file_config))
    fprintf(settings_file, '%s,', num2str(output_file_config(i)));   % line (20)
end
fprintf(settings_file, '\n');
% write_all_flag
fprintf(settings_file, '%s\n', num2str(write_all_flag));   % line (21)
% bootstrap_flag
fprintf(settings_file, '%s\n', num2str(bootstrap_flag));   % line (22)

% residual , in-sample , out-of-sample (only for parameters estimation)  
% (default is without residuals, but when performing residuals, do both in-sample and out-of-sample) 
fprintf(settings_file, '%s,%s,%s\n', num2str(residual_flag), num2str(in_sample_flag), num2str(out_sample_flag));   % line (23)



end