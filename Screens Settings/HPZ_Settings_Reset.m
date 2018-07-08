function HPZ_Settings_Reset(main_folder)

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
settings_file = fopen(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.settings_file_name, '.csv'), 'w');

% close file when this function end
cleanup_close_file = onCleanup(@() fclose(settings_file));



%% printing the default settings, which are:

% save user choices for next time
fprintf(settings_file, '%s\n', '1');   % line (1)
% fix_endowments - '1' (fix endowments to be equal to exactly 1)
fprintf(settings_file, '%s\n', '1');   % line (2)

% action - '1' (the first in the list)
fprintf(settings_file, '%s\n', '1');   % line (3)

% GARP - perform GARP WARP and SARP, without residuals of any sort,
% but default residuals is both in-sample and out-of-sample, and for GARP only 
fprintf(settings_file, '%s,%s,%s,%s,%s,%s,%s\n', '1', '0', '1', '1', '0', '1', '0');   % line (4)
% AFRIAT - perform AFRIAT, without residuals of any sort,
% but default residuals is only out-of-sample (in-sample is not optional in fact) 
fprintf(settings_file, '%s,%s,%s,%s\n', '1', '0', '0', '1');   % line (5)
% VARIAN - perform VARIAN, without residuals of any sort,
% but default residuals is both in-sample and out-of-sample 
fprintf(settings_file, '%s,%s,%s,%s\n', '1', '0', '1', '1');   % line (6)
% HOUTMAN-MAKS - perform HOUTMAN-MAKS, without residuals of any sort,
% but default residuals is only out-of-sample (in-sample is not optional in fact)
fprintf(settings_file, '%s,%s,%s,%s\n', '1', '0', '0', '1');   % line (7)

% not numeric (analytic is default for both risk and OR)
fprintf(settings_file, '%s,%s\n', '0','0');   % line (8)
% functional form (set to first for both risk and OR)
fprintf(settings_file, '%s,%s\n', '1','1');   % line (9)
% restrictions on the 1st parameter (beta) in risk preferences
fprintf(settings_file, '%s,%s\n', '-1', num2str(HPZ_Constants.infinity));   % line (10)
% restrictions on the 2nd parameter (rho/A) in risk preferences
fprintf(settings_file, '%s,%s\n', '0', num2str(HPZ_Constants.infinity));   % line (11)
% restrictions on the 1st parameter (alpha) in OR preferences
fprintf(settings_file, '%s,%s\n', '0', '1');   % line (12)
% restrictions on the 2nd parameter (rho) in OR preferences
fprintf(settings_file, '%s,%s\n', num2str(-HPZ_Constants.infinity), num2str(HPZ_Constants.infinity));   % line (13)
% fix corners (CFGK) or not - default is not
fprintf(settings_file, '%s\n', '0');   % line (14)

% aggregation_flag for MMI (Max / Mean / AVGSSQ) default is AVGSSQ
fprintf(settings_file, '%s\n', '3');   % line (15)
% metric_flag for NLLS (euclidean or CFGK) default is euclidean
fprintf(settings_file, '%s\n', '1');   % line (16)
% max_time_estimation in minutes ('0' if there is no limit)
fprintf(settings_file, '%s\n', '0');   % line (17)
% min_counter default is the highest possible except for inf
fprintf(settings_file, '%s\n', num2str(cell2mat(HPZ_Constants.min_counter_values(end))));   % line (18)
% parallel_flag default is false
fprintf(settings_file, '%s\n', '0');   % line (19)

% output_file_config (NLLS-eaclidean , NLLS-CFGK , MMI-Max , MMI-Mean , MMI-AVGSSQ, BI)
% default is to print none of them, except for the one being estimated 
fprintf(settings_file, '%s,%s,%s,%s,%s,%s\n', '0', '0', '0', '0', '0', '0');   % line (20)
% write_all_flag (default is not to write all, but only the best one)
fprintf(settings_file, '%s\n', '0');   % line (21)
% bootstrap_flag (default is without bootstrap)
fprintf(settings_file, '%s\n', '0');   % line (22)

% residual , in-sample , out-of-sample (only for parameters estimation)  
% (default is without residuals, but when performing residuals, do both in-sample and out-of-sample) 
fprintf(settings_file, '%s,%s,%s\n', '0', '1', '1');   % line (23)



end

