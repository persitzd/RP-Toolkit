function HPZ_Consistency_Indices_Settings_Reset(main_folder)

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
settings_file = fopen(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.consistency_indices_settings_file_name, '.csv'), 'w');

% close file when this function end
cleanup_close_file = onCleanup(@() fclose(settings_file));



%% printing the default settings, which are:
% Revealed Preference Graphs (.png) - by default we print none, and printing without edge weights 
fprintf(settings_file, '%s,%s,%s,%s,%s\n', '0', '0', '0', '0', '0');   % line (1)
% Power Test - by default don't perform power test, and if perform - with 1000 simulations,   
% and if it is risk preferences - standard randomization (not given FOXD)  
fprintf(settings_file, '%s,%s,%s\n', '0', '1000', '0');   % line (2)
% GARP - perform GARP WARP and SARP, without residuals of any sort,
% but default residuals is both in-sample and out-of-sample, and for GARP only 
fprintf(settings_file, '%s,%s,%s,%s,%s,%s,%s\n', '1', '0', '1', '1', '0', '1', '0');   % line (3)
% AFRIAT - perform AFRIAT, without residuals of any sort,
% but default residuals is only out-of-sample (in-sample is not optional in fact) 
fprintf(settings_file, '%s,%s,%s,%s\n', '1', '0', '0', '1');   % line (4)
% VARIAN - perform VARIAN, without residuals of any sort,
% but default residuals is both in-sample and out-of-sample,
% and by default calculate only mean and AVGSSQ aggregators, and not max aggregator 
fprintf(settings_file, '%s,%s,%s,%s,%s,%s,%s\n', '1', '0', '1', '1', '1','1','0');   % line (5)
% HOUTMAN-MAKS - perform HOUTMAN-MAKS, without residuals of any sort,
% but default residuals is only out-of-sample (in-sample is not optional in fact)
fprintf(settings_file, '%s,%s,%s,%s\n', '1', '0', '0', '1');   % line (6)
% MPI - perform MPI, without residuals of any sort,
% but default residuals is only out-of-sample (in-sample is not optional in fact)
fprintf(settings_file, '%s,%s,%s,%s\n', '1', '0', '0', '1');   % line (7)



end

