function HPZ_Advanced_Options_Settings_Reset(main_folder)

% This function is used to restore the Settings file, in case it was
% lost/deleted for some reason, or if it was corrupted.
% This function is called when an error occurs while attempting to read Settings.csv, 
% but it can also run directly if you desire to reset all the settings.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 


% get the default values
[bootstrap_sample_sizes, bootstrap_significance_level, BI_threshold, max_starting_points, possible_num_convergence_points, one_residuals_file, debugger_mode, print_single_subject, waitbar_settings, Varian_algorithm_settings] = HPZ_Advanced_Options_Settings_Default_Values();

% write the default values to the settings file
HPZ_Advanced_Options_Settings_Write(bootstrap_sample_sizes, bootstrap_significance_level, BI_threshold, max_starting_points, possible_num_convergence_points, one_residuals_file, debugger_mode, print_single_subject, waitbar_settings, Varian_algorithm_settings, main_folder);


end

