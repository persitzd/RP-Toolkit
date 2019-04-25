function [bootstrap_sample_sizes, bootstrap_significance_level, BI_threshold, max_starting_points, possible_num_convergence_points, one_residuals_file, debugger_mode, waitbar_settings, Varian_algorithm_settings] = HPZ_Advanced_Options_Settings_Default_Values()

% this function returns the default values of the advanced options

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 

% assign initial values to variables in accordance with the settings
bootstrap_sample_sizes = [1000 , 100 , 500];
bootstrap_significance_level = 0.05;
BI_threshold = 10^(-5);
max_starting_points = [100 , 30 , 100];
possible_num_convergence_points = [3,4,5,6,7,8,9,10,12,14,16,18,20,25,30];
one_residuals_file = 1;
debugger_mode = 0;
waitbar_settings = [3 , 1 , 0 , 1];
Varian_algorithm_settings = [1000 , 2];

end