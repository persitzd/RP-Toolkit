function [bootstrap_sample_sizes, bootstrap_significance_level, BI_threshold, max_starting_points, possible_num_convergence_points, one_residuals_file, debugger_mode, print_single_subject, waitbar_settings, Varian_algorithm_settings] = HPZ_Advanced_Options_Settings_Read(main_folder)

% this function reads the user-interface screens saved settings

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% read settings from file 
try
    % reading settings (except varian aggregators) from the file
    [settings] = csvread(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.advanced_options_settings_file_name, '.csv'));
    
    % in the next line, we call the 9th row and 3rd column in the settings
    % matrix, because we have 10 rows and the max number of columns is 3.
    % if it doesn't exists, it means that the settings file is corrupted,
    % cause it doesn't have all required data.
    % if in the future you will add more rows (>10) or use longer rows (num
    % of columns > 3), you should increase these numbers respectively.
    check_settings = settings(10 , 3); %#ok<NASGU>
    
catch
    % A) if we failed to read the file,
    % we assume it might be because it was accidently deleted, or because
    % it is not formatted in a way that allows it to be read.
    % B) if the file's data was corrupted, we assume someone might've
    % played with it and made changes to it by mistake.
    % A+B) either way, we reset the file, then read it again.
    
    HPZ_Advanced_Options_Settings_Reset(main_folder);
    
    % reading settings (except varian aggregators) from the file
    [settings] = csvread(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.advanced_options_settings_file_name, '.csv'));
    
end



% assign initial values to variables in accordance with the settings
bootstrap_sample_sizes = settings(1,1:3);
bootstrap_significance_level = settings(2,1);
BI_threshold = settings(3,1);
max_starting_points = settings(4,1:3);
possible_num_convergence_points = settings(5,:);
one_residuals_file = settings(6,1);
debugger_mode = settings(7,1);
waitbar_settings = settings(8,1:4);
Varian_algorithm_settings = settings(9,1:2);
print_single_subject = settings(10,1);


% for (almost) every variable, we check if it has a valid value, otherwise we reset it

% bootstrap_sample_sizes
for i=1:length(bootstrap_sample_sizes)
    if bootstrap_sample_sizes(i) < 10
        % we don't allow bootstrap with less than 10 iterations
        bootstrap_sample_sizes(i) = 10;
    elseif bootstrap_sample_sizes(i) > 1000000
        % we don't allow bootstrap with more than 1,000,000 iterations
        bootstrap_sample_sizes(i) = 1000000;
    end
end

% we keep bootstrap_significance_level in the reasonable range
if bootstrap_significance_level < 0
    bootstrap_significance_level = 0;
elseif bootstrap_significance_level > 0.5
    bootstrap_significance_level = 0.5;
end

% we keep BI_threshold in the reasonable range
if BI_threshold < 0
    BI_threshold = 0;
elseif BI_threshold > 1
    BI_threshold = 1;
end

% max_starting_points
for i=1:length(max_starting_points)
    if max_starting_points(i) < 1
        max_starting_points(i) = 1;
    end
end
% possible_num_convergence_points
possible_num_convergence_points = possible_num_convergence_points(mod(possible_num_convergence_points,1) == 0 & possible_num_convergence_points > 0);
possible_num_convergence_points = sort(possible_num_convergence_points);
if isempty(possible_num_convergence_points)
    possible_num_convergence_points = [3,5,10,15,20,30];
end

% variables that are either 0 or 1 - we make sure they are that way
if one_residuals_file ~= 0 && one_residuals_file ~= 1
    one_residuals_file = 1;
end
if debugger_mode ~= 0 && debugger_mode ~= 1
    debugger_mode = 1;
end
if print_single_subject ~= 0 && print_single_subject ~= 1
    print_single_subject = 1;
end
if waitbar_settings(2) ~= 0 && waitbar_settings(2) ~= 1
    waitbar_settings(2) = 1;
end
if waitbar_settings(3) ~= 0 && waitbar_settings(3) ~= 1
    waitbar_settings(3) = 1;
end
if waitbar_settings(4) ~= 0 && waitbar_settings(4) ~= 1
    waitbar_settings(4) = 1;
end

% minimum number of subjects so that we use a single waitbar for all of them 
if waitbar_settings(1) < 1
    waitbar_settings(1) = 1;
end

% Varian_algorithm_settings
if Varian_algorithm_settings(1) < 2
    Varian_algorithm_settings(1) = 2;
end
if Varian_algorithm_settings(2) < 1
    % it cannot be less than 1; can result an infinite loop in Varian
    Varian_algorithm_settings(2) = 1;
end



end