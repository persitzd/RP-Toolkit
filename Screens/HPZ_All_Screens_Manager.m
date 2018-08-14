function [ok, data_matrix, subjects_index, action_flag, GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags, pref_class, function_flag, numeric_flag, param1_restrictions, param2_restrictions, fix_corners, metric_flag, aggregation_flag, max_time_estimation, max_starting_points, min_counter, parallel_flag, output_file_config, write_all_flag, bootstrap_flag, file_val_str, residual_flag, in_sample_flag, out_sample_flag] = HPZ_All_Screens_Manager(main_folder, runs_counter)

% This function calls to all relevant user-interface screens one after
% another, and returns all the choices and user-preferences received from all of them 

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% getting the saved settings
[save_user_choices, fix_endowments, action_choice, GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags, risk_numeric_flag, OR_numeric_flag, risk_function_flag, OR_function_flag, risk_param1_restrictions, risk_param2_restrictions, OR_param1_restrictions, OR_param2_restrictions, fix_corners, aggregation_flag, metric_flag, max_time_estimation, min_counter, parallel_flag, output_file_config, write_all_flag, bootstrap_flag, residual_flag, in_sample_flag, out_sample_flag] = HPZ_Settings_Read(main_folder);

% initialization of all returned variables
% (most of them are read from the settings file (such as function_flag that
% is read but there are risk_function_flag and OR_function_flag),
% the rest are initialized just to avoid errors)  
subjects_index = 0;
action_flag = 0;
function_flag = 0;
numeric_flag = 0;
param1_restrictions = [-inf , inf];
param2_restrictions = [-inf , inf];
max_starting_points = 1000;
fix_corners_save = fix_corners;
output_file_config_save = output_file_config;
file_val_str = {'Value', 'Value', 'Value', 'Value', 'Value', 'Value'};





%% Select a data set.
% Also select whether a fix of the endowments in the data to exactly 1 is
% required or not.

[data_matrix, pref_class, fix_endowments, ok] = HPZ_Screen_Data_Set_Selection(fix_endowments, main_folder, runs_counter);
% if the user pressed "cancel" - end the program
if (ok == 0)
    return
end





% list of subjects' numbers in the required treatment

% 'stable' is essential to keep the original order of subjects as it is
% in the file data
subjects = unique(data_matrix(:,1), 'stable');  % array of subjects' IDs
subjects_str = num2str(subjects);               % as strings
% number of subjects
num_subjects = length(subjects);
% an array of location indices - the i'th element 		
% is the index of the first observation of the subject i 
[rows,~] = size(data_matrix);
first_obs_row = zeros(num_subjects, 1);
counter_subjects = 1;
for i=1:rows
    if data_matrix(i,2) == 1
        first_obs_row (counter_subjects) = i;
        counter_subjects = counter_subjects + 1;
    end
end





%% Select a calculation to be performed on the data

% dimensions of the list window
pc_screen_size = get(0,'ScreenSize');
pc_screen_height = pc_screen_size(4);
list_height =  min(max(100 , 13.7*max(size(HPZ_Constants.all_actions_names))) , HPZ_Constants.max_height_percent * pc_screen_height - HPZ_Constants.listdlg_extra_height); 
list_size = [500 , list_height];

while isempty(action_flag) || ~((action_flag == HPZ_Constants.Consistency_action) || (action_flag == HPZ_Constants.NLLS_action) || ...
                                    (action_flag == HPZ_Constants.MMI_action) || (action_flag == HPZ_Constants.BI_action))
    
    [action_choice, ok] = listdlg('PromptString','Action Selection', 'InitialValue',action_choice,...
        'SelectionMode','single', 'ListString',HPZ_Constants.all_actions_names,...
        'Name','Action Selection', 'ListSize',list_size, 'uh',30, 'fus',8, 'ffs',8);
    
    % if the user pressed "cancel" - end the program
    if (ok == 0)
        return
    end
    
    % it may seem tedious, but this approach will prevent the whole program 
    % from collapsing/misfunctioning if someone changes the values of some
    % of the constants
    if strcmp(HPZ_Constants.all_actions_names{action_choice} , HPZ_Constants.Consistency_action_name)
        action_flag = HPZ_Constants.Consistency_action;   % Consistency indices
    elseif strcmp(HPZ_Constants.all_actions_names{action_choice} , HPZ_Constants.NLLS_action_name)
        action_flag = HPZ_Constants.NLLS_action;   % NLLS estimation
    elseif strcmp(HPZ_Constants.all_actions_names{action_choice} , HPZ_Constants.MMI_action_name)
        action_flag = HPZ_Constants.MMI_action;   % MMI estimation
    elseif strcmp(HPZ_Constants.all_actions_names{action_choice} , HPZ_Constants.BI_action_name)
        action_flag = HPZ_Constants.BI_action;   % BI estimation
    end
end





%% Select the subjects from the data that the calculation will be performed on

% dimensions of the list window
pc_screen_size = get(0,'ScreenSize');
pc_screen_height = pc_screen_size(4);
list_height = min(max(100 , 13.7*max(size(subjects_str))) , HPZ_Constants.max_height_percent * pc_screen_height - 2*HPZ_Constants.listdlg_extra_height); 
list_size = [400 , list_height];

legal_assignment = 0;

while ~(legal_assignment == 1)
    
    legal_assignment = 1;

    [subjects_index, ok] = listdlg('PromptString','Subjects Selection (the same action applies to all)',...
                'SelectionMode','multiple', 'ListString',subjects_str,...
                'Name','Subjects', 'ListSize',list_size, 'uh',30, 'fus',8, 'ffs',8);
    
    % if the user pressed "cancel" - end the program
    if (ok == 0)
        return
    end     
    
    chosen_subjects_num = length(subjects_index);
    
    % if the user pressed "ok" without choosing any subject,
    % don't let him escape the loop, show him this listdlg again
    if (chosen_subjects_num == 0)
        legal_assignment = 0;
    end
    
    for i=1:chosen_subjects_num

        % This code is in order to prevent a bug that may occur if a 
        % subject that does not exist is chosen. This should not ever 
        % occur, but in case it occurs the bug will be stopped at this stage. 
        if isempty(find(subjects == data_matrix(first_obs_row(subjects_index(i)),1), 1))
            
            legal_assignment = 0;
            
        end
    end
end





%% consistency and inconsistency measures settings
if action_flag == HPZ_Constants.Consistency_action   % Consistency indices
    
    % here the user can choose which consistency and inconsistency measures
    % to calculate, and whether to calculate residuals for them, and which
    % method of residuals (in sample or out of sample)
    [GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags, ok] = HPZ_Screen_Consistency_Indices(GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags, runs_counter);  
    
    % if the user clicked "cancel", stop the program
    if (ok == 0)
        return
    end
end





%% if estimation was chosen we need some more information
if ((action_flag == HPZ_Constants.NLLS_action) || (action_flag == HPZ_Constants.MMI_action) || (action_flag == HPZ_Constants.BI_action))    
    
    if pref_class == HPZ_Constants.risk_pref
        % the following is specific information needed for the estimation  
        % of risk preferences
        [risk_numeric_flag, risk_function_flag, risk_param1_restrictions, risk_param2_restrictions, fix_corners, ok] = HPZ_Screen_Functional_Form_Settings_Risk(action_flag, risk_numeric_flag, risk_function_flag, risk_param1_restrictions, risk_param2_restrictions, fix_corners, runs_counter);
        % fix_corners_save is the value that will be saved in settings
        fix_corners_save = fix_corners;
        % if the user clicked "cancel", stop the program
        if (ok == 0)
            return
        end
        param1_restrictions = risk_param1_restrictions;
        param2_restrictions = risk_param2_restrictions;
        function_flag = risk_function_flag;
        numeric_flag = risk_numeric_flag;
    elseif pref_class == HPZ_Constants.OR_pref
        % the following is specific information needed for the estimation 
        % of Kurtz et al. (2016) data set.
        fix_corners = false;   % allow corner choices
        [OR_numeric_flag, OR_function_flag, OR_param1_restrictions, OR_param2_restrictions, ok] = HPZ_Screen_Functional_Form_Settings_OR(action_flag, OR_numeric_flag, OR_function_flag, OR_param1_restrictions, OR_param2_restrictions, runs_counter);        
        % if the user clicked "cancel", stop the program
        if (ok == 0)
            return
        end
        param1_restrictions = OR_param1_restrictions;
        param2_restrictions = OR_param2_restrictions;
        function_flag = OR_function_flag;
        numeric_flag = OR_numeric_flag;
    end
    
	% if estimation was chosen we need some more information
	[aggregation_flag, metric_flag, max_time_estimation, min_counter, parallel_flag, max_starting_points, ok] = HPZ_Screen_Optimization_Settings(aggregation_flag, metric_flag, max_time_estimation, min_counter, parallel_flag, action_flag, numeric_flag, runs_counter);
    % if the user clicked "cancel", stop the program
    if (ok == 0)
        return
    end
    
    % set corners:
    if metric_flag == HPZ_Constants.CFGK_metric   % CFGK metric (not Euclidean metric)
        % if we use CFGK metric, we must not allow corner choices, 
        % even if the user forgot to choose to not allow corner choices
        fix_corners = true;
    end
    
    [output_file_config, output_file_config_save, write_all_flag, bootstrap_flag, file_val_str, ok] = HPZ_Screen_Output_File_Format(output_file_config, write_all_flag, bootstrap_flag, action_flag, aggregation_flag, metric_flag, runs_counter);
    % if the user clicked "cancel", stop the program
    if (ok == 0)
        return
    end
    
    %% residual settings
    [residual_flag, in_sample_flag, out_sample_flag, ok] = HPZ_Screen_Residual_Calculation_Settings(residual_flag, in_sample_flag, out_sample_flag, runs_counter);
    % if the user clicked "cancel", stop the program
    if (ok == 0)
        return
    end
end



if true%(save_user_choices)
    % saving the user's decisions and settings for next time
    HPZ_Settings_Write(save_user_choices, fix_endowments, action_choice, GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags, risk_numeric_flag, OR_numeric_flag, risk_function_flag, OR_function_flag, risk_param1_restrictions, risk_param2_restrictions, OR_param1_restrictions, OR_param2_restrictions, fix_corners_save, aggregation_flag, metric_flag, max_time_estimation, min_counter, parallel_flag, output_file_config_save, write_all_flag, bootstrap_flag, residual_flag, in_sample_flag, out_sample_flag, main_folder);
end



end