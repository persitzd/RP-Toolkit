function HPZ_Debug_Interface

% this function is used to perform calculations for the purpose of
% debugging, and possibly for other purposes as well
% it allows to easily calculate the various criteria for subjects, assuming
% any desired parameters, and without the need to perform a full estimation



%% general preparations

% before we start meddling with the warnings ('on' and 'off'), we first
% save the warning settings, and we call onCleanup in order to restore
% these settings when the function finishes or stops for any other reason
original_state = warning;
cleanup_warnings = onCleanup(@() warning(original_state));

% adding all the subfolders to the code path:
% supress warning. the warning is due to function 'assert' in the 'test'
% folder of 'matlab_bgl', that has the same name as MATLAB builtin
warning('off','MATLAB:dispatcher:nameConflict');
% determine where your m-file's folder is (we need it in other functions as well) 
debug_folder = fileparts(which(mfilename)); 
mydir = debug_folder; % mydir  = pwd;
idcs   = strfind(mydir,'/');
if isempty(idcs)
    idcs   = strfind(mydir,'\');
end
main_folder = mydir(1:(idcs(end)-1)) %newdir = mydir(1:idcs(end)-1); 
% add that folder plus all subfolders to the path
addpath(genpath(main_folder));
% as default, the results files will be in a "Results" folder under the
% main folder which contains the code. but if you wish to have a different
% folder for the results, you can change the following line as needed
% ( e.g. - main_folder_for_results = pwd; )
main_folder_for_results = main_folder;



%% Select a data set

% getting the saved fix_endowment settings
[~, fix_endowments, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = HPZ_Settings_Read(main_folder);
% Also select whether a fix of the endowments in the data to exactly 1 is
% required or not.
[data_matrix, pref_class, ~, ok] = HPZ_Screen_Data_Set_Selection(fix_endowments, main_folder, 0);
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



%% if the chosen data set is of risk preferences: ask the user what functional form is desired

if pref_class == HPZ_Constants.risk_pref
    
    function_flag = 0;

    % dimensions of the list window
    list_size = [500 , 120];

    while ~((function_flag == HPZ_Constants.CRRA_func) || (function_flag == HPZ_Constants.CARA_func))

        [function_flag, ok] = listdlg('PromptString','Function Form Selection', 'InitialValue',1,...
            'SelectionMode','single', 'ListString',{'CRRA','CARA'},...
            'Name','Action Selection', 'ListSize',list_size, 'uh',30, 'fus',8, 'ffs',8);

        % if the user pressed "cancel" - end the program
        if (ok == 0)
            return
        end

    end
    
elseif pref_class == HPZ_Constants.OR_pref
    
    function_flag = HPZ_Constants.CES_func;
    
end



%% select calculation solution approach (analytic, semi-numeric or numeric)
numeric_flag = -1;

% dimensions of the list window
list_size = [500 , 120];

while ~((numeric_flag == HPZ_Constants.analytic) || numeric_flag == HPZ_Constants.semi_numeric || (numeric_flag == HPZ_Constants.numeric))

    [numeric_flag, ok] = listdlg('PromptString','Solution Approach Selection', 'InitialValue',1,...
        'SelectionMode','single', 'ListString',{'analytic','numeric','semi-numeric'},...
        'Name','Action Selection', 'ListSize',list_size, 'uh',30, 'fus',8, 'ffs',8);
    numeric_flag = numeric_flag - 1;
    
    % if the user pressed "cancel" - end the program
    if (ok == 0)
        return
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



%% Select a file that contains the subjects and their parameters
% in this file, the columns should be EXACTLY as follows:
%   #1 : subject ID (number)
%   #2 : 1st parameter
%   #3 : 2nd parameter
% browse for a data file
[FileName, PathName, ok] = uigetfile('../*.csv', 'Select a CSV Data File', debug_folder);
% if the user chose a file (didn't cancel)
if (ok ~= 0)
    % reading the data file
    [data] = csvread(strcat(PathName,FileName));
    % formatting the matrix as desired
    param_mat = [data(:,1) , data(:,2) , data(:,3)];
else
    return
end



%% Which of the Criteria to calculate
which_criteria = HPZ_Debug_Criteria_To_Calculate();









%% Now we can finally start the calculations

% how many sets of parameters we have (at most)
[num_of_operations, ~] = size(param_mat);

% the chosen subjects
subjects = unique(data_matrix(:,1), 'stable');  % array of subjects' IDs
chosen_subjects = subjects(subjects_index);
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



% prepare results file/s, including column headers
[file_handle, ~, ~, ~] = HPZ_Write_Result_File_Initializer(which_criteria, HPZ_Constants.NLLS_action, pref_class, function_flag, {'Criterion','Criterion','Criterion','Criterion','Criterion','Criterion','Criterion'}, 0, 0, main_folder_for_results);

% make sure the file closes when execution stops for any reason 
cleanup_close_file = onCleanup(@() fclose(file_handle));



% loop over these sets of parameters
for i=1:num_of_operations
    
    current_chosen_subject = find(chosen_subjects == param_mat(i,1), 1);
    
    % if this is one of the chosen subjects, only then we calculate
    if ~isempty(current_chosen_subject)
        
        current_chosen_subject_index = (find(subjects == param_mat(i,1), 1));
        
        % extracting the data that is necessary for all methods and calculations
        if current_chosen_subject_index < num_subjects
            data = data_matrix(first_obs_row(subjects_index(current_chosen_subject_index)):((first_obs_row(subjects_index(current_chosen_subject_index)+1))-1),:);
            %obs_num = data_matrix((first_obs_row(subjects_index(current_chosen_subject_index)+1))-1,2);
        else
            data = data_matrix(first_obs_row(subjects_index(current_chosen_subject_index)):rows,:);
            %obs_num = data_matrix(rows,2);
        end

        % ...
        asymmetric_flag = 1;    % As is
        treatment = 1;          % As is
        fix_corners = false;    % we don't bother to support this option
        
        % the current parameters
        param = param_mat(i,2:3);
        
        % extracting the observations
        observations = data(:,3:6);
        % calculating the endowments (assumed to be 1)
        endowments = sum(observations(:,1:2) .* observations(:,3:4), 2);
        
        
        
        % calculations of criteria takes place here
        final_output = zeros(1,9);
        final_output(1:2) = param;
        % NLLS Criterion (Euclidean metric, Choi et al. (2007) metric)
        [final_output(3), final_output(4), final_output(5), ~] = HPZ_NLLS_Metrics (param, endowments, observations, treatment, function_flag, fix_corners, asymmetric_flag, pref_class, numeric_flag);
        
        % MMI Criterion (Max Waste, Mean Waste, Sum of Squares Wastes)
        [final_output(6), final_output(7), final_output(8), ~] = HPZ_MMI_Aggregates(param, endowments, observations, treatment, function_flag, pref_class, numeric_flag);
        
        % BI Criterion
        [final_output(9), ~] = HPZ_BI_Criterion(param, endowments, observations, treatment, function_flag, pref_class, numeric_flag);
        
        
        
        % writing the results to the results file
        HPZ_Write_Result_File_Finalizer(file_handle, which_criteria, final_output, ...
                                        1, num2str(param_mat(i,1)), 0, HPZ_Constants.print_precision);
    end
    
end



% Address the folder where the result files were saved
msgbox(char(strcat('All files are saved under the following path:', {' '}, main_folder_for_results, '\', HPZ_Constants.results_files_dir)), 'Output File Path','help');
% Opening the folder (works only in Windows)
if strcmp(computer('arch'),'win64') || strcmp(computer('arch'),'win32')
    winopen(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir));
end



end