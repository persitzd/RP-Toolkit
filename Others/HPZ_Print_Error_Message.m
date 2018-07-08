function HPZ_Print_Error_Message(error_var, main_folder_for_results, current_run, total_runs)

% when an error occurs during a run, the program will print the error
% message to a file

% error_var - the exception that was thrown by the error

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% (just to be sure it won't stop the next runs because of an error - we use here try/catch as well) 
try
    % make sure the Results directory exists, if not - create it
    dir_exists = exist(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir) , 'dir');
    if ~dir_exists
        mkdir(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir));
    end
    
    % name of results file
    file_name_str = char(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir, '/', 'Error in Run No.', {' '}, num2str(current_run), ' Out of', {' '}, num2str(total_runs), {' , '}, date, '-' , generate_random_str));
    % opening the results file
    file_handle = fopen(strcat(file_name_str, '.csv'), 'w'); 
    
    % make sure the file closes when execution stops for any reason 
    cleanup_close_file = onCleanup(@() fclose(file_handle));

    % printing the error messages and details to the file 
    fprintf(file_handle, '%s,%s\n', error_var.identifier, error_var.message);
    stack_size = size(error_var.stack);
    for i=1:stack_size
        fprintf(file_handle, '%s,%s,%s\n', error_var.stack(i).file, error_var.stack(i).name, char(strcat('Line', {' '}, num2str(error_var.stack(i).line))));
    end
    
catch error_print
    
    msgbox(char(strcat('An error occured in Run No.', {' '}, num2str(current_run), {' '}, ' Out of', {' '}, num2str(total_runs), {', '}, 'and due to another error, printing the error details to a file have failed.')), 'An Error Occured', 'help');
    
    try
        fprintf('%s,%s\n', error_var.identifier, error_var.message);
        stack_size = size(error_var.stack);
        for i=1:stack_size
            fprintf('Error in the Run: File: %s, Name: %s, Line: %s\n', error_var.stack(i).file, error_var.stack(i).name, char(strcat('Line', {' '}, num2str(error_var.stack(i).line))));
        end
    catch
        error_var
    end
    
    try
        fprintf('%s,%s\n', error_print.identifier, error_print.message);
        stack_size = size(error_print.stack);
        for i=1:stack_size
            fprintf('Error in Printing the Error in the Run: File: %s, Name: %s, Line: %s\n', error_print.stack(i).file, error_print.stack(i).name, char(strcat('Line', {' '}, num2str(error_print.stack(i).line))));
        end
    catch
        error_print
    end
end

end