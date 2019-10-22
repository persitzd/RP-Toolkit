function HPZ_Print_Intlinprog_Internal_Error_Details(main_folder_for_results, f, intcon, A,b, Aeq,beq, lb,ub)

% when an internal error occurs during a call to intlinprog, 
% the program will print the details to a file

% (just to be sure it won't stop the next runs because of an error - we use here try/catch as well) 
try
    % make sure the Results directory exists, if not - create it
    dir_exists = exist(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir) , 'dir');
    if ~dir_exists
        mkdir(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir));
    end
    
    % name of results file
    file_name_str = char(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir, '/', 'Intlingprog Internal Error', {' , '}, version('-release'), {' , '}, date, '-' , generate_random_str));
    file_name = strcat(file_name_str, '.xlsx');
    % opening the results file
    %file_handle = fopen(strcat(file_name_str, '.csv'), 'w'); 
    
    % make sure the file closes when execution stops for any reason 
    %cleanup_close_file = onCleanup(@() fclose(file_handle));

    % printing the details to the file 
    %fprintf(file_handle, '%s,%s\n', error_var.identifier, error_var.message);
    try
        to_print_cell = {f, intcon, A,b, Aeq,beq, lb,ub};
        for i = 1:length(to_print_cell)
            element_to_print = to_print_cell{i};
            if ~isempty(element_to_print)
                xlswrite(file_name, element_to_print, i);
            end
        end
    catch error_var
       error_identifier = error_var.identifier
       error_message = error_var.message 
    end
    
catch error_print %#ok<NASGU>
    
    msgbox('An error occured while trying to print the details of the Intlinprog Internal Error.', 'An Internal Error Occured', 'help');
    
end

end