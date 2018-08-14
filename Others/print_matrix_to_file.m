function print_matrix_to_file (file_name, matrix, precision_string, col_headers, row_headers)

% this function prints a specified 2D matrix with specified column/row headers  
% to a new CSV file

% file_name - name for the new CSV file (not including the '.csv' in the end) 
% matrix - the data matrix (2 dimensions) to be printed
% precision_string - a string in the format: %10.XXf or %10.XXg that determines how 
%           many digits after the point (f) or how many significant digits (g)
%           should be printed.
%           (if set to any non-cell argument, it prints as the matlab default dictates)
% col_headers - a cell of strings that will be used as column headers (optional) 
% row_headers - a cell of strings that will be used as row headers (optional) 
%       (if no column and/or no row headers are required, enter any non-cell argument to these parameters) 



% creating the new file
% (if the file already exists, it will overwrite it; if it exists and is
% currently open, this function will fail; please note that when using this
% function)
file_var = fopen(strcat(file_name, '.csv'), 'w');
% make sure the file closes when execution stops for any reason 
cleanup_close_file = onCleanup(@() fclose(file_var));



if iscell(col_headers)
    % if column headers are required - print them to the file
    
    % we first take one empty step if there are also row headers
    if iscell(row_headers)
        fprintf(file_var, ',');
    end
    % now we print the column headers
    num_of_col_headers = length(col_headers);
    for i=1:num_of_col_headers
        fprintf(file_var, '%s,', col_headers{i});
    end
    fprintf(file_var, '\n');
end



% size of matrix
[rows , cols] = size(matrix);



% number of row headers (if there are)
if iscell(row_headers)
    num_of_row_headers = length(row_headers);
end



for i=1:rows
    % for each row:
    % print the row header (if required)
    if iscell(row_headers)
        if (num_of_row_headers < i)
            %fprintf(file_var, ',');
        else
            fprintf(file_var, '%s,', row_headers{i});
        end
    end
    % print this row of the matrix
    for j=1:cols
        if precision_string == 0
            fprintf(file_var, '%s,', num2str(matrix(i,j)));
        else 
            fprintf(file_var, '%s,', num2str(matrix(i,j), precision_string));
        end
    end
    % move to next row
    fprintf(file_var, '\n');
end

end