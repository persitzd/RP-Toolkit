function print_table_to_file (file_name, table)
% this function prints a specified 2D table to a new CSV file

% file_name - name for the new CSV file (not including the '.csv' in the end) 
% table - the table (2 dimensions) to be printed

% creating the new file
% (if the file already exists, it will overwrite it; if it exists and is
% currently open, this function will fail; please note that when using this
% function)
file_var = fopen(strcat(file_name, '.csv'), 'w');
% make sure the file closes when execution stops for any reason 
cleanup_close_file = onCleanup(@() fclose(file_var));



% size of table
[rows , cols] = size(table);



% print the table headers
table_headers = table.Properties.VariableNames;
num_of_table_headers = length(table_headers);
for i=1:num_of_table_headers
    fprintf(file_var, '%s,', table_headers{i});
end
fprintf(file_var, '\n');



for i=1:rows
    % for each row:
    % print this row of the table
    for j=1:cols
        current_element = table{i,j};
        if ~isempty(current_element)
            if iscell(current_element)
                current_element = current_element{1};
            end
            if isnumeric(current_element)
                current_element = num2str(current_element);
            end
            fprintf(file_var, '%s,', current_element);
        end
    end
    % move to next row
    fprintf(file_var, '\n');
end

end