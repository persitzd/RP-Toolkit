function Print_Subject_Data_To_File (main_folder_for_results, data, precision_string, choice_set_type)

if choice_set_type == HPZ_Constants.choice_set_budget_line
    % calculation of number of goods
    [obs_num, num_of_columns] = size(data); %#ok<ASGLU>
    num_of_goods = (num_of_columns - 2) / 2;
    % the data file contained maximum quantities, but we transferred them
    % to prices (assuming income=1), now we transfer it back.
    data(:, (2 + num_of_goods + 1):(2 + 2*num_of_goods)) = 1 ./ data(:, (2 + num_of_goods + 1):(2 + 2*num_of_goods));
end

% make sure the Results directory exists, if not - create it
dir_exists = exist(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir) , 'dir');
if ~dir_exists
    mkdir(strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir));
end

% name of main results file - we make sure that there isn't an existing file by the same name    
while (true)
    file_name_str = char(strcat('Subject Data - Date', {' '}, date, {' - '}, generate_random_str));
    file_path_str = strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir, '/', file_name_str);
    if ~exist(file_path_str, 'file')
        break
    end
end

print_matrix_to_file (file_path_str, data, precision_string, [], []);

end