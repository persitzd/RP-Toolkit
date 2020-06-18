function [mat, success, is_valid] = HPZ_Data_Format_Finite_Set (file_path)

% The function formats a data file to a matrix

% file_path is the path (including name file) to the data file to be turned to a matrix 
%   column 1 is the Subject ID in the matrix.
% 	column 2 is observation number of the subject.
%   columns 3-4 are the actual choice that was made.
% 	starting column 5 and onwards, every pair of of columns (i.e. (5,6), (7,8), (9,10), ...)     
%   represent a bundle that the consumer could have chosen.
%   the number of optional bundles can differ between observations.

% It returns a matrix of data corresponding to the required treatment.
% The first column is the subject ID.
% The second column is the observation number.
% The remaining columns represent the bundles the consumer could choose from.  

try
    data_not_read_yet = true;
	% reading the data file
    [data] = csvread(file_path);
    data_not_read_yet = false;
    % formatting the matrix as desired
    %mat = [data(:,locations(1)) data(:,locations(2)) data(:,locations(3)) data(:,locations(4)) 1./data(:,locations(5)) 1./data(:,locations(6))];
    mat = data;
    [rows, cols] = size(data);
    total_num_obs = rows;
    % the number of columns should be even, not odd
    if mod(cols,2) == 1
        warning('Number of columns must be even, not odd!');
        error('Number of columns must be even, not odd!');
        %cols = cols-1; 
        %mat = mat(:, 1:cols);
    end
    % we want to make sure that every empty space will be turned to NaN.
    % we also want to make sure that the actual choice is one of the specified options.  
    for obs = 1:total_num_obs
        actual_choice_was_optional = false;
        for c = 5:2:cols
            if (isempty(mat(obs,c)) || mat(obs,c) == 0) && (isempty(mat(obs,c+1)) || mat(obs,c+1) == 0)
                mat(obs, c:(c+1)) = nan;
            elseif mat(obs,c) == mat(obs,3) && mat(obs,c+1) == mat(obs,4)
                actual_choice_was_optional = true;
            end
        end
        if actual_choice_was_optional == false
            warning('Actual Choice is not one of the options presented to the consumer');
            error('Actual Choice is not one of the options presented to the consumer'); 
        end
    end
    success = true;
catch
    if data_not_read_yet
        warning(char(strcat(HPZ_Constants.could_not_read_file_1, {' '}, file_path, HPZ_Constants.could_not_read_file_2)));
    end
    mat = 0;
    success = false;
    is_valid = false;   % just to prevent unwanted error
end



if success
    % check if the file is formatted correctly - if not it will print an error and will crash 
    is_valid = HPZ_check_file_format_correctness (mat);
end



end


