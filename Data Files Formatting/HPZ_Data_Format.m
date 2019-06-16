function [mat, success, is_valid] = HPZ_Data_Format (file_path, locations)

% The function formats a data file to a matrix

% file_path is the path (including name file) to the data file to be turned to a matrix 
% locations is a 6-length vector:
%   1 - column number of Subject ID in the matrix
% 	2 - column number of observation number of the subject
% 	3 - column number of the quantity of good 1 chosen by the subject.
% 	4 - column number of the quantity of good k chosen by the subject.
%       the columns of quantities of goods 1 to k are in locations(3):locations(4).   
% 	5 - column number of the max quantity of good 1 (=1/price).
% 	6 - column number of the max quantity of good k (=1/price).
%       the columns of max quantities of goods 1 to k are in locations(5):locations(6).   

% It returns a matrix of data corresponding to the required treatment.
% The matrix has six columns:
% The first column is the subject ID.
% The second column is the observation number - 50 observations per subject
% Columns 3:(2+k) are the quantities of goods chosen by the subject.
% Columns (3+k):(2+2k) are the price of goods. 

try
	% reading the data file
    [data] = csvread(file_path);
    % formatting the matrix as desired
    %mat = [data(:,locations(1)) data(:,locations(2)) data(:,locations(3)) data(:,locations(4)) 1./data(:,locations(5)) 1./data(:,locations(6))];
    mat = [data(:, locations(1:2)), data(:, locations(3):locations(4)), 1./data(:, locations(5):locations(6))];
    success = true;
catch
    warning(char(strcat(HPZ_Constants.could_not_read_file_1, {' '}, file_path, HPZ_Constants.could_not_read_file_2)));
    mat = 0;
    success = false;
end



if success
    % check if the file is formatted correctly - if not it will print an error and will crash 
    is_valid = HPZ_check_file_format_correctness (mat);
end



end


