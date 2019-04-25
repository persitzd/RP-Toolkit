function [adjusted_data] = HPZ_Fix_Endowments_To_One(data_matrix, start_index)

% This function takes the input data, and "adjust" its bundles so that the
% endowment of the bundle will be equal to 1
% NOTE: This function should NOT be used to adjust a set of data that has
% endowments not equal to 1 on purpose. It is meant to handle a set of data
% where endowments should be equal to 1, but due to some inaccuracy they
% slightly differ from 1 (0.99/0.98/0.97).

% We assume the input data_matrix is such that:
% column no. start_index represents the subject's ID
% column no. start_index+1 represents the number of observation of this
% subject
% column no. start_index+2 represents the quantity of x1
% column no. start_index+3 represents the quantity of x2
% column no. start_index+4 represents the price of x1
% column no. start_index+5 represents the price of x2
% For example: if the data_matrix has only 4 columns representing 2 quantities and
% 2 prices, then set start_index = 1; if the first two columns are subject
% ID and observation number, and then 4 columns of quantities and prices,
% set start_index = 3.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% initializing the adjusted data
adjusted_data = data_matrix;

% calculate number of goods
[~, num_of_columns] = size(data_matrix);
num_of_goods = (num_of_columns - 2) / 2;

% the endowments (may not equal to 1)
%endowments = data_matrix(:,start_index+2) .* data_matrix(:,start_index+4) + data_matrix(:,start_index+3) .* data_matrix(:,start_index+5);
endowments = sum( data_matrix(:, (start_index+1+1):(start_index+1+num_of_goods)) .* data_matrix(:, (start_index+num_of_goods+1+1):(start_index+num_of_goods+1+num_of_goods)) , 2 );
endowments_times_goods = repmat(endowments, 1, num_of_goods);

% dividing both x1 and x2 by the endowment, thus creating a new bundle that
% has an endowment exactly equal to 1
%adjusted_data(:,start_index+2) = data_matrix(:,start_index+2) ./ endowments(:);
%adjusted_data(:,start_index+3) = data_matrix(:,start_index+3) ./ endowments(:);
adjusted_data(:, (start_index+1+1):(start_index+1+num_of_goods)) = data_matrix(:, (start_index+1+1):(start_index+1+num_of_goods)) ./ endowments_times_goods;

end

