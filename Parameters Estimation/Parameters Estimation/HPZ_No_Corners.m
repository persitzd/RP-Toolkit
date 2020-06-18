function [altered_data] = HPZ_No_Corners (data, obs_num, first)

% The function constructs a data matrix that has no corners, by using the 
% algorithm in CFGK (2007).

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 

% input:
% data is a matrix with some number of columns, that contains:
%   a column that holds the quantity of good 1 chosen by the subject.
%   a column that holds the quantity of good 2 chosen by the subject.
% obs_num is the number of observations.
% first is the column that holds the quantity of the first commodity.
%   this function assumes that the quantity of the second
%   commodity is held in the following column (first+1).


% A "small" distance/change/difference/gap is measured by w (CFGK (2007))
w = 0.001;

% disp(['Correcting for corners with omega equals ' , num2str(w)]);

% creating a copy of the data
altered_data = data;

for i=1:obs_num
    
    % if one parameter is less than one of a thousand
    % (assuming w=0.001) of the other parameter,
    % we consider it too close to the corner,
    % therefore we set it to be exactly 
    % one of a thousand of the other observation
    
    % if x < 0.001y
    if data(i,first) < w*data(i,(first+1))
        
        altered_data(i,first) = w*data(i,(first+1));
        
    end
    
    % if y < 0.001x
    if data(i,(first+1)) < w*data(i,first)
        
        altered_data(i,(first+1)) = w*data(i,first);
        
    end

end

end


