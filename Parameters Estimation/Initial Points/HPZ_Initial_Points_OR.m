function [initial_points_report] = HPZ_Initial_Points_OR (max_starting_points)

% The function constructs a matrix of initial values for the parameters 
% optimization.
% The function returns a matrix with initial_points_num rows and two columns: 
% The first column is the initial value for alpha.
% The second column is the initial value for rho.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% Initializing the matrix
initial_points_report = zeros(max_starting_points, 2);

% A matrix of random numbers (random within the interval (0,1))
% each row represents the coordinates of a random point
temp = rand(max_starting_points-2, 2);
%temp = rand(max_starting_points, 2);

% we first test the 2 most extreme cases - alpha = 0 and alpha = 1
% (the value of rho in these 2 cases doesn't matter - we arbitrarily decided rho = 0) 
initial_points_report(1,:) = [0 , 0];   % extreme altruism
initial_points_report(2,:) = [1 , 0];   % extreme selfishness

% for alpha (between 0 and 1)
initial_points_report(3:max_starting_points,1) = temp(:,1);

% at least one initial point with a negative rho
initial_points_report(3,2) = (-0.5)*temp(1,2)-0.5;
% at least one initial point with a positive rho
initial_points_report(4,2) = 0.5*temp(2,2)+0.5;
% the rest are random between -2 and 1
initial_points_report(5:max_starting_points,2) = 3*temp(3:(max_starting_points-2),2) - 2;

end


