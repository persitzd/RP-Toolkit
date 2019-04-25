function [initial_points_report] = HPZ_Initial_Points_Risk (param1_restrictions, param2_restrictions, rho_flag, zero_rho_initial, max_starting_points) %#ok<INUSL>

% The function constructs a matrix of initial values for the parameters' 
% optimization
% The function returns a matrix with initial_points_num rows and two columns: 
%   The first column is the initial value for beta.
%   The second column is the initial value for the additional parameter 
%   (rho if CRRA, A if CARA).

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 

% rho_flag indicates what restrictions should be put on rho:
%   true   –   1 > rho >= 0
%   false  –   exp(4)-exp(-2) > rho >= 0



% Initializing the matrix
initial_points_report = zeros(max_starting_points, 2);

% The first initial point is the choice of CFGK (2007)
initial_points_report(1,1) = exp(exp(-3)) - 1;
initial_points_report(1,2) = exp(-2);

% The second initial point is some point with beta=-1 (rho=0 is arbitrary) 
initial_points_report(2,1) = -1;
initial_points_report(2,2) = 0;

if zero_rho_initial == true
    % The second initial point [ 0 0 ]
    initial_points_report(3,:) = [0 0];
else
    % a random initial point
    initial_points_report(3,:) = rand(1,2);
end

% A matrix of random numbers (random within the interval (0,1))
% each row represents the coordinates of a random point
temp = rand(max_starting_points-3, 2);

% assigning initial values for rho
if rho_flag == true
    % for the additional parameter the value is taken from (0,1)
    % if the subject has corner solutions and the CRRA is chosen as the
    % utility (not implemented)
    initial_points_report(4:max_starting_points, 2) = temp(:,2);
else
    % we take a random number from (-2 , 4), then we use exponent on it, 
    % so we get a random number from (exp(-2) , exp(4)), then we reduce 
    % exp(-2) to get a random number from (0 , ~54.5).
    % but it is not random from uniform distribution; lower numbers have 
    % higher chances.
    initial_points_report(4:max_starting_points, 2) = exp(6*temp(:,2)-2) - exp(-2);
end



if param1_restrictions(1) == 0
    % if beta is constrained to be non-negative, the value is taken
    % from (0,2)
    initial_points_report(4:max_starting_points,1) = 2*temp(:,1);
else
    % if beta is not constrained to be non-negative, the value is taken
    % from (-1,2)
    initial_points_report(4:max_starting_points,1) = 3*temp(:,1) - 1;
end  



end


