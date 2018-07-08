function [ new_observations , divide_numbers ] = HPZ_CRRA_Quantities_Transformation( observations )

% This function transforms the observations to prevent a computational problem 
% that may arise in the MMI analytical estimation in CRRA. 
% The problem, the transformation, its validity and its reasoning are explained below.

% the function returns:
%   1) the transformed observations
%   2) an array of divisors used in creating the transformed observations.
%      this array can be used to scale back to the original, pre-trasformed, observations. 

% Input
% -----
% observations is a matrice with 4 columns:
%   1 - observed quantity of x
%   2 - observed quantity of y
%   3 - price of x
%   4 - price of y
%   (endowment is assumed to be fixed at 1)

% The computational problem
% -------------------------
% The computational problem arises when Matlab's optimization routine 
% tries to calculate the criterion for an extreme value of rho (goes to inf). 
% Since the utility function includes rho as power, then either:
%   1)  if both quantities are bigger than 1, the arguments of
%       the utility function might be so small that they'll be
%       rounded to zero, resulting with zero utility when in fact the utility is positive.
%   2)  if one (or both) quantity is smaller than 1, the
%       argument of the utility function for this quantity
%       might be so big that it will be rounded to Inf, resulting -Inf utility
%       when in fact the utility is finite.

% The transformation
% ------------------
% we modify the original quantities and prices, in the following way:
% Denote A=max(x_1,x_2) and B=min(x_1,x_2).
% if B=0, then we divide both quantities by A and multiply both prices 
% by A.
% Otherwise, we divide both quantities by B and multiply both prices by B.

% Validity
% --------
% We proved in the DA2 Implementation document that such a modification will result in an 
% equivalent problem that has the exact same solution as the original.

% Reasoning
% ---------
% This modification ensures that at most one quantity is bigger than 1, 
% and that none of the quantities is positive and smaller than 1, hence 
% these issues do not occur.

% initialization
new_observations = observations;

% number of observations
[num_obs,~] = size(observations);

% min quantitiy and max quantity in each observation
min_x = min(observations(:,1:2), [], 2);
max_x = max(observations(:,1:2), [], 2);

% which number we divide in, in each observation (either min_x or max_x)
% (initialization):
divide_numbers = zeros(length(observations), 1);

for i=1:num_obs
    
    % p = Max_Y / Max_X => p = p1 / p2
    %p = observations(i, 3) / observations(i,4);
    
    % min quantitiy and max quantity
    %min_x  = min (observations(i, 1), observations(i, 2));
    %max_x  = max (observations(i, 1), observations(i, 2));
    
    % if the smaller observed quantity is 0, we choose the other
    % quantity, otherwise we choose the minimal quantity.
    % (we assume at least one quantity is positive)
    if (min_x(i) == 0)
        divide_numbers(i) = max_x(i);
    else
        divide_numbers(i) = min_x(i);
    end
    
    % modifying the quantities
    new_observations(i, 1) = observations(i, 1) / divide_numbers(i);
    new_observations(i, 2) = observations(i, 2) / divide_numbers(i);
    
    % modifying the prices
    new_observations(i, 3) = observations(i, 3) * divide_numbers(i);
    new_observations(i, 4) = observations(i, 4) * divide_numbers(i);
    
end

end

