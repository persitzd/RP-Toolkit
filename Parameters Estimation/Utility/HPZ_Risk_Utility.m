function [utility] = HPZ_Risk_Utility (prob_x, x, beta, rho, A, function_flag)

% The function calculates the utility function of a subject in risk 
% preferences studies (such as CFGK (2007)) for a given bundle using MMI 
% (Money Metric Index).

% The function calculates the utility function in CFGK (2007) for a given 
% bundle.

% prob_x is the probability of the account x.
% x is a vector of quantities - x(1) is the quantity of x in the given 
%   lottery (bundle), while x(2) is the quantity of y in the given lottery (bundle). 
% beta is the disappointment aversion parameter (Gul (1991)).
% rho is the parameter of the CRRA function.
% A is the parameter of the CARA function.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% By Gul (1991), in the case of two prizes, alpha is the probability of the 
% higher prize.
if x(1) > x(2)
    alpha = prob_x;
else
    alpha = 1 - prob_x;
end

% By Gul (1991), the weights are calculated using gamma(alpha)
gamma = alpha / ( 1 + ((1-alpha)*beta) );

% By Gul (1991), the disappointment aversion utility function for the case
% of two prizes is the weighted sum: 
% gamma*v(max(x,y))+(1-gamma)*v(min(x,y))

% First, let us calculate v(max(x,y)) and v(min(x,y)) 
if function_flag == HPZ_Constants.CRRA_func
    v_max = CRRA(max(x),rho);
    v_min = CRRA(min(x),rho);
elseif function_flag == HPZ_Constants.CARA_func
    v_max = CARA(max(x),A);
    v_min = CARA(min(x),A);
end

% We need these if's because sometimes v_max or v_min might be equal to
% Inf or -Inf, then we might duplicate 0 with Inf, which results NaN.
% When gamma or (1-gamma) is 0, we assume that gamma*Inf or (1-gamma)*Inf,
% respectively, are 0 as well.
if gamma == 1
    utility = gamma*v_max;
elseif gamma == 0
    utility = (1-gamma)*v_min;
else
    utility = gamma*v_max + (1-gamma)*v_min;
end



%% Redundant Code:

% if x is zero and rho is greater than 1, in CRRA we get -inf. The
% optimization algorithms sometimes collapse due to this value. Therefore 
% we replace it with a very big number (-4.5036e+015).
% if isinf(utility)
%     if utility < 0
%         warning('Calculated utility was equal to -Inf, was replaced by -(1/eps).');
%         utility = -(1/eps);
%     else
%         warning('Calculated utility was equal to Inf, was replaced by (1/eps).');
%         utility = (1/eps);
%     end
% end

end
