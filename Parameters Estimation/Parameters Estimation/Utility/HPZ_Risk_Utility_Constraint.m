function [c, ceq] = HPZ_Risk_Utility_Constraint (prob_x, x, beta, rho, A, function_flag, utility)

% The function generates the equality constraint for the mincon function
% in the MMI procedure.
% prob_x is the probability of the account x.
% x is a vector of quantities - x(1) is the quantity of x in the given 
% lottery while x(2) is the quantity of y in the given lottery.
% beta is the disappointment aversion parameter (Gul (1991)).
% rho is the parameter of the CRRA function.
% A is the parameter of the CARA function.
% utility is the level of utility in the observation.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



%% Equality constraint!

% for computational reasons, for both CRRA and CARA we calculate
% a log order-preserving-transforamtion of the utility
ceq = utility - HPZ_Risk_Log_Utility (prob_x, x, beta, rho, A, function_flag);


if (isinf(ceq) || isnan(ceq))
    disp('oi-HPZ_Utility_Constraint');
end

c = []; % set nonequality constraint to null

end