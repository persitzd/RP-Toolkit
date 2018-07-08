function [c,ceq] = HPZ_OR_Utility_Constraint (x, alpha, rho, function_flag, utility)

% The function generates the inequality constraint for the mincon function
% in the MMI procedure.
% prob_x is the probability of the account x.
% x is a vector of quantities - x(1) is the quantity of x in the given 
% lottery while x(2) is the quantity of y in the given lottery.
% alpha is the share parameter.
% rho is the parameter of the CES function.
% utility is the level of utility in the observation.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



%% Equality constraint!

ceq = utility - HPZ_OR_Utility (x, alpha, rho, function_flag);

if (isinf(ceq) || isnan(ceq))
    disp('oi-HPZ_OR_Utility_Constraint');
end

c = []; % set nonequality constraint to null

end