function [utility] = HPZ_OR_Utility (x, alpha, rho, function_flag)

% The function calculates the utility function for other-regarding 
% preferences (as in Kurtz et al. (2016)), for a given bundle.

% x - a vector of quantities - x(1) is the quantity of x in the given 
%   bundle while x(2) is the quantity of y in the given bundle.
% alpha - the share parameter.
% rho - the parameter of the CES function.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



if function_flag == HPZ_Constants.CES_func
    % calculate the utility (the CES-function value) using the formula:
    % u(x,y) = [alpha*x^rho + (1-alpha)*y^rho]^(1/rho)
    %self = alpha*(x(1,1)^rho);
    %other = (1-alpha)*(x(2,1)^rho);
    %utility = (self+other)^(1/rho);
    
    % we use these alternative equivalent set of formulas to avoid Infs and NaNs 
    if (alpha == 0)
        % alpha = 0, the formula becomes simple like this:
        utility = x(2);
        
    elseif (alpha == 1)
        % alpha = 1, the formula becomes simple like this:
        utility = x(1);
        
    % (in all remaining cases ahead, alpha is in (0,1))
    
    elseif (rho == 0)
        % when rho = 0, the function is not defined. the following
        % definition is equal to the limit of the utility to 0+
        utility = x(1)^alpha + x(2)^(1-alpha);
        
    elseif (x(1) == 0 || x(2) == 0) && (rho < 0)
        % if one of the quantities is 0 and rho is negative, the utility is 0 - 
        % the minimal possible utility for this DM.
        utility = 0;
        
    elseif x(1) == 0
        % if x is 0 and rho > 0, the formula becomes simple like this:
        utility = (1-alpha)^(1/rho) * x(2);
        
    elseif x(2) == 0
        % if y is 0 and rho > 0, the formula becomes simple like this:
        utility = alpha^(1/rho) * x(1);
        
    elseif ((x(2) > x(1)) && (rho > 0)) || ((x(2) < x(1)) && (rho < 0))
        % in this case, we divide and multiply by x(2)
        utility = x(2) * ( ((alpha * ((x(1)/x(2)) ^ rho)) + (1 - alpha)) ^ (1 / rho) );
    else
        % in this case, we divide and multiply by x(1)
        utility = x(1) * ( ((alpha + (1 - alpha)*((x(2)/x(1)) ^ rho))) ^ (1 / rho) );
    end
    
end



% The code is build in a way the the utility should never be NaN.
% If it is NaN, there might be a serious problem, hence we prompt a
% detailed warning to the user
if isnan(utility)
    warning('Calculated utility is equal to NaN. Please check the source of this miscalculation.');
    fprintf('alpha = %.20g , rho = %.20g\n', alpha, rho);
    fprintf('x(1) = %.20g , x(2) = %.20g\n', x(1), x(2));
    fprintf('(alpha * ((x(1)/x(2)) ^ rho) = %.20g\n', (alpha * ((x(1)/x(2)) ^ rho)));
    fprintf('(1 - alpha)) ^ (1 / rho) = %.20g\n', (1 - alpha) ^ (1 / rho));
end



%% Redundant Code:

% The optimization algorithms sometimes collapse due to infinite values. Therefore 
% we replace it with a very big number (-4.5036e+015). 
% if isinf(utility)
%     if utility < 0
%         %warning('Calculated utility was equal to -Inf, was replaced by -(1/eps).');
%         warning('Calculated utility is equal to -Inf.');
%         %utility = -(1/eps);
%     else
%         %warning('Calculated utility was equal to Inf, was replaced by (1/eps).');
%         warning('Calculated utility is equal to +Inf.');
%         %utility = (1/eps);
%     end
% end



end

