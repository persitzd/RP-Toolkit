function [log_utility] = HPZ_Risk_Log_Utility(prob_x, x, beta, rho, A, function_flag)

% The function calculates a log transformation of the utility function for 
% a given subject in risk preferences studies (such as CFGK (2007)) 
% for a given bundle.
% for ranges of parameters for which the utility is negative 
%   (or non-positive), the function calculates  -ln(-utility).
% for ranges of parameters for which the utility is positive 
%   (or non-pegative), the function calculates  ln(utility).
% for parameter values for which the utility itself is a log-based
%   function, this function will return the regular, original, utlility.
% Therfore:
% for CRRA, the function calculates the utility for rho=1,
%   the ln(utility) for rho<1, and the -ln(-utility) for rho>1.
% for CARA, the function always calculates the -ln(-utility).

% This function is needed, because the original utility functions can get
% too easily to extremely big or small numbers, that will be mistakenly 
% rounded to 0 or Inf, resulting erroneous results in the final calculations. 
% This function is therefore used in both CRRA and in CARA in numeric
% estimation, and is partially also used in CARA in MMI analytical estimation.

% prob_x is the probability of the account x.
% x is a vector of quantities - x(1) is the quantity of x in the given 
%   lottery (bundle), while x(2) is the quantity of y in the given lottery (bundle). 
% beta is the disappointment aversion parameter (Gul (1991)).
% rho is the parameter of the CRRA function.
% A is the parameter of the CARA function.
% flag is the parameter that specifies v(x) - in CFGK (2007) - is either 
%   CRRA (flag=1) or CARA (flag=2)

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



if function_flag == HPZ_Constants.CRRA_func   % CRRA
    
    if (rho == 1)
        % we calculate the utility and not ln(utility), since the utility
        % itself is a ln function
        log_utility = HPZ_Risk_Utility(prob_x, x, beta, rho, A, function_flag);
    elseif (rho < 1)
        if (min(x) == 0)
            % when min(x) == 0 and rho < 1 :
            % ln(utility) =
            % = ln( (gamma * max_x^(1-rho) + (1-gamma) * min_x^(1-rho)) / (1-rho) ) = 
            % = ln(gamma * max_x^(1-rho) + (1-gamma) * min_x^(1-rho)) - ln(1-rho) = 
            % = ln(gamma * max_x^(1-rho)) - ln(1-rho) = 
            % = ln(gamma) + (1-rho)*ln(max_x) - ln(1-rho)
            log_utility = log(gamma) + (1-rho)*log(max(x)) - log(1-rho);
        else
            % when min(x) ~= 0 and rho < 1 :
            % ln(utility) =
            % = ln( (gamma * max_x^(1-rho) + (1-gamma) * min_x^(1-rho)) / (1-rho) ) = 
            % = ln(gamma * max_x^(1-rho) + (1-gamma) * min_x^(1-rho)) - ln(1-rho) = 
            % = ln( (gamma * (max_x/min_x)^(1-rho) + (1-gamma) * 1^(1-rho)) * min_x^(1-rho) ) - ln(1-rho) =
            % = ln(gamma * (max_x/min_x)^(1-rho) + (1-gamma)) + (1-rho)*ln(min_x) - ln(1-rho)
            log_utility = log(gamma * (max(x)/min(x))^(1-rho) + (1-gamma)) + (1-rho)*log(min(x)) - log(1-rho);
        end
    else
        % if rho > 1 , we might have a problem when rho -> Inf, if one of
        % the quantities is smaller than 1. this is in fact the case that
        % this function is needed for.
        if (gamma == 1)
            % -ln(-utility) = 
            % = -ln(-max(x)^(1-rho) / (1-rho)) =
            % = -( ln(max(x)^(1-rho)) - ln(-(1-rho)) ) =
            % = -(1-rho)*ln(max(x)) + ln(-(1-rho))
            log_utility = -(1-rho)*log(max(x)) + log(-(1-rho));
        elseif (gamma == 0)
            % -ln(-utility) = 
            % = -ln(-min(x)^(1-rho) / (1-rho)) =
            % = -( ln(min(x)^(1-rho)) - ln(-(1-rho)) ) =
            % = -(1-rho)*ln(min(x)) + ln(-(1-rho))
            log_utility = -(1-rho)*log(min(x)) + log(-(1-rho));
        else
            if (min(x) == 0)
                % when min(x) == 0 and rho > 1 and gamma:
                % 1-rho is negative, therefore 0^(1-rho) = inf, and when
                % divided by (1-rho) we get -inf. -ln(-utility) is also -inf
                log_utility = -inf;
            else
                % when min(x) ~= 0 and rho > 1 :
                % ln(-utility) =
                % = ln( -(gamma * max_x^(1-rho) + (1-gamma) * min_x^(1-rho)) / (1-rho) ) = 
                % = ln(gamma * max_x^(1-rho) + (1-gamma) * min_x^(1-rho)) - ln(rho-1) = 
                % = ln( (gamma * (max_x/min_x)^(1-rho) + (1-gamma) * 1^(1-rho)) * min_x^(1-rho) ) - ln(rho-1) =
                % = ln(gamma * (max_x/min_x)^(1-rho) + (1-gamma)) + (1-rho)*ln(min_x) - ln(rho-1)
                % -ln(-utility) =
                % = - [ ln(gamma * (max_x/min_x)^(1-rho) + (1-gamma)) + (1-rho)*ln(min_x) - ln(rho-1) ] 
                log_utility = - ( log(gamma * (max(x)/min(x))^(1-rho) + (1-gamma)) + (1-rho)*log(min(x)) - log(rho-1) );
            end
        end
        
    end
    
elseif function_flag == HPZ_Constants.CARA_func   % CARA
    
    % we do not calculate the utility as normally defined,
    % but we rather calculate -log(-utility) from computational reasons.
    % it is a monotone transformation of the original utility function,
    % therefore we can use it in the same manner
    
    % We need these if's because sometimes v_max or v_min might be equal to
    % Inf or -Inf, then we might duplicate 0 with Inf, which results NaN.
    % When gamma or (1-gamma) is 0, we assume that gamma*Inf or (1-gamma)*Inf,
    % respectively, are 0 as well.
    if gamma == 1
        
        % log_minus_utility = log( exp(-A*max(x) ) = -A*max(x)
        % we want the minus of the log of the minus utility:
        log_utility = A*max(x);
        
    elseif gamma == 0
        
        % log_minus_utility = log( exp(-A*min(x) ) = -A*min(x)
        % we want the minus of the log of the minus utility:
        log_utility = A*min(x);
        
    else
        A_threshold = eps * 2^22;   % 2^(-30)
        if A == 0
            % we don't want to allow A=0,
            % so we give it the worst utility possible
            log_utility = -inf;
        elseif A > 0 && A < A_threshold   %( (min(x) > 0 && A*min(x) < A_threshold) || (min(x) == 0 && A*max(x) < A_threshold) )
            % when A -> 0, then we can't use the formula in the "else" part
            % of this statement, since "exp(-A*(max(x)-min(x)))" will be
            % rounded to 1, while after using logarithm it should not be
            % rounded to 1. 
            % when A -> 0, we use the approximation: w(x)=ax-1
            % (instead of the exact value: w(x) = -e^(-ax))
            %
            % utility = 
            % = - gamma*exp(-A*max(x)) - (1-gamma)*exp(-A*min(x)) =
            % = gamma*(A*max(x)-1) + (1-gamma)*(A*min(x)-1) =
            % = gamma*A*max(x) + (1-gamma)*A*min(x) - 1 =
            % = A * (gamma*max(x) + (1-gamma)*min(x)) - 1
            %
            % Note that we are not really interested in getting the real
            % minus log minus utility; we are only interested in any
            % order-preserving function of the utility function.
            % In this case, since "A * (gamma*max(x) + (1-gamma)*min(x)) - 1" 
            % will usually be rounded to 1, we prefer to use the monotone
            % transformation: f(x) = x + 1, meaning:
            % f(utility) = utility + 1 = A * (gamma*max(x) + (1-gamma)*min(x))  
            log_utility = A * (gamma*max(x) + (1-gamma)*min(x));
        else
            % log_minus_utility =
            % = log( gamma*exp(-A*max(x)) + (1-gamma)*exp(-A*min(x)) ) =
            % = log( ( gamma*exp(-A*max(x))/exp(-A*min(x)) + (1-gamma) ) * exp(-A*min(x)) ) = 
            % = log( gamma*exp(-A*max(x)+A*min(x)) + (1-gamma) ) + log( exp(-A*min(x)) ) = 
            % = log( gamma*exp(-A*(max(x)-min(x))) + (1-gamma) ) - A*min(x) 
            % we want the minus of the log of the minus utility:
            log_utility = - ( log( gamma*exp(-A*(max(x)-min(x))) + (1-gamma) ) - A*min(x) );
        end
        
    end
    
end



end

