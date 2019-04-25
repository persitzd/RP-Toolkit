function [Criterion] = HPZ_MMI_Grid_Search (prob_x, first_param, second_param, third_param, function_flag, max_x1, max_x2, utility, treatment, pref_class, debugger_mode)

% Lion in the desert

% This function receives information regarding a single observation of a
% subject (utility, and prices in the shape of max purchasable quantities
% of each good), and general information about the subject, such as the
% proposed parameters. 
% The function finds the money metric value (v*i(D,u)) for this
% observation, assuming the given functional form and the proposed parameters.
%
% The function uses a "Lion in the desert" algorithm, that is - it performs  
% a binary search in which it starts in the middle (criterion = 0.5), and in
% every step it decides whether to go up (e.g. in the beginning: to 0.75),
% or down (e.g. 0.25). In the current version of the code, the algorithm 
% goes so for 40 times, reaching a precision of +-4.5 x 10^-13, or 
% +-0.00000000000045, which we consider good enough for our needs.
% The algorithm changes the endowment at each point; but since all the 
% other code files assume constant endowment equals to 1, it actually
% changes the prices to adjust the endowment. For example, if we want to
% check for endowment = 0.25, we just need to double prices by 4 (or divide
% them by 0.25). the variable "temp_price" is the one used for these
% adjustments.
%
% parameters:
%   prob_x - the probability of the account x (relevant for CFGK)
%   first_param -
%       for risk preferences            - beta
%       for other regarding preferences - alpha
%   second_param -
%       for CRRA -  rho
%       for CARA - (irrelevant)
%       for CES  -  rho
%   third_param -
%       for CRRA - (irrelevant)
%       for CARA -  A
%       for CES  - (irrelevant)
%   max_x1 - maximum quantity purchasable of x1 given the prices (=1/p1)
%   max_x2 - maximum quantity purchasable of x2 given the prices (=1/p2)
%   utility - the utility of the chosen bundle, given the functional form
%       and the estimated parameters of the subject.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% number of iterations in the binary search loop. 
% if nubmer of iterations = k, then the precision of the 
% result will be +-(1/2)^(k+1)
num_of_iterations = 40;

% initial upper bound of range
Criterion_UPPER = 1;

% initial lower bound of range
Criterion_LOWER = 0;

% initial middle point of range
% (Criterion_UPPER + Criterion_LOWER) / 2
Criterion = 1/2;



% the 2 parameters of the utility function - will be assigned next
param = zeros (1,2);

if pref_class == HPZ_Constants.risk_pref   % risk preferences

    % the first parameter is beta
    param(1,1) = first_param;

    if function_flag == HPZ_Constants.CRRA_func   % CRRA
    
        % the second parameter is rho
        param(1,2) = second_param;
    
    elseif function_flag == HPZ_Constants.CARA_func   % CARA
    
        % the second parameter is A
        param(1,2) = third_param;
    
    end

elseif pref_class == HPZ_Constants.OR_pref   % other regarding preferences (CES)
    
    % the first parameter is alpha
    param(1,1) = first_param;

    % the second parameter is rho
    param(1,2) = second_param;    
    
end



for i=1:num_of_iterations
    
    % change the prices (but leave their ratio unchanged) to control endowment 
    temp_price = zeros(1,2);
    temp_price(1,1) = 1 / (Criterion*max_x1);
    temp_price(1,2) = 1 / (Criterion*max_x2);

    % the optimal bundle for these prices assuming these
    % functional form and estimated parameters
    x = HPZ_NLLS_Choices_Numeric (param, temp_price, 1, treatment, function_flag, 1, pref_class, debugger_mode);

    if pref_class == HPZ_Constants.risk_pref   % risk preferences
    
        % the utility of the optimal bundle assuming these functional
        % form and estimated parameters
        % for computational reasons, for both CRRA and CARA we calculate
        % a log order-preserving-transforamtion of the utility
        u = HPZ_Risk_Log_Utility (prob_x, x, first_param, second_param, third_param, function_flag);
        
    elseif pref_class == HPZ_Constants.OR_pref   % other regarding preferences
        
        % the utility of the optimal bundle assuming these
        % functional form and estimated parameters
        u = HPZ_OR_Utility (x, first_param, second_param, function_flag);
        
    end

    
    % if utility of optimal choice (u) was lower than the
    % utility of the chosen bundle (utility), than we next try to increase
    % the endowment, to increase u (we aim for u=utility). if it is bigger,
    % we next try to decrease the endowment, from the same reasoning.
    if (utility - u) > 0
        Criterion_LOWER = Criterion;
    else
        Criterion_UPPER = Criterion;
    end 
    
    Criterion = (1/2)*(Criterion_LOWER + Criterion_UPPER);
    
end   % end of loop



end