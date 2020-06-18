function [criterions, param] = HPZ_MMI_Semi_Numeric(param, observations, function_flag, pref_class, debugger_mode)

% Lion in the desert

% The function finds the money metric value (v*i(D,u)) for every
% observation, assuming the given functional form and the proposed parameters.
%
% For each observation:
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

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% equal probability - 50%:50%
prob_x = 1/2;
        


% these thresholds are the same as the ones in HPZ_NLLS_Choices_Analytic
if pref_class == HPZ_Constants.risk_pref
    
    % we need these threshold in order to avoid precision of print problems. 
    % for more detailed, see the documentation in "HPZ_Constants".
    if abs(param(1) - (-1)) < HPZ_Constants.print_threshold
        param(1) = -1;
    end
    if abs(param(2) - 1) < HPZ_Constants.print_threshold
        param(2) = 1;
    end
elseif pref_class == HPZ_Constants.OR_pref   % other regarding preferences
    
    % we need these threshold in order to avoid precision of print problems. 
    % for more detailed, see the documentation in "HPZ_Constants".
    if abs(param(1) - 1) < HPZ_Constants.print_threshold
        param(1) = 1;
    end
    if abs(param(2) - 1) < HPZ_Constants.print_threshold
        param(2) = 1;
    elseif abs(param(2) - (-1)) < HPZ_Constants.print_threshold
        param(2) = -1;
    end
    
    % we need these threshold for computational reasons; in the NLLS_choices
    % function these threshold are used, so for consistency we must use them
    % when calculating the utility as well
    epsilon = eps * 2^12;   % 2^(-40)
    if (param(1) < epsilon)
        param(1) = 0;
    elseif (1-param(1) < epsilon)
        param(1) = 1;
    end
    if (abs(param(2)) < epsilon)
        param(2) = 0;
    elseif (abs(param(2)-1) < epsilon)
        param(2) = 1;
    end
end



% number of observations for this subject
[num_obs,~] = size(observations);

% initialization of criterions vector
criterions = zeros(num_obs, 1);



% loop over all the observations
for i=1:num_obs
    
%     max_x1 = 1 / observations(i,3); % compute the intersection point (1/p1)
%     max_x2 = 1 / observations(i,4); % compute the intersection point (1/p2)
    
    % compute the utility level
    if pref_class == HPZ_Constants.risk_pref
        utility = HPZ_Risk_Log_Utility (prob_x, observations(i,1:2), param(1), param(2), param(2), function_flag); %#ok<*PFBNS>
    elseif pref_class == HPZ_Constants.OR_pref
        utility = HPZ_OR_Utility (observations(i,1:2), param(1), param(2), function_flag);
    end
    
    
    
    % number of iterations in the binary search loop. 
    % if nubmer of iterations = k, then the precision of the 
    % result will be +-(1/2)^(k+1)
    num_of_iterations = 40;

    % initial upper bound of range
    Endowment_UPPER = 1;

    % initial lower bound of range
    Endowment_LOWER = 0;

    % initial middle point of range
    % (Endowment_UPPER + Endowment_LOWER) / 2
    % but in the first iteration, we check Endowment = 1
    % to check if it is perfectly consistent
    Endowment = 1; %1/2;
    
    % the binary search
    for j=1:(num_of_iterations+1)

        % change the prices (but leave their ratio unchanged) to control endowment 
        % (we also change in respective the choices (the quantities), since
        % in some rare cases (e.g. CRRA when beta=rho=0 & p1=p2), all
        % choices on the budget line are optimal, and the choices function
        % therefore returns the chosen bundle)
        temp_observations = [observations(i,1:2)*Endowment , observations(i,3:4)/Endowment];
%         temp_quantities = zeros(1,2);
%         temp_quantities(1,1) = observations(i,1) * Endowment;
%         temp_quantities(1,2) = observations(i,2) * Endowment;
%         temp_prices = zeros(1,2);
%         temp_prices(1,1) = observations(i,3) / Endowment;
%         temp_prices(1,2) = observations(i,4) / Endowment;
        
        % the optimal bundle for these prices assuming these
        % functional form and estimated parameters
        x = HPZ_NLLS_Choices_Analytic(param, temp_observations, function_flag, pref_class, debugger_mode);
        
        if pref_class == HPZ_Constants.risk_pref   % risk preferences

            % the utility of the optimal bundle assuming these functional
            % form and estimated parameters
            % for computational reasons, for both CRRA and CARA we calculate
            % a log order-preserving-transforamtion of the utility
            u = HPZ_Risk_Log_Utility (prob_x, x, param(1), param(2), param(2), function_flag);
            
        elseif pref_class == HPZ_Constants.OR_pref   % other regarding preferences

            % the utility of the optimal bundle assuming these
            % functional form and estimated parameters
            u = HPZ_OR_Utility (x, param(1), param(2), function_flag);

        end
        

        if j == 1
            % in the first interation, we check Endowment = 1, to see if
            % there is any need to reduce the endowment at all.
            if (utility - u) >= 0
                % then we found that Endowment = 1 and there is no need 
                % for any further calculations
                break
            else
                % then we prepare for the next iteration:
                Endowment = 1/2; %(1/2)*(Endowment_LOWER + Endowment_UPPER);
            end
        else
            % in all iterations but the fisrt one:
            %
            % if the utility of optimal choice (u) was lower than the
            % utility of the chosen bundle (utility), then we next try to increase
            % the endowment, in order to increase u (we aim for u=utility). if it is bigger,
            % we next try to decrease the endowment, from the same reasoning.
            %
            % we intentionally ask whehter "utility - u) > 0" and not
            % "utility - u) >= 0", because if they are equal, we want to check
            % if there is a lower endowment where the equality still holds.
            if (utility - u) > 0
                Endowment_LOWER = Endowment;
            else
                Endowment_UPPER = Endowment;
            end  

            Endowment = (1/2)*(Endowment_LOWER + Endowment_UPPER);

        end
        
    end   % end of binary search
    
    %fprintf('obs %d utility : %.20g\n', i, utility);
    %fprintf('obs %d choice : ( %.20g , %.20g ) \n', i, x(1), x(2));
    
    % assigning the result for this observation to the criterions vector
    criterions(i) = 1 - Endowment;
    
    
end   % end of loop over observations





end