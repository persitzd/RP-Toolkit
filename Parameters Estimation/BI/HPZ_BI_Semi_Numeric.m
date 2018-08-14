function [criterions, param] = HPZ_BI_Semi_Numeric(param, observations, function_flag, pref_class)

% For each observation, the function finds the optimal choice given the
% given endowment, then checks if this optimal choice has a better utility
% than the chosen one; if it is better - than the criterion for this
% observation is 1, otherwise it is 0.
% for numerical reasons, we do not take the actual endowment = 1,
% but we take endowment = 1 - HPZ_Constants.BI_threshold
% in order to make sure that when we determine the observations to have a
% criterion of 1, it is really so.

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
    end
end



% number of observations for this subject
[num_obs,~] = size(observations);

% initialization of criterions vector
criterions = zeros(num_obs, 1);



 % the endowment we check
Endowment = 1 - HPZ_Constants.BI_threshold;

% change the prices (but leave their ratio unchanged) to control endowment
% (we also change in respective the choices (the quantities), since
% in some rare cases (e.g. CRRA when beta=rho=0 & p1=p2), all
% choices on the budget line are optimal, and the choices function
% therefore returns the chosen bundle)
temp_observations = [observations(:,1:2)*Endowment , observations(:,3:4)/Endowment];

% the optimal bundle for these prices assuming these
% functional form and estimated parameters
optimal_choices = HPZ_NLLS_Choices_Analytic(param, temp_observations, function_flag, pref_class);



% loop over all the observations
for i=1:num_obs
    
%     max_x1 = 1 / observations(i,3); % compute the intersection point (1/p1)
%     max_x2 = 1 / observations(i,4); % compute the intersection point (1/p2)
%     % change the prices (but leave their ratio unchanged) to control endowment 
%     temp_prices = zeros(1,2);
%     temp_prices(1,1) = 1 / (Endowment*max_x1);
%     temp_prices(1,2) = 1 / (Endowment*max_x2);
%     temp_observations = [observations(i,1:2) , temp_prices];

%     % the optimal bundle for these prices assuming these
%     % functional form and estimated parameters
%     %x = HPZ_NLLS_Choices_Numeric (param, temp_prices, 1, treatment, function_flag, 1, pref_class);
%     x = HPZ_NLLS_Choices_Analytic(param, temp_observations(i,:), function_flag, pref_class);
    

    % compute % the utility of the optimal bundle assuming these 
    % functional form and estimated parameters
    if pref_class == HPZ_Constants.risk_pref
        % for computational reasons, for both CRRA and CARA we calculate
        % a log order-preserving-transforamtion of the utility
        
        % for the chosen bundle
        utility = HPZ_Risk_Log_Utility (prob_x, observations(i,1:2), param(1), param(2), param(2), function_flag); %#ok<*PFBNS>
        % for the optimal bundle
        u = HPZ_Risk_Log_Utility (prob_x, optimal_choices(i,1:2), param(1), param(2), param(2), function_flag);
        
    elseif pref_class == HPZ_Constants.OR_pref
        % for the chosen bundle
        utility = HPZ_OR_Utility (observations(i,1:2), param(1), param(2), function_flag);
        % for the optimal bundle
        u = HPZ_OR_Utility (optimal_choices(i,1:2), param(1), param(2), function_flag);
    end


    % if the utility of optimal choice (u) was lower than the
    % utility of the chosen bundle (utility), than we know that the BI
    % criterion is 0, otherwise it is 1
    if utility > u
        criterions(i) = 0;
    else
        criterions(i) = 1;
    end 

    
end   % end of loop over observations





end