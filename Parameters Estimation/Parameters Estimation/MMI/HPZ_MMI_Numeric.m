function [criterions, param] = HPZ_MMI_Numeric(param, endowments, observations, treatment, function_flag, pref_class, debugger_mode)

% The function calculates  and returns the MMI criterion per observation. 
% Given a specific functional form and prices, we look for the lowest 
% expenditure level that yields at least the same level of utility as does 
% the observed choice.

% The function returns the value of the level of expenditure for 
% each observaion, therefore criterion is a vector with length equal to the 
% number of observations.

% "param" is returned because it may change (be rounded) during the
% calculations

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% number of observations for this subject
[num_obs,~] = size(observations);

% choose proper algorithm for the fminsearchbnd optimization procedure
% the algorithm is for interior points
% not displaying error messages
options = optimset('Algorithm','interior-point','Display','off');

% prices data
prices = observations(:,3:4);

% this code is required for the CFGK (2007) data, that has
% experiments in which the probability for the two states of the world
% weren't necessarily 1:1 but it was either 50%:50%, 33%:67% or 67%:33%.
if pref_class == HPZ_Constants.risk_pref   % risk preferences
    if ((treatment == 1) || (treatment == 4))
        % equal probability - 50%:50%
        prob_x = 1/2;
    elseif treatment == 2
        % non-equal probability - 33.3%:66.7%
        prob_x = 1/3;
    elseif treatment == 3
        % non-equal probability - 66.7%:33.3%
        prob_x = 2/3;
    end
end


% a vector to store the MMI criterion values for each observation 
MMI_values = zeros(num_obs,1);


if pref_class == HPZ_Constants.risk_pref   % risk preferences

    % we need these threshold in order to avoid precision of print problems. 
    % for more detailed, see the documentation in "HPZ_Constants".
    if abs(param(1) - (-1)) < HPZ_Constants.print_threshold
        param(1) = -1;
    end
    if abs(param(2) - 1) < HPZ_Constants.print_threshold
        param(2) = 1;
    end
    
    if function_flag == HPZ_Constants.CRRA_func   % CRRA
        
        parfor i=1:num_obs
            % compute the utility level
            utility = HPZ_Risk_Log_Utility (prob_x, observations(i,1:2), param(1), param(2), 0, function_flag); %#ok<*PFBNS>
            
            if ((observations(i,1) == 0 || observations(i,2) == 0))
                
                if param(2) >= 1
                    
                    % 1-rho is negative,
                    % therefore 0^(1-rho) equals Inf,
                    % therefore the utility is -Inf (minus Inf),
                    % cause we divide by 1-rho that is negative
                    % obviously, endowments=0 is enough to achieve -Inf utility... 
                    MMI_values(i,1) = 0;
                    
                else
                    
                    % if CRRA is chosen as the utility function, then
                    % go through the possible corner solutions,
                    % and use the grid search to find the minimal 
                    % expenditure for that corner solution:
                    max_x1 = 1 / observations(i,3); % compute the intersection point (1/p1)
                    max_x2 = 1 / observations(i,4); % compute the intersection point (1/p2)
                    % compute minimal expenditure using the grid search
                    MMI_values(i,1) = HPZ_MMI_Grid_Search (prob_x, param(1), param(2), 0, function_flag, max_x1, max_x2, utility, treatment, pref_class, debugger_mode);
                    
                end
                
            else
                % for non-corner solutions:
                % use the numeric approach in order to find the minimal
                % expenditure:
                
                init_grid = zeros(2,5);
                
                x_1_init = endowments(i) / prices(i,1);
                
                % first point
                init_grid(1,1) = (102/(1000+1)) * x_1_init;
                init_grid(1,2) = ( endowments(i) - (init_grid(1,1)*prices(i,1)) ) / prices(i,2);
                [init_grid(1,3:4), init_grid(1,5), eflag1] = fmincon (@(x) prices(i,1:2)*x', init_grid(1,1:2), prices(i,1:2), endowments(i), [], [], [0 0], [], @(x) HPZ_Risk_Utility_Constraint(prob_x, x, param(1), param(2), 0, function_flag, utility), options);
                
                % second point
                init_grid(2,1) = (898/(1000+1)) * x_1_init;
                init_grid(2,2) = ( endowments(i) - (init_grid(2,1)*prices(i,1)) ) / prices(i,2);
                [init_grid(2,3:4), init_grid(2,5), eflag2] = fmincon (@(x) prices(i,1:2)*x', init_grid(2,1:2), prices(i,1:2), endowments(i), [], [], [0 0], [], @(x) HPZ_Risk_Utility_Constraint(prob_x, x, param(1), param(2), 0, function_flag, utility), options);
                
                if eflag1 == 1 && eflag2 == 1
                    % if both fmincon estimations were
                    % successful, that is they ended because they reached
                    % the minimum (and not becuase of another reason such
                    % as reached max number of iterations and etc.)
                    
                    % compute minimal expenditure
                    MMI_values(i,1) = min(init_grid(:,5));
                else
                    max_x1 = 1 / observations(i,3); % compute the intersection point (1/p1)
                    max_x2 = 1 / observations(i,4); % compute the intersection point (1/p2)
                    % compute minimal expenditure using the grid search
                    MMI_values(i,1) = HPZ_MMI_Grid_Search (prob_x, param(1), param(2), 0, function_flag, max_x1, max_x2, utility, treatment, pref_class, debugger_mode);
                end
                
            end
            
            % warnings to consol for debugging purposes (only when debugger mode is activated) 
            if debugger_mode
                if isnan(MMI_values(i,1))
                    warning('Criterion is NaN when Rho is %.20g, Beta is %.20g, p1 is %.20g, p2 is %.20g, utility is %.20g, x1 is %.20g, x2 is %.20g and i is %d.', param(2), param(1), prices(i,1), prices(i,2), utility, observations(i,1), observations(i,2), i);
                elseif isinf(MMI_values(i,1))
                    % note that the MMI Value is the opposite of the
                    % criterion: Criterion = 1 - MMI_Value
                    if (MMI_values(i,1) > 0)
                        inf_str = '-';
                    else
                        inf_str = '+';
                    end
                    warning('Criterion is %sInf when Rho is %.20g, Beta is %.20g, p1 is %.20g, p2 is %.20g, utility is %.20g, x1 is %.20g, x2 is %.20g and i is %d.', inf_str, param(2), param(1), prices(i,1), prices(i,2), utility, observations(i,1), observations(i,2), i);
                elseif (MMI_values(i,1) < -HPZ_Constants.MMI_threshold) || (MMI_values(i,1) > 1+HPZ_Constants.MMI_threshold)
                    warning('Criterion is not in the range [0,1]: Criterion is %.20g when Rho is %.20g, Beta is %.20g, p1 is %.20g, p2 is %.20g, utility is %.20g, x1 is %.20g, x2 is %.20g and i is %d.', 1-MMI_values(i,1), param(2), param(1), prices(i,1), prices(i,2), utility, observations(i,1), observations(i,2), i);
                end
            end
            
        end   % end of loop (parfor)
        
    elseif function_flag == HPZ_Constants.CARA_func   % CARA
        
        parfor i=1:num_obs
            % compute the utility level
            utility = HPZ_Risk_Log_Utility (prob_x, observations(i,1:2), param(1), 0, param(2), function_flag);
            
            % use the numeric approach in order to find the minimal expenditure:
            % for both corner and non-corner solutions
            init_grid = zeros(2,5);
            
            x_1_init = endowments(i) / prices(i,1);
            
            % first point
            init_grid(1,1) = (102/(1000+1)) * x_1_init;
            init_grid(1,2) = ( endowments(i) - (init_grid(1,1)*prices(i,1)) ) / prices(i,2);
            [init_grid(1,3:4), init_grid(1,5), eflag1] = fmincon (@(x) prices(i,1:2)*x', init_grid(1,1:2), prices(i,1:2), endowments(i), [], [], [0 0], [], @(x) HPZ_Risk_Utility_Constraint(prob_x, x, param(1), 0, param(2), function_flag, utility), options);
            
            % second point
            init_grid(2,1) = (898/(1000+1)) * x_1_init;
            init_grid(2,2) = ( endowments(i) - (init_grid(2,1)*prices(i,1)) ) / prices(i,2);
            [init_grid(2,3:4), init_grid(2,5), eflag2] = fmincon (@(x) prices(i,1:2)*x', init_grid(2,1:2), prices(i,1:2), endowments(i), [], [], [0 0], [], @(x) HPZ_Risk_Utility_Constraint(prob_x, x, param(1), 0, param(2), function_flag, utility), options);
            
            if eflag1 == 1 && eflag2 == 1
                % compute minimal expenditure
                MMI_values(i,1) = min(init_grid(:,5));
            else
                max_x1 = 1 / observations(i,3); % compute the intersection point (1/p1)
                max_x2 = 1 / observations(i,4); % compute the intersection point (1/p2)
                % compute minimal expenditure using the grid search
                MMI_values(i,1) = HPZ_MMI_Grid_Search (prob_x, param(1), 0, param(2), function_flag, max_x1, max_x2, utility, treatment, pref_class, debugger_mode);
            end
            
            % warnings to consol for debugging purposes (only when debugger mode is activated) 
            if debugger_mode
                if isnan(MMI_values(i,1))
                    warning('Criterion is NaN when A is %.20g, Beta is %.20g, p1 is %.20g, p2 is %.20g, utility is %.20g, x1 is %.20g, x2 is %.20g and i is %d.', param(2), param(1), prices(i,1), prices(i,2), utility, observations(i,1), observations(i,2), i);
                elseif isinf(MMI_values(i,1))
                    % note that the MMI Value is the opposite of the
                    % criterion: Criterion = 1 - MMI_Value
                    if (MMI_values(i,1) > 0)
                        inf_str = '-';
                    else
                        inf_str = '+';
                    end
                    warning('Criterion is %sInf when A is %.20g, Beta is %.20g, p1 is %.20g, p2 is %.20g, utility is %.20g, x1 is %.20g, x2 is %.20g and i is %d.', inf_str, param(2), param(1), prices(i,1), prices(i,2), utility, observations(i,1), observations(i,2), i);
                elseif (MMI_values(i,1) < -HPZ_Constants.MMI_threshold) || (MMI_values(i,1) > 1+HPZ_Constants.MMI_threshold)
                    warning('Criterion is not in the range [0,1]: Criterion is %.20g when A is %.20g, Beta is %.20g, p1 is %.20g, p2 is %.20g, utility is %.20g, x1 is %.20g, x2 is %.20g and i is %d.', 1-MMI_values(i,1), param(2), param(1), prices(i,1), prices(i,2), utility, observations(i,1), observations(i,2), i);
                end
            end
            
        end   % end of loop (parfor)
        
    end   % end of: if (CRRA) else (CARA)

elseif pref_class == HPZ_Constants.OR_pref   % other regarding preferences
    
    % we need these threshold in order to avoid precision of print problems. 
    % for more detailed, see the documentation in "HPZ_Constants".
    if abs(param(1) - 1) < HPZ_Constants.print_threshold
        param(1) = -1;
    end
    if abs(param(2) - 1) < HPZ_Constants.print_threshold
        param(2) = 1;
    end
    
    if function_flag == HPZ_Constants.CES_func   % CES 
    
        parfor i=1:num_obs
            % compute the utility level
            utility = HPZ_OR_Utility (observations(i,1:2), param(1), param(2), function_flag);
            
            % use the numeric approach in order to find the minimal
            % expenditure:
            
            init_grid = zeros(2,5);
            
            x_1_init = endowments(i) / prices(i,1);
            
            % first point
            init_grid(1,1) = (102/(1000+1)) * x_1_init;
            init_grid(1,2) = ( endowments(i) - (init_grid(1,1)*prices(i,1)) ) / prices(i,2);
            [init_grid(1,3:4), init_grid(1,5), eflag1] = fmincon (@(x) prices(i,1:2)*x', init_grid(1,1:2), prices(i,1:2), endowments(i), [], [], [0 0], [], @(x) HPZ_OR_Utility_Constraint(x, param(1), param(2), function_flag, utility), options);
            
            % second point
            init_grid(2,1) = (898/(1000+1)) * x_1_init;
            init_grid(2,2) = ( endowments(i) - (init_grid(2,1)*prices(i,1)) ) / prices(i,2);
            [init_grid(2,3:4), init_grid(2,5), eflag2] = fmincon (@(x) prices(i,1:2)*x', init_grid(2,1:2), prices(i,1:2), endowments(i), [], [], [0 0], [], @(x) HPZ_OR_Utility_Constraint(x, param(1), param(2), function_flag, utility), options);
            
            if eflag1 == 1 && eflag2 == 1
                % compute minimal expenditure
                MMI_values(i,1) = min(init_grid(:,5));
            else
                max_x1 = 1 / observations(i,3); % compute the intersection point (1/p1)
                max_x2 = 1 / observations(i,4); % compute the intersection point (1/p2)
                % compute minimal expenditure using the grid search
                MMI_values(i,1) = HPZ_MMI_Grid_Search (0, param(1), param(2), 0, function_flag, max_x1, max_x2, utility, treatment, pref_class, debugger_mode);
            end
            
            % warnings to consol for debugging purposes (only when debugger mode is activated) 
            if debugger_mode
                if isnan(MMI_values(i,1))
                    warning('Criterion is NaN when Rho is %.20g, Alpha is %.20g, p1 is %.20g, p2 is %.20g, utility is %.20g, x1 is %.20g, x2 is %.20g and i is %d.', param(2), param(1), prices(i,1), prices(i,2), utility, observations(i,1), observations(i,2), i);
                elseif isinf(MMI_values(i,1))
                    % note that the MMI Value is the opposite of the
                    % criterion: Criterion = 1 - MMI_Value
                    if (MMI_values(i,1) > 0)
                        inf_str = '-';
                    else
                        inf_str = '+';
                    end
                    warning('Criterion is %sInf when Rho is %.20g, Alpha is %.20g, p1 is %.20g, p2 is %.20g, utility is %.20g, x1 is %.20g, x2 is %.20g and i is %d.', inf_str, param(2), param(1), prices(i,1), prices(i,2), utility, observations(i,1), observations(i,2), i);
                elseif (MMI_values(i,1) < -HPZ_Constants.MMI_threshold) || (MMI_values(i,1) > 1+HPZ_Constants.MMI_threshold)
                    warning('Criterion is not in the range [0,1]: Criterion is %.20g when Rho is %.20g, Alpha is %.20g, p1 is %.20g, p2 is %.20g, utility is %.20g, x1 is %.20g, x2 is %.20g and i is %d.', 1-MMI_values(i,1), param(2), param(1), prices(i,1), prices(i,2), utility, observations(i,1), observations(i,2), i);
                end
            end
            
        end   % end of loop (parfor)
        
    end   % end of CES

end   % end of other regarding



% MMI_values is a vector of minimal endowment (required to achieve 
% the same utility) for each observation, therefore for each observation
% the criterion is the real endowment (1) minus that minimal endowment

% in the function so far, we calculated the minimal endowment required to
% achieve the same utility, e.g. if the subject can't get this utiity with
% a smaller endowment, then MMI_values(i) == 1. but our criterion is one
% that is set to 0 when the subject chose in accordance to its preferences,
% therefore: criterions(i) = 1 - MMI_values(i)
criterions = endowments - MMI_values;



end