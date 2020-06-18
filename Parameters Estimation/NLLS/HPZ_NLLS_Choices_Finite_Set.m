function [choice_matrix, param] = HPZ_NLLS_Choices_Finite_Set(param, observations, function_flag, pref_class, debugger_mode)

% This function finds and returns the optimal choice (a bundle) of the DM,
% given his/hers preferences, a finite set that the DM had to choose from, and other things.  
% In case of multiple optimal choices - it returns the optimal bundle that
% is closest (in euclidean metric) to the observed bundle.

% The output is a matrix. The number of rows is the number of observations, 
% while the number of columns is the number of goods (currently 2).
% Each row is a vector of optimal choices (one quantity for each good).

% "param" is returned because it may change (be rounded) during the
% calculations

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% number of observations
[num_obs,~] = size(observations);

% a finite set of bundles the DM could have chosen from
optional_bundles = observations(:,3:end);

% extract data of observed bundles (quantity of each good - (x1,x2)) 
% from the observation matrix
observed_bundle = observations(:,1:2);

% a threshold that determines when do different utility values are considered the same 
utility_threshold = 2^(-48);


% the output matrix initialization
choice_matrix = zeros(num_obs,2);

%% risk preferences
if pref_class == HPZ_Constants.risk_pref
    
    % we need these threshold in order to avoid precision of print problems. 
    % for more detailed, see the documentation in "HPZ_Constants".
    if abs(param(1) - (-1)) < HPZ_Constants.print_threshold
        param(1) = -1;
    end
    if abs(param(2) - 1) < HPZ_Constants.print_threshold
        param(2) = 1;
    end
    
    % prob_x is the probability of the account x (probability of y is: 1-prob_x)  
    prob_x = 0.5;
    
    %% CRRA
    if function_flag == HPZ_Constants.CRRA_func
        
        % it is more convenient to use these notations
        beta = param(1);
        rho  = param(2);
        
        % update "param" to be returned (currently not needed in CRRA)
        %param = [beta , rho];

    %% CARA
    elseif function_flag == HPZ_Constants.CARA_func
        
        % it is more convenient to use these notations
        beta = param(1);
        A  = param(2);
        
        % update "param" to be returned (currently not needed in CARA)
        %param = [beta , A];
        
%         % threshold for the CARA parameter (A). 
%         % Cases 17, 18 and 19 in the DA-2 Document refer to the case where A approaches 0. 
%         % However, the distinction between A values that are close to zero
%         % and those that are not, depend on the value of x.
%         % Therefore, the threshold in the code is on the multiplication of
%         % A and x.
%         % As a value for x we take min(x) and not max(x) because they are both in a
%         % negative power and therefore min(x) is more dominant than max(x).
%         % If A*min(x)>0 is less than this threshold, we use the cases as stated 
%         % in the document. 
%         % However if A*min(x)=0 we do the same using A*max(x).
%         A_threshold = eps * 2^22;   % 2^(-30)

    end   % end of CARA

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
    
    %% CES
    if function_flag == HPZ_Constants.CES_func

        % it is more convenient to use these notations
        alpha = param(1);
        rho = param(2);
        
        % this epsilon is the threshold for deciding when a value of alpha
        % or rho is close to enough to 0 or to 1, to be considered 0 or 1,
        % respectively
        epsilon = eps * 2^12;   % 2^(-40)
        
        % rounding rho to 0 or to 1 when needed
        if (abs(rho) < epsilon)
            rho = 0;
        elseif (abs(rho-1) < epsilon)
            rho = 1;
        end
        
        % rounding alpha to 0 or to 1 when needed
        if (alpha < epsilon)
            alpha = 0;
        elseif (1-alpha < epsilon)
            alpha = 1;
        end
        
        % update "param" to be returned
        param = [alpha , rho];
        
%         % alpha / (1 - alpha)
%         alpha_ratio = alpha / (1-alpha);
%         % Elasticity of substitution (when rho ~= 0,1)
%         EOS = 1 / (1 - rho); 
        
    end   % end of CES
    
end   % end of other regarding





% this loop goes one set (of prices and endowment) after another,
% and calculates the optimal choice for each of them
% cannot use parfor any more since x2* computation depends on x1*
for i=1:num_obs 

    % HERE CALCULATION OF choice_matrix(i,:)
    optional_bundles_i_with_nan = optional_bundles(i,:);
    optional_bundles_i_without_nan = optional_bundles_i_with_nan(~isnan(optional_bundles_i_with_nan));
    optional_bundles_i = [optional_bundles_i_without_nan(1:2:(end-1))' , optional_bundles_i_without_nan(2:2:end)'];
    utilities = nan(1, size(optional_bundles_i, 1)); % initialization
    for option_index = 1:size(optional_bundles_i, 1)
        % calculate the relevant utility for each of the options
        if pref_class == HPZ_Constants.risk_pref    % Risk Peferences
            if function_flag == HPZ_Constants.CRRA_func         % CRRA
                utilities(option_index) = HPZ_Risk_Log_Utility(prob_x, optional_bundles_i(option_index,:), beta, rho, [], HPZ_Constants.CRRA_func);
            elseif function_flag == HPZ_Constants.CARA_func     % CARA
                %option_index
                %beta
                %A
                %optional_bundles_i
                utilities(option_index) = HPZ_Risk_Log_Utility(prob_x, optional_bundles_i(option_index,:), beta, [], A, HPZ_Constants.CARA_func);
                %u = utilities(option_index)
            end
        elseif pref_class == HPZ_Constants.OR_pref  % other regarding preferences
            if function_flag == HPZ_Constants.CES_func          % CES
                utilities(option_index) = HPZ_OR_Utility (optional_bundles_i(option_index,:), alpha, rho, HPZ_Constants.CES_func);
            end
        end
    end
    best_utility = max(utilities);
    if best_utility == -inf
        % which means, all optional bundles have u = -inf.
        % this occurs in CARA when A = 0.
        % we want these cases to be rebuked, so we choose "bundles" that
        % will not be good for NLLS minimization.
        choice_matrix(i,:) = [inf,inf];
    else
        best_bundles_indexes = find(utilities >= best_utility - utility_threshold);
        if isempty(best_bundles_indexes)
            error('length(best_bundles_indexes) == 0');
        elseif length(best_bundles_indexes) == 1
            % this is the best bundle given these parameters
            choice_matrix(i,:) = optional_bundles_i(best_bundles_indexes(1),:);
        else
            % there are few best bundles - we need to choose the one
            % that is closest to the bundle that was actually chosen  
            best_bundle_with_shortest_distance = nan;
            shortest_distance_squared = inf;
            for best_index = 1:length(best_bundles_indexes)
                current_bundle = optional_bundles_i(best_bundles_indexes(best_index),:);
                current_distance_squared = sum((observed_bundle(i,:) - current_bundle).^2);
                if current_distance_squared < shortest_distance_squared
                    best_bundle_with_shortest_distance = best_index;
                    shortest_distance_squared = current_distance_squared;
                end
            end
            choice_matrix(i,:) = optional_bundles_i(best_bundles_indexes(best_bundle_with_shortest_distance),:);
        end
    end

    % warnings to consol for debugging purposes (only when debugger mode is activated) 
    if debugger_mode
        if (isnan(choice_matrix(i,1)) || isnan(choice_matrix(i,2)))
            warning('At least one of the commodities in the optimal bundle is NaN when Rho is %.20g, Beta is %.20g and i is %d.', param(1), param(2), i);
        elseif (isinf(choice_matrix(i,1)) || isinf(choice_matrix(i,2)))
            warning('At least one of the commodities in the optimal bundle is +Inf or -Inf when Rho is %.20g, Beta is %.20g and i is %d.', param(1), param(2), i);
        end
    end

end   % end of loop over all observations



end   % end of function

