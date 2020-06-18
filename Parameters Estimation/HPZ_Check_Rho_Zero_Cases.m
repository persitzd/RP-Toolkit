function [param, criterion] = HPZ_Check_Rho_Zero_Cases(obs_num, endowments, observations, choice_set_type, treatment, function_flag, beta_restrictions, fix_corners, metric_flag, aggregation_flag, asymmetric_flag, pref_class, action_flag, numeric_flag, BI_threshold, debugger_mode)

% we want to check the Criterion value for every pair of (beta,rho)
% such that beta=p-1 and rho=0, when p is an intermediate of two of the
% prices ratios (p>=1) that the DM was exposed to.
% We want to find the lower (hence better) criterion from these cases. 
% Note: we also check the intermediate between 1 and the lower price ratio
% (if the lower isn't 1); and we also check an "intermediate" between the
% highest ratio and Infinity, by taking arbitrarily a ratio slightly bigger
% than the highest ratio (currently e take: highest ratio + 1).
% Also note that in the case of CRRA, when beta > p-1 for every 
% obserbed p, the criterion does not depend on the value of rho at all
% (that is, rho doesn't have to be 0).

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% initialization of param
param = zeros(1,2);

% a threshold - if the ratio between two prices ratios minus 1 is less than
% this threshold, the intermediate prices ratio between them will not be addressed
% Oriel Notes - Note - this threshold was chosen specifically to handle the HPZ data. a
% smaller threshold would probably fail in HPZ (would separate prices
% ratios that should be handled as identical)
price_ratio_threshold = 0.00015;

% the restrictions on the prices ratio as result of the restrictions on beta 
prices_ratios_restrictions = beta_restrictions + 1;

if choice_set_type ~= HPZ_Constants.choice_set_finite_set
    
    % list of all different prices ratios that the subject was revealed to
    prices_ratios = nan(1, obs_num+2);
    for i=1:obs_num
        % we always take the prices ratio as a number >= 1
        % e.g. we always take 1.25 and not 0.8
        if observations(i,3) > observations(i,4)
            %prices_ratios(i+1) = data(i,5) / data(i,6);
            prices_ratios(i+1) = observations(i,3) / observations(i,4);
        else
            %prices_ratios(i+1) = data(i,6) / data(i,5);
            prices_ratios(i+1) = observations(i,4) / observations(i,3);
        end
    end
    % in case the prices ratio 1 wasn't observed, we still need to check for a
    % prices ratio between the lowest ratio and 1
    prices_ratios(1) = max(1, prices_ratios_restrictions(1));
    % we also need to check for a prices ratio bigger than the largest prices
    % ratio observed. we arbitrarily take the maximum ratio + 1 (we write here
    % +2 so the average between this and the max ratio itself will be +1)
    prices_ratios(obs_num+2) = min(max(prices_ratios)+2, prices_ratios_restrictions(2));
    
else % i.e. choice_set_type == HPZ_Constants.choice_set_finite_set
    
    % in this case, we may not be able to decide properly which beta values
    % we should check, so we arbitrarily check:
    % beta = 0.1,0.2,0.5,1,1.2,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10
    prices_ratios = 1+[0.1,0.2,0.5,1,1.2,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10];
    
end

% now we sort it, and delete double entries
prices_ratios_by_order = sort(unique(prices_ratios));

% we still need to take the restrictions on beta into account, by deleting the prices that are out of range 
prices_ratios_by_order = prices_ratios_by_order(prices_ratios_by_order >= prices_ratios_restrictions(1) & prices_ratios_by_order <= prices_ratios_restrictions(2));

% number of different prices ratios
num_of_prices_ratios = length(prices_ratios_by_order);


% list of prices ratios that are intermediates of the prices ratios that
% the subject was revealed to

intermediate_prices_ratios_by_order = zeros(1, num_of_prices_ratios-1);

for k=1:(num_of_prices_ratios-1)
    % we only want to address intermediate prices ratios that differ from
    % another enough to be considered "different"; by at least some threshold
    if (prices_ratios_by_order(k+1) / prices_ratios_by_order(k)) - 1 > price_ratio_threshold
        % the intermediate prices ratio
        intermediate_prices_ratios_by_order(k) = (prices_ratios_by_order(k) + prices_ratios_by_order(k+1)) / 2;
    else
        % we arbitrarily set it to 0, it can be 1 or another number in (0,1) 
        intermediate_prices_ratios_by_order(k) = 0; 
    end
end

% now we get rid of those indexes of intermediate prices ratios that
% were less than a threshold apart
intermediate_prices_ratios_by_order = unique(intermediate_prices_ratios_by_order);
num_of_intermediate_prices_ratios = length(intermediate_prices_ratios_by_order);



% now we calculate and evaluate the Criterion for all (beta,rho) pairs of
% the type (p-1,0) for all intermediate prices ratios (p) in the vector we
% created

criterion = Inf;

for k=1:(num_of_intermediate_prices_ratios)
    % calculating the function value for the point (p-1,0)
    if (action_flag == HPZ_Constants.NLLS_action)
        % NLLS
        [criterion_temp] = HPZ_NLLS_Criterion([intermediate_prices_ratios_by_order(k)-1 , 0], endowments, observations, choice_set_type, treatment, function_flag, fix_corners, metric_flag, asymmetric_flag, pref_class, numeric_flag, debugger_mode);
    elseif (action_flag == HPZ_Constants.MMI_action)
         % MMI
        [criterion_temp] = HPZ_MMI_Criterion([intermediate_prices_ratios_by_order(k)-1 , 0], endowments, observations, treatment, function_flag, aggregation_flag, pref_class, numeric_flag, debugger_mode);
    elseif (action_flag == HPZ_Constants.BI_action)
        % BI
        [criterion_temp] = HPZ_BI_Criterion([intermediate_prices_ratios_by_order(k)-1 , 0], endowments, observations, treatment, function_flag, pref_class, numeric_flag, BI_threshold, debugger_mode);
    end
    
    if (k == 1) || (criterion_temp < criterion)
        % in the first time, and also in the next times if we got a 
        % better estimation, we enter the new result to all
        % the relevant variables and matrices
        param(1:2) = [intermediate_prices_ratios_by_order(k)-1, 0];
        criterion = criterion_temp;
    end
end
    
    
end

