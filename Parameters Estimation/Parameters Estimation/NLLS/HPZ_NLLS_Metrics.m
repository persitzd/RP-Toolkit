function [metric_euclidean, metric_CFGK, metric_normalized_euclidean, param] = HPZ_NLLS_Metrics (param, endowments, observations, choice_set_type, treatment, function_flag, fix_corners, asymmetric_flag, pref_class, numeric_flag, debugger_mode)

% The function calculates the NLLS criterion with both euclidean metric
% and with the metric used in CFGK (2007)) for a given specification of a 
% utility function and a set of choices. 
% The function returns the values of the criterion by these 2 metrics
% for the specified functional form and metric and the given data.

% "param" is returned because it may change (be rounded) during the
% calculations

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% number of observations of this subject
obs_num = length(observations(:,1));

if choice_set_type ~= HPZ_Constants.choice_set_finite_set
    if numeric_flag == HPZ_Constants.numeric
        % numeric approach
        [optimal_bundles, param] = HPZ_NLLS_Choices_Numeric(param, observations(:,3:4), endowments, treatment, function_flag, asymmetric_flag, pref_class, debugger_mode);
    elseif numeric_flag == HPZ_Constants.analytic || numeric_flag == HPZ_Constants.semi_numeric
        % analytic approach, and also semi-numeric approach of MMI
        [optimal_bundles, param] = HPZ_NLLS_Choices_Analytic(param, observations, function_flag, pref_class, debugger_mode);
    end
else
    % finite set, go over all options
    [optimal_bundles, param] = HPZ_NLLS_Choices_Finite_Set(param, observations, function_flag, pref_class, debugger_mode);
end



if (fix_corners == true) && (choice_set_type ~= HPZ_Constants.choice_set_finite_set)
    % Choi et al. (2007) correction for corners should be held
    optimal_bundles = HPZ_No_Corners (optimal_bundles, obs_num, 1);
end



if choice_set_type ~= HPZ_Constants.choice_set_finite_set
    % compute NLLS criterion using Euclidean metric
    metric_euclidean = sum(HPZ_NLLS_Criterion_Euclid(observations, optimal_bundles));
    % compute NLLS criterion using Choi et al. (2007) metric
    metric_CFGK = sum(HPZ_NLLS_Criterion_Ldr(observations, optimal_bundles));
    % compute NLLS criterion using normalized-Euclidean metric
    % we use mean and not sum, because we want it to stay normalized as a whole 
    metric_normalized_euclidean = mean(HPZ_NLLS_Criterion_Euclid_normalized(observations, optimal_bundles));
else
    % compute NLLS criterion using Euclidean metric
    metric_euclidean = sum(HPZ_NLLS_Criterion_Euclid(observations, optimal_bundles));
    % CFGK metric and normalized euclidean cannot be used
    metric_CFGK = nan;
    metric_normalized_euclidean = nan;
end



end