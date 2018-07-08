function [criterion, param] = HPZ_NLLS_Criterion (param, endowments, observations, treatment, function_flag, fix_corners, metric_flag, asymmetric_flag, pref_class, numeric_flag)

% The function calculates the NLLS criterion (either by euclidean metric
% or by the metric used in CFGK (2007)) for a given specification of a 
% utility function and a set of choices. 
% The function returns the value of the criterion for the specified
% functional form and metric and the given data.

% "param" is returned because it may change (be rounded) during the
% calculations

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% number of observations of this subject
obs_num = length(observations(:,1));

if numeric_flag == true 
    % numeric approach
    [optimal_bundles, param] = HPZ_NLLS_Choices_Numeric(param, observations(:,3:4), endowments, treatment, function_flag, asymmetric_flag, pref_class);
else
    % analytic approach
    [optimal_bundles, param] = HPZ_NLLS_Choices_Analytic(param, observations, function_flag, pref_class);
end



if (fix_corners == true)
    % Choi et al. (2007) correction for corners should be held
    optimal_bundles = HPZ_No_Corners (optimal_bundles, obs_num, 1);
end



% all per observation criterions
criterions = HPZ_NLLS_Criterion_Per_Observation(observations, optimal_bundles, metric_flag);

% we use sum aggregation
criterion = sum(criterions);



end