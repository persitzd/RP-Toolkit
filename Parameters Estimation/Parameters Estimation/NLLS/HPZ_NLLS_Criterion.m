function [criterion, param] = HPZ_NLLS_Criterion (param, endowments, observations, choice_set_type, treatment, function_flag, fix_corners, metric_flag, asymmetric_flag, pref_class, numeric_flag, debugger_mode)

% The function calculates the NLLS criterion (either by euclidean metric
% or by the metric used in CFGK (2007), or by normalized euclidean) 
% for a given specification of a utility function and a set of choices. 
% The function returns the value of the criterion for the specified
% functional form and metric and the given data.

% "param" is returned because it may change (be rounded) during the
% calculations

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% calculating all metrics at once
[euclidean_criterion, CFGK_criterion, normalized_euclidean_criterion, param] = HPZ_NLLS_Metrics (param, endowments, observations, choice_set_type, treatment, function_flag, fix_corners, asymmetric_flag, pref_class, numeric_flag, debugger_mode);

% taking the desired metric
if metric_flag == HPZ_Constants.euclidean_metric
    criterion = euclidean_criterion;
elseif metric_flag == HPZ_Constants.CFGK_metric
    criterion = CFGK_criterion;
elseif metric_flag == HPZ_Constants.normalized_euclidean_metric
    criterion = normalized_euclidean_criterion;
end



end