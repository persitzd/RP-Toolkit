function [criterions] = HPZ_NLLS_Criterion_Per_Observation(observations, optimal_bundles, metric_flag)

% The function calculates the NLLS in-sample residual criterion for each observation, 
% given a set of observations and corresponding predictions.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% explanations about these can be found in: HPZ_Global_Variables 
global warnings_nan
global warnings_minus_inf
global warnings_plus_inf
global current_subject



if metric_flag == HPZ_Constants.euclidean_metric
    criterions = HPZ_NLLS_Criterion_Euclid(observations, optimal_bundles);
    %sqrt( sum((observations(:,1:2) - optimal_bundles).^2, 2));
elseif metric_flag == HPZ_Constants.CFGK_metric
    criterions = HPZ_NLLS_Criterion_Ldr(observations, optimal_bundles); 
    %(ldr_chosen - ldr_predicted).^2;
elseif metric_flag == HPZ_Constants.normalized_euclidean_metric
    criterions = HPZ_NLLS_Criterion_Euclid_normalized(observations, optimal_bundles); 
end



% we now want to check if we have any non-legitimate values, that is -
% values that are infs and nans.
% if there are any - we print a warning to the user.
positive_inf = sum(isinf(criterions) & criterions > 0);
negative_inf = sum(isinf(criterions) & criterions < 0);
nan_values = sum(isnan(criterions));
warnings_plus_inf(current_subject) = warnings_plus_inf(current_subject) + positive_inf;
warnings_minus_inf(current_subject) = warnings_minus_inf(current_subject) + negative_inf;
warnings_nan(current_subject) = warnings_nan(current_subject) + nan_values;

% Old Code (Redundant):
% obs_num = length(criterions);   % total number of values ( = number of observations) 
% no_inf_criterions = criterions(~isinf(criterions));
% if sum([positive_inf, negative_inf, nan_values]) > 0
%     warnings_cell = {'Invalid Criterion Warning', param(1), param(2), obs_num, positive_inf, negative_inf, nan_values, bigger_than_1, smaller_than_0, max_deviation};
% else
%     warnings_cell = {};
% end



end