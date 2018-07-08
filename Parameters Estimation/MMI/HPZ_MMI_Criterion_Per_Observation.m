function [criterions, param] = HPZ_MMI_Criterion_Per_Observation (param, endowments, observations, treatment, function_flag, pref_class, numeric_flag)

% The function calculates the MMI criterion per subject. Given a
% specific functional form and prices, we look for the lowest expenditure 
% level that yields at least the same level of utility as does the observed
% choices. We aggregate these differences in a few ways.

% "param" is returned because it may change (be rounded) during the
% calculations

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% explanations about these can be found in: HPZ_Global_Variables 
global warnings_nan
global warnings_minus_inf
global warnings_plus_inf
global warnings_bigger_than_1
global warnings_smaller_than_0
global current_subject



% % turn all warnings on
% warning ('on','all');
% % set off warning when adding a sheet using xlswrite
% warning('off','MATLAB:xlswrite:AddSheet');



if (pref_class == HPZ_Constants.risk_pref) && (function_flag == HPZ_Constants.CRRA_func) && (numeric_flag == false)
    % if we assume CRRA, we perform a certain modification for the sake of
    % avoiding a certain computational problem.
    % this modification is needed in the calculation of the utility in its
    % CRRA form, but it helps only if we know ahead of time the quantities
    % for which we calculate the utility, therefore it is useful only in 
    % the calculations of analytic estimation, that uses the utility of the
    % observed bundle in the formulas.
    % for more detail, see the documentation of this function.
    % (dealing with the problem in numeric calculation is being done by
    % using HPZ_Log_Utility instead of HPZ_Utility)
    observations = HPZ_CRRA_Quantities_Transformation (observations);
end



if (numeric_flag)
    % numeric approach
    [criterions, param] = HPZ_MMI_Numeric(param, endowments, observations, treatment, function_flag, pref_class);
else
    % analytic approach
    [criterions, param] = HPZ_MMI_Analytic(param, observations, function_flag, pref_class);
end



% we now want to check if we have any non-legitimate values, that is -
% values that are bigger than 1 or smaller than 0, as well as infs and nans.
% if there are any - we print a warning to the user.
bigger_than_1 = sum(criterions > 1 + HPZ_Constants.MMI_threshold);   % greater than 1 values
smaller_than_0 = sum(criterions < 0 - HPZ_Constants.MMI_threshold);   % smaller than 0 values
positive_inf = sum(isinf(criterions) & criterions > 0);
negative_inf = sum(isinf(criterions) & criterions < 0);
nan_values = sum(isnan(criterions));
warnings_bigger_than_1(current_subject) = warnings_bigger_than_1(current_subject) + bigger_than_1;
warnings_smaller_than_0(current_subject) = warnings_smaller_than_0(current_subject) + smaller_than_0;
warnings_plus_inf(current_subject) = warnings_plus_inf(current_subject) + positive_inf;
warnings_minus_inf(current_subject) = warnings_minus_inf(current_subject) + negative_inf;
warnings_nan(current_subject) = warnings_nan(current_subject) + nan_values;

% Old Code (Redundant):
% obs_num = length(criterions);   % total number of values ( = number of observations) 
% no_inf_criterions = criterions(~isinf(criterions));
% max_deviation = max(max(no_inf_criterions) - 1 , - min(no_inf_criterions));   % greatest non-inf deviation from legitimate interval ([0,1]) 
% if sum([bigger_than_1, smaller_than_0, positive_inf, negative_inf, nan_values]) > 0
%     warnings_cell = {'Invalid Criterion Warning', param(1), param(2), obs_num, positive_inf, negative_inf, nan_values, bigger_than_1, smaller_than_0, max_deviation};
% else
%     warnings_cell = {};
% end



end

