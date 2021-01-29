function [max_criterion, average_criterion, sqrt_avg_sum_of_squares_criterion, param] = HPZ_MMI_Aggregates (param, endowments, observations, treatment, function_flag, pref_class, numeric_flag, debugger_mode)

% The function calculates the MMI criterion per subject. Given a
% specific functional form and prices, we look for the lowest expenditure 
% level that yields at least the same level of utility as does the observed
% choices. We aggregate these differences in a few ways.

% The function returns three aggregates of waste (one minus the minimal 
% expenditure) - max, mean and sum of squares

% "param" is returned because it may change (be rounded) during the
% calculations

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



[criterions] = HPZ_MMI_Criterion_Per_Observation (param, endowments, observations, treatment, function_flag, pref_class, numeric_flag, debugger_mode);



% the 3 aggregates - max, mean, and sum of squares

% max waste
max_criterion = max(criterions);

% average waste
average_criterion = mean(criterions);

% root square of average sum of squares waste
sqrt_avg_sum_of_squares_criterion = sqrt(Smart_MeanSqr(criterions));
% (we intentionally avoid using meansqr because it emits infs and nans. we
% prefer the program and the user to get clearly unreasonable results,
% rather than results that seem to make sense but are in fact erroneous



end

