function [criterion, param] = HPZ_MMI_Criterion (param, endowments, observations, treatment, function_flag, aggregation_flag, pref_class, numeric_flag)

% The function calculates the MMI criterion per subject. Given a
% specific functional form and prices, we look for the lowest expenditure 
% level that yields at least the same level of utility as does the observed
% choices. We aggregate these differences in a few ways. By 
% aggregation_flag we choose one of those.
% The function returns an one aggregate of waste.

% "param" is returned because it may change (be rounded) during the
% calculations

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% calculating all aggrgators at once
[max_criterion, average_criterion, sum_of_squares_criterion, param] = HPZ_MMI_Aggregates(param, endowments, observations, treatment, function_flag, pref_class, numeric_flag);

% taking the desired aggregator
if aggregation_flag == HPZ_Constants.MMI_Max
    criterion = max_criterion;
elseif aggregation_flag == HPZ_Constants.MMI_Mean
    criterion = average_criterion;
elseif aggregation_flag == HPZ_Constants.MMI_AVGSSQ
    criterion = sum_of_squares_criterion;
end

end