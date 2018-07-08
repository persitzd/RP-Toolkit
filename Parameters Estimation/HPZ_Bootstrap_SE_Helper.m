function [me_1, me_2, se_1, se_2, le_1, le_2, ue_1, ue_2] = HPZ_Bootstrap_SE_Helper(data, treatment, number_of_samples, action, function_flag, param1_restrictions, param2_restrictions, fix_corners, metric_flag, aggregation_flag, asymmetric_flag, pref_class, numeric_flag, significance_level, rows_out, active_waitbar, current_run, total_runs, max_time_estimation, min_counter, max_starting_points)

% this function calls HPZ_Bootstrap_SE in order to perform bootstrap on the
% subjects in the data.
% for more information about the bootstrap, see the documentation there.

% this function takes the single values that were produced by
% HPZ_Bootstrap_SE, and turns each of them to a vector of rows_out
% identicle copies of the original value.
% that is needed in the case the user defined that she wishes to see a few
% of the results and not only the optimal. in that case we show the
% bootstrap results near each of these results that relate to the same
% subject.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% get the 95% upper and lower estimates of the samples.
% (not necessarily 95%, but in general it is (1-significance_level)) 
% the samples take randomly obs_num observations from the
% subject's observations, with repeat, and find the estimation 
% for this new set of observations.
[me_param_1, me_param_2, se_param_1, se_param_2, lower_95_percent_param_1, lower_95_percent_param_2, upper_95_percent_param_1, upper_95_percent_param_2] = ...
                HPZ_Bootstrap_SE(data, treatment, number_of_samples, action, function_flag, param1_restrictions, param2_restrictions, fix_corners, metric_flag, aggregation_flag, asymmetric_flag, pref_class, numeric_flag, significance_level, active_waitbar, current_run, total_runs, max_time_estimation, min_counter, max_starting_points);

% we now need to duplicate these values, for the case the user
% want to get not only the one best estimation, but all the
% saved best estimations (in which case these statistic
% parameters will appear identically near each of them)

% mean of beta/alpha
me_1 = me_param_1 * ones(rows_out,1);

% mean of rho/A
me_2 = me_param_2 * ones(rows_out,1);

% standard deviation of beta/alpha
se_1 = se_param_1 * ones(rows_out,1);

% standard deviation of rho/A
se_2 = se_param_2 * ones(rows_out,1);

% the beta/alpha value that 95% of the samples are 
% higher than it and 5% are lower than it
le_1 = lower_95_percent_param_1 * ones(rows_out,1);

% the rho/A value that 95% of the samples are 
% higher than it and 5% are lower than it
le_2 = lower_95_percent_param_2 * ones(rows_out,1);

% the beta/alpha value that 95% of the samples are 
% lower than it and 5% are higher than it
ue_1 = upper_95_percent_param_1 * ones(rows_out,1);

% the rho/a value that 95% of the samples are 
% lower than it and 5% are higher than it
ue_2 = upper_95_percent_param_2 * ones(rows_out,1);

end

