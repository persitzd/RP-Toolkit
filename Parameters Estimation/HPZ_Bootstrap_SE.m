function [me_param_1, me_param_2, se_param_1, se_param_2, lower_95_percent_param_1, lower_95_percent_param_2, upper_95_percent_param_1, upper_95_percent_param_2] = HPZ_Bootstrap_SE(data_matrix, choice_set_type, treatment, number_of_samples, action_flag, function_flag, param1_restrictions, param2_restrictions, fix_corners, metric_flag, aggregation_flag, asymmetric_flag, pref_class, numeric_flag, significance_level, max_time_estimation, min_counter, max_starting_points, BI_threshold, debugger_mode, active_waitbar, current_run, total_runs)

% This function repeatedly picks randomly n observations out of the n 
% observations of the subject, with replacement. That means, that the common
% set of observations that will result will contain most observations, but
% not all of them.
% (for example: if we have 8 observations of the subject, we might take the
% set: {1,1,2,3,5,7,8,8} or the set: {2,3,3,3,5,6,6,7} and etc.)
% We then take this set of observations of the subject, and estimate the 
% subject's utility function parameters assuming these are the only 
% observations available.
% The function returns the average of the parameters, the maximal and
% minimal values of them, and the standard error. In fact, this function
% can return the whole distribution of the parameters' values.
% In the results files ("Choi at al. (2007) - Results" and "Halevy et al
% (2016) Part 1 - Results") we include descriptive statistics of the
% parameter frequencies in 1000 re-samplings of each individual data set in
% every reported recovery scheme. Potentially, we could have used these
% distributions to evaluate whether the restriction can be rejected.
% However, since we do not provide any proof that these re-samplings indeed
% recover confidence sets for the parameters, we merely interpret them as a
% measure for the sensitivity of the recovered parameters to extreme
% observations.
% The importance of this function is that it allows the user to see if there are 
% few/small number of observations that cause a significant change in the
% estimation, what the regular estimation (taking all observations) cannot
% emphasize. 
% one specific reason it is important is because if the DM made some mistakes in 
% choices, or changed his mind, or from any other reason, had a single or
% few exceptional choices, we might be interested to know that it might be
% the case and maybe the DM's preferences are somewhat different than the
% result that was received from estimating while accounting for all
% observations.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% NOTE: to make the code and the notes easier to read, we assume a
% significance level of 5% in variables names and in documentation.



% we want only the main waitbar to show, and not the seperate waitbar
% for each of the 100 / 1000 (or more...) repeatitive estimations
active_sub_waitbar = false;


% number of observations of the subject
obs_num = size(data_matrix,1);

% a vector to store the results of all the samples
parameters = zeros(number_of_samples,2);

% recover parameters from actual data

% re-sampling w/ replacement

% a (obs_num x bs_samples) matrice of randomized integers between 1 and obs_num. 
% every column (with length of obs_num) will be used for 1 sample:
% the integers stored in the column will be the index numbers of the
% observations that will remain; the rest of the observations will not be
% included in the estimation for that sample.
sample_indices = ceil(obs_num * rand(obs_num, number_of_samples));

% define the waitbar
if (active_waitbar)
    waitbar_name = char(strcat(HPZ_Constants.waitbar_name_estimation, {' '}, '(', HPZ_Constants.current_run_waitbar, {' '}, num2str(current_run), {' '}, HPZ_Constants.total_runs_waitbar, {' '}, num2str(total_runs), ')'));
    waitbar_msg = char(strcat(HPZ_Constants.waitbar_recovery, {' '}, num2str(data_matrix(1,1)), {' '}, HPZ_Constants.waitbar_bootstrap));
    new_bar_val = 0;
    h_wb = wide_waitbar(new_bar_val, {waitbar_msg, ''}, waitbar_name, HPZ_Constants.waitbar_width_multiplier, [0,0.12]);
end

for i=1:number_of_samples
    % we take only some of the rows in the data, only those whose indexes
    % apear at least once in the i'th column of sample_indices
    bs_sample = data_matrix(sample_indices(:,i),:);
    
    % recover and store parameters from bs_sample 
    
    % estimation for the subset of observations
    [param, ~, ~] = HPZ_Estimation (bs_sample, obs_num, choice_set_type, action_flag, treatment, function_flag, param1_restrictions, param2_restrictions, fix_corners, metric_flag, asymmetric_flag, aggregation_flag, pref_class, numeric_flag, false, max_time_estimation, min_counter, max_starting_points, BI_threshold, debugger_mode, active_sub_waitbar, current_run, total_runs);

    % assigning the result of the i'th sample to the vector
    parameters(i,:) = param;
    
    % updating the waitbar
    if (active_waitbar)
        new_bar_val = i / number_of_samples;
        waitbar(new_bar_val, h_wb, {waitbar_msg , char(strcat({'Completed '}, num2str(i), {' Iterations out of '}, num2str(number_of_samples)))});
    end
end



% calculate mean of each column
% (column 1 = beta/alpha, column 2 = rho/A)
me = mean(parameters);
% calculate standard deviation of each column
se = std(parameters);

% mean of beta/alpha of all samples
me_param_1 = me(1,1);
% mean of rho/A of all samples
me_param_2 = me(1,2);

% standard deviation of beta/alpha of all samples
se_param_1 = se(1,1);
% standard deviation of rho/A of all samples
se_param_2 = se(1,2);



% the next code is in order to find lower and upper 95%

% ordered vector of all beta/alpha values from all samples
ordered_param_1 = sort (parameters(:,1));
% ordered vector of all rho/A values from all samples
ordered_param_2 = sort (parameters(:,2));


% "lower" and "upper" are the indexes that we will pull from the matrice we
% created, that will represent the (significance level) quantile
% and the (1-significance_level) quantile (e.g. 0.05 and 0.95 quantiles).

% in calculating the indexes for the quantiles, we use the "1 + (n-1)p"
% formula, which always results a number between 1 and n.
% we use linear interpolation between the two adjacent indexes 
% ( e.g. for 5.7, we take: x(5)*0.3 + x(6)*0.7 ).
% (what we perform is exactly R-7 in https://en.wikipedia.org/wiki/Quantile) 

% the index corresponding to the 5% quantile
lower = 1 + significance_level * (number_of_samples - 1);
% the index corresponding to the 95% quantile
upper = 1 + (1-significance_level) * (number_of_samples - 1);

% the beta/alpha value that 95% of the samples are higher than it and 5% are lower than it 
lower_95_percent_param_1 = ordered_param_1 (floor(lower)) + (lower - floor(lower)) * (ordered_param_1 (floor(lower)+1) - ordered_param_1 (floor(lower)));
% the rho/A value that 95% of the samples are higher than it and 5% are lower than it 
lower_95_percent_param_2 = ordered_param_2 (floor(lower)) + (lower - floor(lower)) * (ordered_param_2 (floor(lower)+1) - ordered_param_2 (floor(lower)));
% the beta/alpha value that 95% of the samples are lower than it and 5% are higher than it 
upper_95_percent_param_1 = ordered_param_1 (ceil(upper)) - (ceil(upper) - upper) * (ordered_param_1 (ceil(upper)) - ordered_param_1 (ceil(upper)-1));
% the rho/A value that 95% of the samples are lower than it and 5% are higher than it
upper_95_percent_param_2 = ordered_param_2 (ceil(upper)) - (ceil(upper) - upper) * (ordered_param_2 (ceil(upper)) - ordered_param_2 (ceil(upper)-1));


% close the waitbar
if (active_waitbar)
    close(h_wb);
end


end
