function HPZ_Variables_Documentation

% following is a documentation with explanations about the nature and
% meaning of various variables that appear in the code in multiple
% functions.
% in order not to copy-paste it in many functions, we gathered all these
% explanations here.

%% data / data_matrix
% a matrix with six columns (we assume that the endowment is fixed at 1):
%   The first column is the subject ID.
%   The second column is the observation number. 
%   The third column is the quantity of good 1 chosen by the subject.
%   The fourth column is the quantity of good 2 chosen by the subject.
%   The fifth column is the price of good 1.
%   The sixth column is the price of good 2.
% it could refer to a matrix of all the subjects, or to a sub-matrix
% containing the information for only a single subject.

%% endowment / endowments
% basically, we assume all the endowments are equal to precisely 1,
% but still in some function there is a vector of the endowment of each
% observation

%% observations
% simply a matrix containing only columns 3-6 from the data_matrix above

%% Choices
% same as observations above

%% observed_choices / observed_bundles
% % simply a matrix containing only columns 3-4 from the data_matrix above

%% prices
% % simply a matrix containing only columns 5-6 from the data_matrix above

%% obs_num 
% the number of observations of this subject



%% param
% the currently examined, the estimated, or to be estimated parameters

%% criterion
% the criterion for the given subject given the parameters and method

%% criterions
% all per-observation criteria for the given subject



%% pref_class 
% indicates the type of preferences being examined:
%   1 – risk preferences
%   2 – other regarding (OR) preferences

%% function_flag
% the type of function that turns prizes into utility:
%   if pref_class = 1 (risk preferences), then:
%       1 - CRRA.
%       2 - CARA.
%   if pref_class = 2 (OR preferences), then:
%       1 - CES.



%% metric_flag 
% indicates the aggregator for the differences between the
%   observed and predicted bundles.
%   1 - Euclidean norm.
%   2 - Geometric mean (as in choi et al(2007)). 
%   3 - normalized Euclidean norm.

%% aggregation_flag
% indicates the aggregator for the MMI criterion
%   1 - Max
%   2 - Mean
%   3 - AVGSSQ (Average Sum of Squares)



%% treatment 
% the number of treatment in CFGK (2007).

%% asymmetric_flag 
% should equal to 0 for treatments that 
%   involve equal probability of the two states.
%   in case the probability is not 50-50 
%   (as in treatments 2 and 3 of CFGK),
%   this flag should be 1 (to reconstruct CFGK results use 0).

%% fix_corners 
% indicates whether correction for corner choices should be held:
%   true  – Yes (as in choi et al(2007)).
%   false – No



%% param1_restrictions
% the restrictions on the first parameter (beta/alpha)
% e.g. [-7,20] would mean that the estimation should limit the value of the
% first parameter to be not less than -7 and not more than 20.

%% param2_restrictions
% the restrictions on the second parameter (rho/A)
% e.g. [-7,20] would mean that the estimation should limit the value of the
% second parameter to be not less than -7 and not more than 20.



%% write_all_flag
% preface: 
%   when performing an estimation, the program tries to find the 
%   absolute minimum, but the fminsearchbnd can only guarantee a local 
%   minimum, until we find min_counter different points with the same 
%   (or almost same, by some threshold) criterion.
% explanation:
%   if write_all_flag is true, we print all these points to the results
%   file. if it is false, we print only the best (smallest) criterion among
%   them.

%% active_waitbar
% whether to show a waitbar per subject (true) or not (false).
% this variable is set to true without the user given an option to cancel
% it.
% this variable's purpose is to help code developers that may want to make
% alterations and are not interested in the per-subject waitbar.

%% max_time_estimation
% time (minutes) that the user specified to spend  (per subject) on running  
% the code on his machine. 
% note that this applies only for a single estimation, meaning that when
% performing bootstrap or out-of-sample residuals, each single estimation
% in those procedures will be limited by this time, and not the whole
% bootstrap/residuals process.
% also note that it is not perfectly accurate, both because the program
% will not stop in the middle of checking a single initial point, and also
% because the final part where all criterions are calculated for the
% parameters that were found may take little more time (especially when
% running NLLS with numeric approach, since MMI takes longer)

%% min_counter
% number of convergence points that the user specified. when performing an
% estimation, the program tries to find the absolute minimum, but the
% fminsearchbnd can only guarantee a local minimum. therefore we try many
% different initial points. if min_counter different initial points reached
% to the same (or almost same, by some threshold) criterion, we conclude
% that this is probably the best possible criterion and end conclude the
% estimation for that subject

%% max_starting_points
% the maximum amount of attempts to estimate the parameters. each such
% attempt starts from a different initial point, so even if the
% optimization process is analytic, it can still reach a different result
% each time.
% the process may and actually it usually stops before all these points are
% used, because the time limitaion and the convergence points limitation
% are usually more tight.
% by default this variable is 100 for analytic approach and 20 for numeric
% approach.



%% output_file_config
% determines which criterions to print in the results file. for each of 
% the 6 criterions, 1 means to print, and 0 means to avoid printing.
% E.g. if the vector is [1,0,1,0,0,1], then it will print NLLS Euclid, 
% MMI SSQ and BI criterions, but it will not print the NLLS CFGK metric,
% the MMI mean and the MMI max criterions


%% file_handle
% The Matlab object of the results file, assigned when creating the file,
% and used when writing to the file and when closing it

%% file_name_str
% A string that is used in the name of the result file, depending on the
% method chosen (Consistency Indexes, NLLS, MMI, BI)

%% file_residuals_str
% A string that is used in the name of the residuals result files, 
% depending on the method chosen (NLLS, MMI, BI)

%% one_residuals_file
% print residuals of different subjects to one file with different sheets (true),
% or to a seperate for each subject (false). default is true.

%% file_val_str / criterion_value_str
% a cell array of length 6, each corresponding to one method of estimation:  
%   1. NLLS euclidean   2. NLLS CFGK
%   3. MMI Max          4. MMI Mean
%   5. MMI AVGSSQ       6. BI (Mean)
% all of its values are 'Value', except for the value corresponding to the
% currently chosen method, which will be set to "Criterion".
% the values of this array use in the headers of the result file -
% in the result file one can distinguish the criterion that was minimized
% from the other criterions by noticing whether it is "Value" or "Criterion". 



%% fix_endowments
% whether to fix endowments to be precisely 1, since all the calculations
% ahead assume the endowment is precisely 1. it is recommended in case due
% to rounding or other ussues, the endowments might be slightly different
% from 1 (e.g. 0.99, 1.01), it is NOT meant to fix a data set that is
% not at all meant to have an endowment of 1 (in which case, first fix the
% prices and endowments so the endowments will all be 1, only then run
% this code)

%% fix_endowments_error
% if while fixing endowments, the program notices that the original
% endowment is more than (1+fix_endowments_error) or less than 
% 1/(1+fix_endowments_error), then it prints a warning



%% bootstrap_flag
% whether to perform bootstrap (true) or not (false), depending on the
% user's decision

%% number_of_samples
% sample size when performing bootstrap. currently (11.2017) it is set to
% 1000 for analytic approach and 100 for numeric approach



%% runs_counter
% when the user defines/programs the runs one after another, this variable
% keeps track of the number of runs that were defined already + the
% currently being defined run (e.g. when definig the first run it equals 1). 
% after the user click the "run now" button, this variable does not change
% anymore, and its role changes; from that point on it is a constant that
% keeps the total number of runs that should be performed (in other
% functions except for HPZ_Interface, it is named: total_runs).

%% total_runs
% the final total number of runs that the user defined to be performed.

%% current run
% keeps track of the current run (out of total_runs)while the program 
% performs its calculations and estimations.



%% GARP_flags 
% is the user choices regarding GARP WARP and SARP,
% and it is a 7-length vector, with the following boolean elements:
%   1. calculate GARP WARP and SARP (1) or not (0)
%   2. calculate residuals for them (1) or not (0)
%   3. calculate in sample residuals (1) or not (0)
%   4. calculate out of sample residuals (1) or not (0)
%   5. calculate residuals for WARP (1) or not (0)
%   6. calculate residuals for GARP (1) or not (0)
%   7. calculate residuals for SARP (1) or not (0)

%% AFRIAT_flags 
% is the user choices regarding AFRIAT index,
% and it is a 4-length vector, with the following boolean elements:
%   1. calculate AFRIAT (1) or not (0)
%   2. calculate residuals for it (1) or not (0)
%   3. calculate in sample residuals (1) or not (0)   (always 0)
%   4. calculate out of sample residuals (1) or not (0)   (always 1 when residuals is 1)

%% VARIAN_flags 
% is the user choices regarding VARIAN index,
% and it is a 5-length vector, with the following boolean elements:
%   1. calculate VARIAN (1) or not (0)
%   2. calculate residuals for it (1) or not (0)
%   3. calculate in sample residuals (1) or not (0)
%   4. calculate out of sample residuals (1) or not (0)

%% HOUTMAN_flags 
% is the user choices regarding HOUTMAN-MAX index,
% and it is a 5-length vector, with the following boolean elements:
%   1. calculate HOUTMAN-MAX (1) or not (0)
%   2. calculate residuals for it (1) or not (0)
%   3. calculate in sample residuals (1) or not (0)   (always 0)
%   4. calculate out of sample residuals (1) or not (0)   (always 1 when residuals is 1) 


end