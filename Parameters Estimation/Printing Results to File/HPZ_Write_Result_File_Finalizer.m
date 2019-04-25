function HPZ_Write_Result_File_Finalizer(file_handle, output_file_config, final_output, index, subject_index, obs_num, bootstrap_flag, print_precision)

% This function prints one row (specified by "index") of the "final_output" 
% it recieves as argument, to the results file.
% final_output is a matrix with 16 columns and it's number of rows is as the
% number of convergence points, 
% the columns contain the following data:
%   1 - estimated beta/alpha (beta for CRRA & CARA, alpha for CES)
%   2 - estimated rho/A      (rho for CRRA & CES,   A for CARA)
%   3 - NLLS Euclidean Criterion
%   4 - NLLS CFGK Criterion
%   5 - NLLS normalized-Euclidean Criterion
%   6 - MMI Max Criterion
%   7 - MMI Mean Criterion
%   8 - MMI Average Sum of Squares Criterion
%   9 - BI Criterion
%   10-17 are bootstrap results:
%       10 - mean of beta/alpha
%       11 - mean of rho/A
%       12 - standard deviation of beta/alpha
%       13 - standard deviation of rho/A
%       14 - the beta/alpha value that 95% of the samples are higher than it and 5% are lower than it
%       15 - the rho/A value that 95% of the samples are higher than it and 5% are lower than it
%       16 - the beta/alpha value that 95% of the samples are lower than it and 5% are higher than it
%       17 - the rho/a value that 95% of the samples are lower than it and 5% are higher than it

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% we create a new vector with the about-to-be-printed results, but we turn
% the results to strings with the desired level of precision
subject_final_output = final_output(index,:);


    
%% printing the results of this subject column by column to the results file

% string such as '%.7f', that determines the number of digits after the
% point that will be printed (7 digits in the example)
precision_string = strcat(',%.', num2str(print_precision), 'g');

% the first is a string, cause it is the subject's index, the others 
% are the estimated values of the 2 parameters, which are always printed
print_string = strcat('%s,%s', precision_string, precision_string);

% there are at least three strings - the first is the subject's index,
% and the second and third are the estimated parameters
fprintf(file_handle, print_string, subject_index, ...
                                    num2str(obs_num), ...
                                    subject_final_output(1), ...
                                    subject_final_output(2));

% if bootstrap was asked for, we need to add 8 more columns
if (bootstrap_flag)

    % we assume that in subject_final_output,
    % bootstrap results are 8 numbers placed in indexes 10-17

    start_index = 10;

    for i=start_index:(start_index+7)
        fprintf(file_handle, precision_string, subject_final_output(i));
    end

end


% for each of the criterions, we check whether we were asked to print
% that criterion or not
% reminder: the criterions are:
%   1. NLLS - Euclidean metric
%   2. NLLS - CFGK metric
%   3. NLLS - normalized-Euclidean metric
%   4. MMI - Max aggregator
%   5. MMI - Mean aggregator
%   6. MMI - Average Sum of Squares aggregator
%   7. BI
for i=1:length(output_file_config)

    % we assume that in subject_final_output,
    % the criterions are 7 numbers placed in indexes 3-9

    if (output_file_config(i) == 1)
        fprintf(file_handle, precision_string, subject_final_output(2+i));
    end
end


% after printing the headers line, we go down to the next line when
% we start to print the results for the first subject
fprintf(file_handle, '\n');



end