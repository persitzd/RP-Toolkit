function col_headers = HPZ_Estimation_Residuals_Headers (in_sample_flag, out_sample_flag, pref_class, action_flag, function_flag, metric_flag, aggregation_flag)

% this function returns a cell that contains all the column headers required 
% for the residuals files for consistency and inconsistency indices

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



%% step 1 : first we count how many headers/columns headers are needed

% column 1 is subject number, 
% column 2 is estimated value of the 1st parameter,  
% column 3 is estimated value of the 2nd parameter, 
% column 4 is the criterion when estimation with all observation held place, 
% column 5 is observation number.
num_of_columns = 5;

% counting the number of headers
if (in_sample_flag)
    num_of_columns = num_of_columns + 2;
end
if (out_sample_flag)
    num_of_columns = num_of_columns + 4;
end



% create a row cell array for headers, now that we know what length it should be
col_headers = cell(1,num_of_columns);



%% step 2 : now we enter the headers to the cell

col_headers{1} = 'Subject';
if pref_class == HPZ_Constants.risk_pref
    if function_flag == HPZ_Constants.CRRA_func
        col_headers{2} = 'Beta';
        col_headers{3} = 'Rho';
    elseif function_flag == HPZ_Constants.CARA_func
        col_headers{2} = 'Beta';
        col_headers{3} = 'A';
    end
elseif function_flag == HPZ_Constants.CRRA_func
    if function_flag == HPZ_Constants.CES_func
        col_headers{2} = 'Alpha';
        col_headers{3} = 'Rho';
    end
end
if action_flag == HPZ_Constants.NLLS_action
    if metric_flag == HPZ_Constants.euclidean_metric
        col_headers{4} = 'NLLS Euclidean Criterion';
    elseif metric_flag == HPZ_Constants.CFGK_metric
        col_headers{4} = 'NLLS CFGK Criterion';
    elseif metric_flag == HPZ_Constants.normalized_euclidean_metric
        col_headers{4} = 'NLLS normalized-Euclidean Criterion';
    end
elseif action_flag == HPZ_Constants.MMI_action
    if aggregation_flag == HPZ_Constants.MMI_Max
        col_headers{4} = 'MMI Max Criterion';
    elseif aggregation_flag == HPZ_Constants.MMI_Mean
        col_headers{4} = 'MMI Mean Criterion';
    elseif aggregation_flag == HPZ_Constants.MMI_AVGSSQ
        col_headers{4} = 'MMI AVGSSQ Criterion';
    end
elseif action_flag == HPZ_Constants.BI_action
    col_headers{4} = 'BI Criterion';
end

col_headers{5} = 'Observation';

% current (next) column to put a header to
current_col = 6;
% setting all the headers
if (in_sample_flag)
    col_headers{current_col} = 'In-Sample Component Residual';
    col_headers{current_col+1} = 'In-Sample Difference Residual';
    current_col = current_col + 2;
end
if (out_sample_flag)
    % E.g. "Alternative Beta" (the beta value that was estimated when
    % estimation ran without this observation)
    col_headers{current_col} = char(strcat('Alternative', {' '} , col_headers{2}));
    col_headers{current_col+1} = char(strcat('Alternative', {' '}, col_headers{3}));
    col_headers{current_col+2} = 'Alterntive Criterion';
    col_headers{current_col+3} = 'Out-of-Sample Residual';
    current_col = current_col + 4; %#ok<NASGU>
end



end