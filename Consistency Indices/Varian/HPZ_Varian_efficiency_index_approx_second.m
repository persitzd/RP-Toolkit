function [VARIAN, type, Var_in_sample_residuals, one_minus_v] = HPZ_Varian_efficiency_index_approx_second (expenditure, identical_choice, index_threshold) %#ok<INUSL>

% this approximation is the second approximation, denoted as "type 3",
% while the first approximation is "type 2" and the exact calculation is "type 1". 
type = 3;

% implementing the algorithm in Varian (1993) using the version in Alcantud
% et al (2010) (Algorithm 3). Note that the example in this paper is wrong.

% number of observations
[rows,cols] = size (expenditure);

% initialization
var = ones(rows,1);

GARPE = 1;

while GARPE

    % Line 4
    exp_var = expenditure - diag((diag(expenditure)) .* (ones(rows,1)-var));

    % The matrix RATIO has at the cell in the i'th row and the j'th column, the ratio between 
    % the value of the bundle that was chosen in observation j given the prices of observation 
    % i and the value of bundle i.    
    RATIO = exp_var ./ (diag(exp_var)*ones(rows,1)');
    
    % To test for GARP-v we will do the following: for every pair of
    % choices x and y if xRy then not yP0x. We will take RP and the transpose of
    % SDRP and multiply element by element. Every 1 correspondes to 
    % xRy and yP0x. The final matrix is the zero matrix if and only if 
    % GARP-v is satisfied. 
    [GARP_v, ~, ~ , ~] = GARP_based_on_expenditures(exp_var, [], index_threshold);
    
    
    GARP_ERRORS_v = sum(sum(GARP_v));

    %If the data satisfies GARP 
    if GARP_ERRORS_v == 0
        GARPE = 0;        
    else
        % Line 7 - Gv(x_j) is the jth column of GARP

        % Line 8 - Pert is a row vector. if the jth column is zero then Gv(x_j) is empty. 
        % otherwise it includes the number required in line 8.

            %part_mat(j,k)= (p^k*x^j/p^k*x^k)*(p^j*x^k). 
            % (p^k*x^j)/(p^k*x^k)<=1 (<1) if x^kRx^j (P) 
            % and 1=(p^j*x^j)>=(p^j*x^k) (>) if x^jRx^k (P)
            % therefore, if x^k, x^j violate GARP, then part_mat(j,k)<1
        Pert_mat = GARP_v .* (RATIO');

        for j=1:rows
            for k=1:cols
                if Pert_mat(j,k) >= 1
                    Pert_mat(j,k) = 0;
                end
            end
        end

        % Line 9
        r = max(Pert_mat);

        % v is the vector of maximum value of every column in Pert_mat(meaning the maximum
        %value for every observation). w is the vector of the row numbers in which the maximum 
        %was found for every observation 
        [v,w] = max(r);

        % Line 10
        var(w) = v * var(w);
    end    
end

one_minus_v = 1 - var;

% assigning to the result variables:
% main index
VARIAN = [max(one_minus_v), mean(one_minus_v) , (sqrt(meansqr(one_minus_v)))];
% residuals
min_var_residuals = HPZ_Consistency_Indices_In_Sample_Residuals_Calc (one_minus_v, @max);
average_var_residuals = HPZ_Consistency_Indices_In_Sample_Residuals_Calc (one_minus_v, @mean);
meanssq_var_residuals = HPZ_Consistency_Indices_In_Sample_Residuals_Calc (one_minus_v, @(x) sqrt(meansqr(x)));
Var_in_sample_residuals = [min_var_residuals , average_var_residuals , meanssq_var_residuals];

% the best 1-v for each of the 3 aggregators
one_minus_v = [one_minus_v, one_minus_v, one_minus_v];

end

