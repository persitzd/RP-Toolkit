function AFRIAT_Index = HPZ_Afriat_efficiency_index (expenditure, index_threshold)

% this function calculates the Afriat inconsistency index,
% by using a binary search approach ("Lion in the Desert")



num_of_iterations = 50;

AFRIAT_UPPER = 1; % The maximum value for Afriat index.

AFRIAT_LOWER = 0; % The minimum value for Afriat index.

AFRIAT = 1/2;

%[rows,cols] = size(expenditure); % rows=cols=number of observations 

for i=1:num_of_iterations
    
    % diag(diag(expenditure)) is a matrix with the same diagonal 
    % as expenditure, and zero at the cell in the i'th row and 
    % the j'th row column (i~=j).
    % af_exp has at the cell in the i'th row and the 
    % j'th column, where i doesn't equal j, expenditure(i,j), 
    % otherwise (i=j),AFRIAT*expenditure(i,i).        
    
    % fixed expenditure for this v vector (for this potential AFRIAT)
    af_exp = expenditure - (diag(diag(expenditure))*(1-AFRIAT));
   
    % calculate GARP_v violations
    [GARP_v, ~, ~, ~, ~] = GARP_based_on_expenditures(af_exp, [], index_threshold);
    
    
    %sums the values of all matrix GARP's cells. Counting the number of violations (notice 
    %that a pair of bundles x^i,x^j might be counted as two violations; one for (i,j), and 
    %the other for (j,i)).

    GARP_v_ERRORS = sum(sum(GARP_v));

    %If the GARP matrix is only zeros then the data satisfies GARP
    
    if GARP_v_ERRORS == 0
        
        % if the data satisfies GARP_v, then it narrows the range to 
        % [AFRIAT, AFRIAT_UPPER], so in the next iteration, in case it's not the 
        % last one, the algorithm will try to look for a smaller parallel shifting 
        % (which means a bigger AFRIAT where the data satisfy GARP_v).    
        
        AFRIAT_LOWER = AFRIAT;
        
    else
        
        % if the data doesn't satisfy GARP_v, then it narrows the range to 
        % [AFRIAT_LOWER, AFRIAT] because if the data doesn’t satisfy GARP_v, 
        % it won't satisfy GARP_u for every u>v
        
        AFRIAT_UPPER = AFRIAT;
        
        % improving efficiency (added 20.04.2020):
        % now when we go to lower v values (AFRIAT values),
        % any observation that wasn't involved in any violation of GARPv,
        % will not be involved in GARPv for lower values of v,
        % hence we can discard this observation.
        GARP_per_obs = sum(GARP_v,1) + sum(GARP_v,2)';
        obs_involved_in_GARP_violations = (GARP_per_obs > 0);
        expenditure = expenditure(obs_involved_in_GARP_violations , obs_involved_in_GARP_violations);
        
    end 

    AFRIAT = (1/2)*(AFRIAT_LOWER + AFRIAT_UPPER);        
end

%AFRIAT_LOWER is the lowest rate found where there are no  
%violations of GARP under a parallel shifting in that rate. 

AFRIAT = AFRIAT_LOWER;      

AFRIAT_Index = 1 - AFRIAT;

end

