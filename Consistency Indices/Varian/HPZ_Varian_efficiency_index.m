function [VARIAN, Var_exact, Var_in_sample_residuals] = HPZ_Varian_efficiency_index (expenditure, identical_choice, index_threshold)

% this function calculates the Varian inconsistency index
% if the calculation seems to be taking too long, it will return an
% approximation of the Varian inconsistency index. Var_exact is a flag that
% tells whether it was an exact calculation (1) or not (2 or 3)
% either way, the function also returns the in-sample residuals of the index.



% in_sample_residuals = 0;
% average_var_residuals = 0;
% meanssq_var_residuals = 0;
% min_var_residuals = 0;

% matrix identical_choice locate identical choices. The value 
% of cell (j,k) equals 1 if the choices are identical, and 0 otherwise. 
% global identical_choice

max_exact = 26;

% rows = number of observations 
[rows, cols] = size(expenditure);

average_var = 1;
         
meanssq_var = 1;
    
min_var = 1;
    
% The matrix REF has at the cell in the i'th row and the j'th column, the difference between the value of the bundle that 
% was chosen in observation i and the bundle that was chosen in observation j given the prices of observation i
REF = diag(expenditure) * ones(rows,1)' - expenditure;

% The matrix RATIO has at the cell in the i'th row and the j'th column, the ratio between 
% the value of the bundle that was chosen in observation j given the prices of observation 
% i and the value of the bundle that was chosen in observation i 
RATIO = expenditure ./ (diag(expenditure) * ones(rows,1)');
 
% DRP - represents the relation R^0. The matrix DRP has at the cell in the i'th row and the j'th
% column 1 if and only if the bundle that was chosen in observation i is directly revealed 
% preferred to the bundle that was chosen in observation j. otherwise it equals 0.
DRP = ceil((REF+index_threshold) / (max(max(abs(REF+index_threshold)))+1));

[~ , Components] = Johnson(DRP'); % Use's Johnson's algorithm to extract list cycles (cycles without subcycles    

% SDRP - represents the relation P^0. 
% The matrix SDRP has at the cell in the i'th row and the j'th column, 1 if and only if the bundle
% that was chosen in observation i is strictly directly revealed preferred to the bundle that
% was chosen in observation j. otherwise it equals 0.
SDRP = ceil((REF-index_threshold) / (max(max(abs(REF-index_threshold)))+1));

set_matlab_bgl_default(struct('full2sparse',1));

% A matrix that in the i'th row and j'th column has a positive number represents the shortest 
%path from vertex i to vertex j, and inf if there is no path between i and j.
NS_RP = all_shortest_paths(DRP);

%RP represents the binary relation R.
RP = zeros(rows);

for j=1:rows
    for k=1:rows
        if ~isinf(NS_RP(j,k))
            RP(j,k) = 1;
        end
    end
end

% GARP has 1 in the i'th row and the j'th column iff (x^iRx^j and x^j P^0 x^i), 
% otherwise, it equals 0. 
GARP = RP .* (SDRP');

% rows_comp is the number of GARP cycles without subcycles. 
[rows_comp,~] = size(Components);

% the number of potential observations for shifting. 
counter_adj = 0;

% potential_adj is a matrix with two columns; at every row of the first (counter_adj)
% rows, the first column is a potential observation for shifting, and the second is 
% the ratio of shifting for that potential observation. The following are rows of zeros. 

% the next definition should be removed:
% potential_adj = zeros (nchoosek(rows,2),2);

% instead:
s = 0;
for i=1:rows_comp
        s = s + sum(Components(i,:) > 0);
end
potential_adj = zeros (s,2);

% goes through every cycle num_elements is the number of observations in cycle i.
for i=1:rows_comp 
    
    num_elements = sum(Components(i,:) > 0);
   
    % pairs is a list of pairs obtained from the observations in cycle i.
    % Note: necessarily wRz and zRw for every pair in the cycle.
    pairs = nchoosek(Components(i,1:num_elements),2);
    
    % goes through every such pair
    for j=1:nchoosek(num_elements,2)
        % sets the matrix GARP so it won't see the pair of observations on row j as violation to GARP. 
       GARP(pairs(j,1),pairs(j,2)) = 0;
       
       GARP(pairs(j,2),pairs(j,1)) = 0;
      % if the observations of pair in row j are not identical:  
       if identical_choice(pairs(j,2),pairs(j,1)) == 0
            % if observation pairs(j,1) is directly strictly revealed preferred to observation pairs(j,2), 
            %( Such that, if we denote them w,z respectively, then(p^w)xz<(p^w)x(w)x(1-(10*index_threshold))  ).
            if RATIO(pairs(j,1),pairs(j,2)) < 1 - (10*index_threshold)
                
                counter_adj = counter_adj + 1;
                
                % observation pairs(j,1) is a potential vertex for parallel shifting of RATIO(pairs(j,1),pairs(j,2)).              
                potential_adj(counter_adj,1) = pairs(j,1);
                
                potential_adj(counter_adj,2) = RATIO(pairs(j,1),pairs(j,2));
                
            end
            
            % if observation pairs(j,2) is directly strictly revealed preferred to observation pairs(j,1), Such 
            % that, if we denote them w,z respectively, then (p^z)x(w)<(p^z)x(z)x(1-(10*index_threshold)) 
            if RATIO(pairs(j,2),pairs(j,1)) < 1 - (10*index_threshold)
                
                counter_adj = counter_adj + 1;
                
                % observation pairs(j,2) is a potential vertex for parallel shifting in ratio of
                % ( 1 - RATIO(pairs(j,2),pairs(j,1)) ).
                potential_adj(counter_adj,1) = pairs(j,2);
                
                potential_adj(counter_adj,2) = RATIO(pairs(j,2),pairs(j,1));
                
            end
          
       end
        
    end   % of pairs in the cycle ("for j" loop)
    
end    % of cycles ("for i" loop) 

% The following is unnecessary we ever get inside the 'if' as it's supposed to be equal to 0
% anyhow. I also checked it on the choi et al. data - it never gets in the 'if'
% if sum(sum(GARP)) > 0
   
    %[GARP_row,GARP_col,~] = find(GARP);
    
    % for i=1:length(GARP_row)
        
       %if identical_choice(GARP_row(i),GARP_col(i)) == 0
           
           %if RATIO(GARP_row(i),GARP_col(i)) < 1 - (10*index_threshold)
           
               % counter_adj = counter_adj + 1;
           
                %potential_adj(counter_adj,1) = GARP_row(i);
          
         %       potential_adj(counter_adj,2) = RATIO(GARP_row(i),GARP_col(i));
                
          % end
           
        %   if RATIO(GARP_col(i),GARP_row(i)) < 1 - (10*index_threshold)
          
            %    counter_adj = counter_adj + 1;
           
               % potential_adj(counter_adj,1) = GARP_col(i);
          
               % potential_adj(counter_adj,2) = RATIO(GARP_col(i),GARP_row(i));
                
        %   end
          
      % end
        
   % end    
    
%end


% if the number of potential budget line adjustments for shifting
% is smaller or equal to 'max_exact', set due to matlab's memory problem 
if counter_adj <= max_exact

% Construct a matrix that represents all the subsets of possible budget line
% adjustments

% 2^counter_adj is the number of options for subsets of the kind (observation,ratio) out of all potential budget line adjustments.

% check is a matrix of ones of zeros.
% Ignoring the last column, every row of the matrix represents 
% a subset of optional shifting - where 1 in the i'th row and the j'th column means the j'th optional potential budget line adjustment (according to matrix 'potential_adj') is in the i'th subset, and 0 otherwise. 

% Henceforth, the observation in the j'th row of the optional 
% budget line adjustments matrix (potential_adj) will be called 
% j'th potential observation, so the set of those observations 
% will be referred as the potential observations. 

% The last column of the matrix was meant to exclude subsets 
% that necessarily doesn't induce the smallest adjustment (a 
% subset can't induce the smallest shifting if it itself has a %subset that satisfies GARP). For the moment, it has no use in 
% the algorithm.   

check = zeros(2^counter_adj, counter_adj+1, 'uint8');

% building the matrix of all potential budget line adjustments 
% subsets: 
for i=1:2^counter_adj
    
    num = i-1;
   
    for k=1:counter_adj
       
        check(i,(counter_adj+1)-k) = rem(num,2);
        
        num = (num-rem(num,2)) / 2;
        
    end
    
end

% reference is a sorted matrix of the potential observations, 
% first by the first column (the originally observation 
% numbers) and then by the parallel shifting ratio. 
reference = sortrows(potential_adj(1:counter_adj,:));

% Go over all the possible subsets and check if they relax the relations
% enough to satisfy GARP. If so, calculate the corresponding indices.

no_adjustments = (1 - (10*index_threshold)) * ones(rows,1);


% start measuring time, in order to estimate total time the process will take 
tic
format shortg;
start_time = clock;
start_time = fix(start_time);
    
for i=1:(2^counter_adj) 
    
    % measuring the time after 10,000 loops, to estimate total time
    if i == 10001
        total_size_of_loop = 2^counter_adj;
        total_size_of_sampled_loop = i-1;
        loop_multiplier = total_size_of_loop / total_size_of_sampled_loop;
        seconds_left = toc * loop_multiplier;
        if seconds_left > HPZ_Constants.estimated_time_to_print
            total_minutes_left = round(seconds_left / 60);
            expected_end_time = addtodate(datenum(start_time), total_minutes_left, 'minute');
            hours_left = floor(total_minutes_left / 60);
            minutes_left = mod(total_minutes_left, 60);
            fprintf('%s : Varian Index is calculated using type %i. Estimated time is %i hours and %i minutes - expected to finish in %s.\n', datestr(start_time), 1, hours_left, minutes_left, datestr(expected_end_time));
        end
    end
    
    adjustments = no_adjustments;
    
    % the following 'if' is unnecessary since (counter_adj+1) was defined as zero, and wasn't changed by 
    % the algorithm. It was originally meant to cut back subsets that necessarily don't 
    % induce the minimum adjustment (as explained earlier). .However, here we 
    % are going through all off the subsets, which might be less efficient.         
    
    % if check(i,(counter_adj+1)) == 0
       
        for j=1:counter_adj
          
            % if the j'th potential observation is in the i'th subset
            if check(i,j) == 1
                
                % the adjustment for the observation in the j'th optional observation,
                % is set to the respective ratio of shifting. notice it might be that there
                % is more than one j that points to the same observation. however, 
                % as 'reference' is sorted, at the end of the 'for' loop for a particular
                % subset i, the adjustment for each observation will be the biggest one
                % suggested in the particular subset, as it should.  
                adjustments(reference(j,1),1) = reference(j,2) - (10*index_threshold);
                
            end
            
        end        

        % the expenditure matrix for checking GARP_v, were v=ones(rows,1)-adjustments 
        % (meaning, the cost of x^t where the prices are p^t, is lowered in a proportion of 
        % adjustments(t), for every observation t).   
        exp_var = expenditure - diag((diag(expenditure)).*(ones(rows,1)-adjustments));
        
        % building the matrices used to examine GARP_v, similarly to the examination of GARP.   
        REF = diag(exp_var)*ones(rows,1)' - exp_var;

        % notice that DRP here doesn't necessarily represent the entire R_v^0 relation because:
        %   1.DRP(i,i) isnt  necessarily 1 for every  i.
        %   2.DRP(i,j) isnt  necessarily 1 for every i,j satisfying (i\=j)and(identical_choice(i,j)=1). 

        DRP = ceil((REF+index_threshold) / (max(max(abs(REF+index_threshold)))+1));
        % SDRP represents P_v^0
        SDRP = ceil((REF-index_threshold) / (max(max(abs(REF-index_threshold)))+1));

        set_matlab_bgl_default(struct('full2sparse',1));

        NS_RP = all_shortest_paths(DRP);

        % we want every bundle x^i to be revealed preferred to itself.  
        RP = eye(rows,rows);

        for j=1:rows
            
            for k=1:cols
                % again, we want every bundle to be revealed preferred to itself.
                % so i added the next 'if':
                if (identical_choice(j,k) == 1) && ~(k == j)
            		RP(j,k) = 1;
                end

                if ~isinf(NS_RP(j,k)) && ~(k == j)
            
                    RP(j,k) = 1;
            
                end
        
            end
    
        end
        
        % GARP represents GARP_v
        GARP = RP.*(SDRP');

        GARP_ERRORS = sum(sum(GARP));
    
        if GARP_ERRORS == 0
    
            % if GARP_v is satisfied we calculate the indices
    
            one_minus_v = (ones(rows,1)-adjustments);
     
            % if the average parallel shifting of budget lines for subset i is the smallest
            % among the subsets checked this far, so set it as average_var  
            if average_var > mean(one_minus_v) 
         
                average_var = mean(one_minus_v);
                
                % calculation of in-samle residuals
                average_var_residuals = HPZ_Consistency_Indices_In_Sample_Residuals_Calc (one_minus_v, @mean);
                
                % average_var_vec = adjustments;
         
            end

             % if the root square of mean of squared parallel shifting of budget lines for
             % subset i is the smallest among the subsets checked this far, so set it as meanssq_var  
            if meanssq_var > sqrt(meansqr(one_minus_v))
                
                meanssq_var = sqrt(meansqr(one_minus_v));
                
                % calculation of in-samle residuals
                meanssq_var_residuals = HPZ_Consistency_Indices_In_Sample_Residuals_Calc (one_minus_v, @(x) sqrt(meansqr(x)));
                
                %meanssq_var_vec = adjustments;
         
            end     
            
            % if the maximum parallel shifting of the budget lines for subset i is the
            % smallest among the subsets checked this fur, so set it as min_var  
            if min_var > max(one_minus_v)
         
                min_var = max(one_minus_v);
                
                % calculation of in-samle residuals
                min_var_residuals = HPZ_Consistency_Indices_In_Sample_Residuals_Calc (one_minus_v, @max);
                
                % min_var_vec = adjustments;         
         
            end 
        
        end
        
  % end (this was meant for the 'if' I erased).
    
end

VARIAN = [min_var, average_var , meanssq_var];

% Varian index is accurate.
Var_exact = 1;

Var_in_sample_residuals = [min_var_residuals , average_var_residuals , meanssq_var_residuals];

else
    
    % if there are more than max_exact potential budget line adjustments,
    % then run an approximation of Varian index.      
    [VARIAN, type, Var_in_sample_residuals] = HPZ_Varian_efficiency_index_approx (expenditure, identical_choice, index_threshold);
    
    Var_exact = type;
end

end

