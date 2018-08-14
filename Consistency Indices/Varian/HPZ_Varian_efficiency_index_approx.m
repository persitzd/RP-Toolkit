function [VARIAN, type, Var_in_sample_residuals, one_minus_v] = HPZ_Varian_efficiency_index_approx (expenditure, identical_choice, index_threshold)

% this function approximates the Varian inconsistency index when exact
% calculation has failed.
% the function also returns the in-sample residuals of the index.



% initializations
max_exact = 12;

% number of observations
[rows, ~] = size(expenditure);

% initializations
average_var = 1;   
meanssq_var = 1;
min_var = 1;
one_minus_v_avg = ones(rows, 1);
one_minus_v_meanssq = ones(rows, 1);
one_minus_v_min = ones(rows, 1);
    


% Calculating the number of cycles in the data (cycles of size 3 for the
% case of two goods)

% The matrix REF has at the cell in the i'th row and the j'th column, the difference
% between the value of the bundle that was chosen in observation i and the bundle
% that was chosen in observation j given the prices of observation i
REF = diag(expenditure) * ones(rows,1)' - expenditure;

% represents the relation R^0.
DRP = ceil((REF+index_threshold) / (max(max(abs(REF+index_threshold)))+1));
% SDRP = ceil((REF-index_threshold)/(max(max(abs(REF-index_threshold)))+1));

% Cycles=Johnson(SDRP'); % Use's Johnson's algorithm to extract list cycles     
[Cycles, ~] = Johnson(DRP'); % Use's Johnson's algorithm to extract list cycles according to relation R^0     

% Create Dependency matrix where every row represents a cycle of length 3. 
% The number of columns is MaxItem. On row i the values in cell(i,j) and (i,k) equal 1 
% and the rest are zeros iff (x^kR^0x^j)&(x^jR^0x^k) (meaning [j,k,j] is a cycle).
% The dependency matrix includes all cycles of length 3, and some of the cycles of length 2.
DepMatrix_sf_sb = DependencyCFGK(Cycles,rows); 

% cycles_num_sf_sb is the number of cycles of length 3 found in the dependency matrix 'DepMatrix_sf_sb'. 
[cycles_num_sf_sb, cols_depMatrix] = size(DepMatrix_sf_sb);

% Discard self loops:
% (sum(DepMatrix_sf_sb,2) - ones(cycles_num_sf_sb,1)) is a vector of length (cycles_num_sf_sb),
% where the i'th row is 0 if the i'th cycle on 'DepMatrix_sf_sb' is a self loop (meaning a cycle of length 2: 
% [x^k,x^k]), and 1 otherwise
DepMatrix_sb = zeros(sum(sum(DepMatrix_sf_sb,2) - ones(cycles_num_sf_sb,1)),cols_depMatrix);

cycles_num_sb = 0;

for i=1:cycles_num_sf_sb
    % if cycle i is not a self loop, add the cycle to DepMatrix_sb 
    if sum(DepMatrix_sf_sb(i,:)) == 2
        
        cycles_num_sb = cycles_num_sb + 1;
        
        DepMatrix_sb(cycles_num_sb,:) = DepMatrix_sf_sb(i,:);
        
    end
    
end

% Discard identical choices
num_of_identical = sum(sum(identical_choice - eye(rows)))/2 ;

% note that every two identical bundles x^j,x^k are counted as exactly 1 cycle. 
DepMatrix = zeros(cycles_num_sb-num_of_identical,cols_depMatrix);

if num_of_identical == 0

    DepMatrix = DepMatrix_sb;
    
    cycles_num = cycles_num_sb;
    
else 

    % identical_pairs will be a matrix with two columns, where  every row is a pair of 
    % identical bundles (chosen from different budget line.
    identical_pairs = zeros(num_of_identical,2);

    identical_index = 0;
    
    % building the matrix identical_pairs. 
    for i=1:rows
    
        for j=(i+1):rows
        
            if identical_choice(i,j) == 1
           
                identical_index = identical_index + 1;
            
                identical_pairs(identical_index,1) = i;
            
                identical_pairs(identical_index,2) = j;
            
            end

        end 

    end
    
    cycles_num = 0;
    
    % marking the cycles of identical bundles. 
    for i=1:cycles_num_sb
        
        no_copy = 0;
       
        % goes over every identical pairs, to check whether the 
        % i'th cycle is of identical bundles
        for j=1:num_of_identical
            
            % if the i'th cycle is of the j'th identical pair, mark it by setting no_copy to 1. 
            if ((DepMatrix_sb(i,identical_pairs(j,1)) == 1) && (DepMatrix_sb(i,identical_pairs(j,2)) == 1))
               
                no_copy = 1;
                
            end
            
        end
        % if the i'th cycle is not of the j'th identical pair, for every j, then count the i'th cycle as a 'legitimate' cycle of length 3, and add it to DepMatrix matrix.   
        if no_copy == 0
            
            cycles_num = cycles_num + 1;
            
            DepMatrix(cycles_num,:) = DepMatrix_sb(i,:);
        
        end
        
    end
    
end

% if the number of cycles of length 3 is smaller/equal to 'max_exact', then set 'type' 
% to 2, and use the first approximation.
if cycles_num <= max_exact
    
    type = 2;

    % The matrix RATIO has at the cell in the i'th row and the j'th column, the ratio between
    % the value of the bundle that was chosen in observation j given the prices of observation
    % i and the value of bundle i. 	
    RATIO = expenditure ./ (diag(expenditure) * ones(rows,1)');


    % Construct a matrix that represents all the subsets of possible budget line adjustments

    % 2^(2*cycles_num) is the number of subsets of possible budget line adjustments, because 
    % for every cycle we can choose which of the 2 lines of the cycle's bundles to move. the value 
    % in i'th row and the j'th column of 'check' matrix is 1, if the j'th observation (first according to the 
    % order of cycles, determined by 'DepMatrix' and then by the original order of the observations), 
    % and 0 otherwise.      

    % The last column of the matrix was meant to exclude subsets that necessarily doesn't
    % induce the smallest adjustment (a subset can't induce the smallest shifting if it itself 
    % has a subset that satisfies GARP). For the moment, it has no use in the algorithm.

    check = zeros(2^(2*cycles_num), 2*cycles_num+1);

    for i=1:(2^(2*cycles_num))

        num = i-1;

        for k=1:(2*cycles_num)

            check(i,(2*cycles_num+1)-k) = rem(num,2);

            num = (num-rem(num,2)) / 2;

        end

    end

    % Construct reference which is a vector of all relevant adjustments where 
    % in each row the first column is the observation number and in the second
    % column the size of the corresponding adjustment

    % 'total' matrix's rows are as the number of cycles, and its columns are as the
    % number of observations.
    total = zeros(cycles_num,cols_depMatrix);

    % recall: rows=cols.
    counter = ones(1,rows);

    for i=1:cycles_num

        for j=1:rows

            if (DepMatrix(i,j) == 1)

                break

            end

        end

        for k=j+1:rows

            if (DepMatrix(i,k) == 1)

               total(counter(1,k),k) = RATIO(k,j);

               total(counter(1,j),j) = RATIO(j,k);

               counter(1,k) = counter(1,k)+1;

               counter(1,j) = counter(1,j)+1;

            end

        end

    end

    % the maximum times an observation is in a cycle of length 3.  
    max_count = max(max(counter));

    % number of observation that are in at least one cycle
    num_obs = rows - sum(counter(:) == 1); 

    relevant = zeros(max_count,num_obs);

    rel_index = 1;

    for i=1:cols_depMatrix
        % if i is in some cycle         
        if counter(1,i) > 1
            % enter the number of observation to the first row of 'relevant'

           relevant(1,rel_index) = i;

           % for every time that observation i is in a cycle of length 3 (exactly (counter(1,i)-1) times).
           for j=1:(counter(1,i)-1)

               relevant(j+1, rel_index) = total(j,i);

           end

           % rel_index counts the observations that are in some cycle of length 3.
           rel_index = rel_index + 1;

       end

    end

    % sorts every row from the second row and forth in a descending order.
    relevant = [relevant(1,:) ; sort(relevant(2:max_count,:), 1, 'descend')];

    [rel_rows, rel_cols] = size(relevant);

    % 2*cycles_num is the number of optional line adjustments. reference is a sorted matrix
    % of the potential observations, first by the first column - the originally observation 
    % numbers, and then by the second column - the parallel shifting ratio.
    reference = zeros((2*cycles_num), 2);

    % building the matrix 'reference'
    counter_ref = 1;

    for k=1:rel_cols

       for m=2:rel_rows

           if relevant(m,k) > 0

               reference(counter_ref,1) = relevant (1,k);

               reference(counter_ref,2) = relevant (m,k);

               counter_ref = counter_ref + 1;

           end

       end

    end

    % Go over all the possible subsets and check if they relax the relations
    % enough to satisfy GARP. If so, calculate the corresponding indices.

    no_adjustments = ones(rows,1);

    % start measuring time, in order to estimate total time the process will take 
    tic
    format shortg;
    start_time = clock;
    start_time = fix(start_time);
    
    for i=1:(2^(2*cycles_num))
        
        % measuring the time after 10,000 loops, to estimate total time
        if i == 10001
            total_size_of_loop = 2^(2*cycles_num);
            total_size_of_sampled_loop = i-1;
            loop_multiplier = total_size_of_loop / total_size_of_sampled_loop;
            seconds_left = toc * loop_multiplier;
            if seconds_left > HPZ_Constants.estimated_time_to_print
                total_minutes_left = round(seconds_left / 60);
                expected_end_time = addtodate(datenum(start_time), round(seconds_left), 'second');
                hours_left = floor(total_minutes_left / 60);
                minutes_left = mod(total_minutes_left, 60);
                fprintf('%s : Varian Index is calculated using type %i. Estimated time is %i hours and %i minutes - expected to finish in %s.\n', datestr(start_time), 2, hours_left, minutes_left, datestr(expected_end_time));
            end
        end

        adjustments = no_adjustments;
        % the following 'if' is unnecessary since (2*cycles_num +1) was defined as zero, and wasn't changed by 
        % the algorithm. It was originally meant to cut back subsets that necessarily don't 
        % induce the minimum adjustment (as explained earlier).However, here we are going 
        % through all of the subsets, which might be less efficient.   
        % if check(i,(2*cycles_num+1)) == 0
        
            for j=1:(2*cycles_num)

                % if the j'th potential observation is in the i'th subset
                if check(i,j) == 1
                    % the adjustment for the observation in the j'th optional observation, is set to the 
                    % respective ratio of shifting. notice it might be that there is more than one j that 
                    % points to the same observation.however, as 'reference' is sorted, at the end 
                    % of the 'for' loop for a particular subset i, the adjustment for each observation will be 
                    % the biggest one suggested in the particular subset, as it should. 
                    adjustments(reference(j,1),1) = reference(j,2) - (10*index_threshold);

                end

            end        

            % the expenditure matrix for checking GARP_v, were v=ones(rows,1)-adjustments (meaning, the 
            % cost of x^t where the prices are p^t, is lowered in a proportion of adjustments(t), for every 
            % observation t).   
            exp_var = expenditure - diag((diag(expenditure)).*(ones(rows,1)-adjustments));

            
            % GARP_v
            [GARP_v, ~, ~ , ~] = GARP_based_on_expenditures(exp_var, [], index_threshold);
            

            GARP_ERRORS_v = sum(sum(GARP_v));

            if GARP_ERRORS_v == 0

                % if GARP_v is satisfied we calculate the indices

                one_minus_v = (ones(rows,1) - adjustments);

                % if the average parallel shifting of budget lines for subset i is the smallest among 
                % the subsets checked this fur, so set it as average_var
                if average_var > mean(one_minus_v) 

                    one_minus_v_avg = one_minus_v;
                    
                    average_var = mean(one_minus_v);
                    
                    % average_var_vec = adjustments;

                end

                % if the root square of mean of squared parallel  shifting of budget lines for subset i 
                % is the smallest among the subsets checked this fur, so set it as meanssq_var  
                if meanssq_var > (sqrt(meansqr(one_minus_v))) 

                    one_minus_v_meanssq = one_minus_v;
                    
                    meanssq_var = (sqrt(meansqr(one_minus_v)));
                    
                    % meanssq_var_vec = adjustments;

                end     

                % if the maximum parallel shifting of the budget  lines for subset i is the smallest
                % among the subsets checked this fur, so set it as min_var  
                if min_var > max(one_minus_v)

                    one_minus_v_min = one_minus_v;
                    
                    min_var = max(one_minus_v);
                    
                    % min_var_vec = adjustments;         

                end 

            end
        
      %  end (for the erased 'if') 
    
    end
    
    % assigning to the result variables:
    VARIAN = [min_var, average_var , meanssq_var];
    
    % calculation of in-samle residuals
    min_var_residuals = HPZ_Consistency_Indices_In_Sample_Residuals_Calc (one_minus_v_min, @max); 
    average_var_residuals = HPZ_Consistency_Indices_In_Sample_Residuals_Calc (one_minus_v_avg, @mean);
    meanssq_var_residuals = HPZ_Consistency_Indices_In_Sample_Residuals_Calc (one_minus_v_meanssq, @(x) sqrt(meansqr(x)));
    Var_in_sample_residuals = [min_var_residuals , average_var_residuals , meanssq_var_residuals];

    % the best 1-v for each of the aggregators
    one_minus_v = [one_minus_v_min, one_minus_v_avg, one_minus_v_meanssq];
    
    
    
else    % If there are more than 'max_exact' cycles of length 3 then run the  
        % second approximation (type 3)
    
    [VARIAN, type, Var_in_sample_residuals, one_minus_v] = HPZ_Varian_efficiency_index_approx_second (expenditure, identical_choice, index_threshold);

end

end

