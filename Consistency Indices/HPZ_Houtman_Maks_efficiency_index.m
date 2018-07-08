function [HM, HM_exact, HM_residuals] = HPZ_Houtman_Maks_efficiency_index (Choices, index_threshold, residuals_flag)

% this function calculates the Houtman-Maks inconsistency index
% if the calculation seems to be taking too long, it will return an
% approximation of the Houtman-Maks inconsistency index. HM_exact is a flag 
% that tells whether it was an exact calculation (1) or not (0)
% either way, the function also returns the residuals of the index, if specified. 



% Choices is a mX4 matrix where the rows are observations and the columns
% are X Y PX PY.



% Boodaghians and Vetta (2015) prove that for two goods, the computation of
% the HM index can be found in polynomial time, while for three goods or
% more it is an NP-Hard problem.

% The proof of Boodaghians and Vetta (2015) for the case of two goods 
% contains three steps.

% First, they recall Rose (1958) that shows that every cycle of
% observations that violates SARP, must contain a pair of observations that
% violates WARP. Therefore they deduce that the HM index is the minimal
% number of observations needed to be dropped for all the size 2 cycles to
% be broken.

% Second, to do that they construct a graph where the nodes are observations and
% two nodes are linked if and only if the corresponding observations
% construct a size 2 cycle. Then they show that this graph is a member of a
% set of graphs called "perfect graphs".

% Third, they point out that calculating the HM index on the graph described 
% above is identical to a problem in computer science called the vertex cover 
% problem (finding the minimal subset of nodes such that every edge has at 
% least one end in this subset). While generally the Vertex cover problem
% is NP-HARD, Grotschel et al (1984, 1988) showed that if the graph is
% perfect then there is a polynomial time algorithm that solves it exactly. 
% Therefore, they conclude that the computation of the HM index in the case 
% of two goods can be found in polynomial time 

% The problem is that currently I was not able to find an implementation of 
% this algorithm and therefore I was not able to exploit this finding for our 
% needs. This is a work for the near future. 

% However, I did made the adjustment of this code to the two goods case,
% usually used in lab experiments. Thus, the current code implements a two 
% stage algorithm. 

% First, an exhaustive search is implemented, where a huge number of subsets   
% is checked directly (GARP is tested on a data set that does not include
% the observations in these subsets). For example if the data set includes
% 50 observations then all subsets of size 46 or more are tested directly 
% (provided that max_subsets = 1000000).

% If no subset that satisfies GARP is found, a complementary toolbox that
% contains a general code for the linear programing solution of the HM index 
% problem is adjusted for the case of two goods. The adjustment turns the 
% linear program much simpler and therefore increases its accuracy.   

% To avoid confusion the function now returns, on top of the index, a flag 
% that reports whether the index is accurate or an upper bound approximation. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First we go with the exhaustive search.
% number_trials is the number of trials in the largest subset we'll check.
% subsets is the number of subsets that we check directly.

max_subsets = 10000000;

% trials is set to be the number of rows of matrix Choices, meaning m, the 
% number of observations for a single subject.

[trials,~] = size(Choices);

subsets = 0;

% loop for i decreasing by 1 from (m-1) to 1, going through all of the subsets at size i, in 
% order to find the number of observations which satisfies that:
% a. the number of subsets of size strictly higher than i and smaller or equal to (m-1), is not 
% more than the max_subsets.  
% b. the number of subsets of size higher or equal to i and smaller or equal to (m-1), is 
% more than the max_subsets.  

for i=(trials-1):-1:1
    
    % trials choose i is the number of options for subsets of size i chosen from a set with 'trial' members.
    subsets_temp = subsets + nchoosek(trials,i);
    
    % if the sum is overstepping the defined maximum number of subsets - max_subsets, then the loop breaks  
    if subsets_temp > max_subsets
        break
    end
    
    subsets = subsets_temp;

    % number_trials is the the number of observations which satisfies the conditions described above.
    number_trials = i;
    
end

% subsets_matrix holds in each row one subset of the require size

% vector_of_trials=[1,2,3,4,…,trials]
vector_of_trials = 1:1:trials;

% start measuring time, in order to estimate total time the process will take 
tic
format shortg;
start_time = clock;
start_time = fix(start_time);

% the following loop over the number of trials will continue until it is 
% finished, or until this flag will be turned to 1 (which will happen when
% we find a truncated data that satisfies GARP)
HM_found_flag = 0;

% initialization of residuals
HM_residuals = zeros(trials, 1);

% loop for j decreasing by 1 from (m-1) to number_trials.  
% Every time it will "deal" with subsets of the same size (j), and look for a subset that satisfies GARP
j = trials-1;
while (j >= number_trials) && (~HM_found_flag)
    
    % measuring the time after the first loop, to estimate total time
    if j == (trials-2)
        total_size_of_loop = subsets;
        total_size_of_sampled_loop = nchoosek(trials,1);
        loop_multiplier = total_size_of_loop / total_size_of_sampled_loop;
        seconds_left = toc * loop_multiplier;
        if seconds_left > HPZ_Constants.estimated_time_to_print
            total_minutes_left = round(seconds_left / 60);
            expected_end_time = addtodate(datenum(start_time), total_minutes_left, 'minute');
            hours_left = floor(total_minutes_left / 60);
            minutes_left = mod(total_minutes_left, 60);
            fprintf('%s : HOUTMAN-MAKS Index is calculated using type %i. Maximal Estimated time is %i hours and %i minutes - expected to finish up to %s.\n', datestr(start_time), 1, hours_left, minutes_left, datestr(expected_end_time));
        end
    end
    
    % every row in the subsets_matrix is an option for a subset of observations numbers with j members, 
    % from the set of all observations numbers (for example: {1,2,3,4,…j} or {2,3,4,5,…,j+1}). 
    % The matrix includes all of the options for sets of observations numbers of size j.    
    % matrix size: (vector_of_trials choose j) x j. 
    % note: 1. Every subset is ordered in ascending order. 
    %       2. setting the options in the matrix induces a certain order of the options.  
    
    subsets_matrix = nchoosek(vector_of_trials,j);
    
    [rows,~] = size(subsets_matrix);
    
    % The loop goes through every option of a subset of size j and checks if 
    % it satisfies GARP. We will be using similar calculations we made when 
    % we checked if the entire data satisfies GARP, with several adjustments. 
    % Henceforth we will refer to the k'th row in subsets_matrix, 
    % meaning the k'th option for a subset of observation numbers of size j, as option k for j.
    
    for k=1:rows
        
        % new_choices is a j x 4.
        % Row i is the row of the i'th observation number of k option for j, 
        % from the matrix Choices (if for example: i=2,j=7 and the k option 
        % for 7 is: {2,3,4,5,6,7,8} the i'th observation number of k option 
        % for j is observation 3. row 2 in the matrix new_choices is the 
        % third row in the Choices matrix). 
        % The new_choices matrix is a "partial" matrix of Choices, 
        % as it includes only the rows suited for the observations that 
        % belong to the subset.   
        
        new_choices = Choices(subsets_matrix(k,:),:); 
        
        % we will henceforth refer to the bundle of the observation in the i'th 
        % row of "new_choices" matrix (meaning: (new_choices(i,1),new_choices(i,2)) ), 
        % as observation i*, and the respective prices (meaning: (new_choices(i,3),new_choices(i,4)) ) 
        % as i* prices.
        
        % The matrix "expenditure" has at the cell in the i'th row and the j'th column, the value of 
        % the observation j* given the prices of observation i*
        
        new_expenditure = (new_choices(:,1)*new_choices(:,3)' + new_choices(:,2)*new_choices(:,4)')';
        
        % The matrix REF has at the cell in the i'th row and the j'th column, 
        % the difference between the value of the bundle of observation i* and 
        % the bundle of observation j* given the prices of observation i*
        
        new_REF = diag(new_expenditure) * ones(j,1)' - new_expenditure;
        
        % new_DRP(i,j)=1 iff the bundle i* is directly revealed preferred to 
        % the bundle j (for a index_threshold small enough). See detailed explanation
        % in HPZ_Subject_Consistency
        
        new_DRP = ceil((new_REF+index_threshold) / (max(max(abs(new_REF+index_threshold)))+1));
        
        % The matrix new_SDRP has at the cell in the i'th row and the j'th column, 1 if and only if the 
        % bundle of observation i* is strictly directly revealed preferred to the bundle of observation 
        % j* , and 0 otherwise.
        
        new_SDRP = ceil((new_REF-index_threshold) / (max(max(abs(new_REF-index_threshold)))+1));
        
        % statement needed for the graph theory external package to work efficiently
        
        set_matlab_bgl_default(struct('full2sparse',1));
        
        % The matrix new_NS_RP has at the cell in the i'th row and the j'th column, Inf if and only if the 
        % bundle of observation i* is not revealed preferred to the bundle of observation j*. 
        % otherwise it includes a positive integer.
        
        new_NS_RP = all_shortest_paths(new_DRP);
        
        % The matrix new_RP has at the cell in the i'th row and the j'th column, 
        % 1 if and only if the bundle i* is revealed preferred to the bundle j*. 
        % Otherwise, it equals 0.      
        
        new_RP = zeros(j);
        
        for l=1:j

            % going through all new_NS_RP’s cells.
            
            for m=1:j

                if ~isinf(new_NS_RP(l,m))

                    % if the length of the path from j to k is 
                    % finite, meaning there is a path from j to k (or in 
                    % other words bundle j* is revealed preferred to bundle
                    % k*), then RP(j,k)=1, otherwise RP(j,k) stays 0. 
                    
                    new_RP(l,m) = 1;

                end
              
            end
          
        end

        % The matrix is the zero matrix if and only if GARP is satisfied.
      
        GARP = new_RP .* (new_SDRP');

        % A flag that keeps the GARP result. It is 0 if and only if the data satisfies GARP.
      
        GARP_FLAG = 1;

        % GARP_ERRORS sums the values of all matrix GARP's cells. Counting 
        % the number of violations (notice that a pair of bundles x^i*,x^j* might be counted 
        % as two violations; one for (i,j), and the other for (j,i)).

        GARP_ERRORS = sum(sum(GARP));

        % If the GARP matrix is only zeros then the data satisfies GARP      
      
        if GARP_ERRORS == 0

            GARP_FLAG = 0;

        end 
      
        if GARP_FLAG == 0
            
            % We found the highest number of observations that satisfies GARP.
            % the HM index is calculated as the number of removed relatively to the number of 
            % observations.
            
            % HM index
            HM = (trials - j) / trials;
            
            % the calculation of the HM index is exact. 
            HM_exact = 1;
            
            if residuals_flag
                % we don't need to check smaller subsets of the data
                HM_found_flag = 1;
                % HM residuals per observation
                % it is enough that the observation will NOT appear in one
                % subset (that is - will be truncated in one subset)
                % , for it to have a residual of 1 instead of 0
                observations_subset = zeros(1,trials);
                observations_subset(subsets_matrix(k,:)) = 1;
                HM_residuals_indexes = ~observations_subset;
                HM_residuals(HM_residuals_indexes) = 1;
            else
                % no need to calculate residuals -
                % the loops are stopped. We return HM and exact=1.           
                return
            end
            
          
        end
        
    end
    
    % decrement the counter
    j = j - 1;
end

% we return the j to its value in the last loop that was actually performed 
j = j + 1;
% we modify the raw residuals ( 0 or 1 ) 
% to actual residuals ( (trials-j)/(trials-1) or (trials-(j+1))/(trials-1) )   
HM_residuals(:) = (trials - (j + HM_residuals(:))) / (trials-1);


% If we were not able to directly compute the index, we have to approximate.
% And this is done here using HoutmanMaksCFGK.m that was written by Daniel Martin. 
% The package of Matlab files that calculates the Houtman Maks index was 
% downloaded from Daniel Martin's personal website on November 5th 2011.
% For the implementation I used the readme file that is located next to the
% package on the same website.

% Comment 1: Some adaptation was introduced for the case of two goods

% Comment 2: the approximation written by Daniel Martin works with a single binary relation. However, 
% testing for GARP requires two relations R and P^0. As P^0 and R^0 are often identical, it almost 
% always doesn't matter. However, it is possible of course that xR^0y but not xP^0y, and so 
% we ought to choose a relation for the approximation algorithm. 
% As the approximation algorithm is as an upper bound, we chose "working" with the relation R. 
% The consequence of that, is that we might add cycles that do not exist. 
% While choosing the relation P^0, on the other hand, might cause at missing cycles 
% (those satisfies xRy, yP^0x, but not xPy). So using  R^0 for the algorithm gives 
% us an upper bound of the number of observations one must remove to eliminate all GARP cycles .  

if (GARP_FLAG) && (~HM_found_flag)
    [HM, HM_exact] = HPZ_Houtman_Maks_efficiency_index_approx (Choices, index_threshold);
end

end

