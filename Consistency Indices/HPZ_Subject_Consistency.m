function [expenditure, DRP, SDRP, RP, SRP, WARP, GARP, SARP, FLAGS, VIO_PAIRS, VIOLATIONS, AFRIAT, VARIAN, Var_exact, HM, HM_exact, MPI, Mat] = HPZ_Subject_Consistency (data_matrix, GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags, active_waitbar, current_run, total_runs)

% this function calculates all the currently implemented consistency and
% inconsistency indices, and their residuals (and depending on the flags 
% it may calculate only some of them), for a single subject.
% the WARP, GARP and SARP measures are implemented within the function
% itself, while the inconsistency indices (Afriat, Varian, Houtman-Maks) are
% being calculated by calling other specific functions.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 

%% output variables :
% expenditure,DRP,SDRP,RP,SRP,WARP,GARP,SARP,FLAGS :
%   all these are not really needed, but are the result of intermediate
%   calculations, and are returned only because someone might want to use
%   them for some purpose.
% VIO_PAIRS,VIOLATIONS :
%   are each 3-length vectors, when the 1st element is for WARP, the 2nd is
%   for GARP, and the 3rd is for SARP.
%   they are measures of inconsistency with WARP/GARP/SARP, by counting
%   number of violations or number of pairs of choices involved in
%   violations.
% AFRIAT,VARIAN,Var_exact,HM,HM_exact :
%   inconsistency indices according to each of these methods.
%   VARIAN is a 3-length vector, since 3 aggregators are implemented:
%   (1) Min, (2) Mean, (3) Average Sum of Squares.
%   Var_exact and HM_exact tell the user whether the calculation was exact
%   or was it only an approximation, and what type of approximation.
% Mat :
%   a matrix with all the residuals results as should be printed to the
%   results file





% are we running on data where the choice of the DM is between only 2
% goods, and not 3 or more (we need it because when it is only 2 goods, 
% we can improve the speed in Houtman-Maks)
% in the current version of the code it is always true cause we haven't yet
% implemented a general framework for more than 2 goods, but it will be
% relevant in the next version (version 3)
is_2_goods = true;



%% Subject_Consistency

% Variables:
% ID - the subjects ID.
% Choices - a mX4 matrix where the rows are observations and the columns
% are X Y PX PY.

% Let A be an amount of tokens. Let B also be an amount of tokens. 
% Then:
%   1. We will consider A > B if and only if A-B > index_threshold.
%      This will be used when we calculate SDRP (strict direct revealed preference).
%   2. We will consider A >= B if and only if A-B >= -index_threshold.
%      This will be used when we calculate DRP (direct revealed preference).
%   In other words, we take a conservative approach to SDRP, and on the
%   other hand a liberal approach to DRP. For the reason we do so, and for
%   how we determined the threshold, see the documentation in "HPZ_Constants".  
% we will use the threshold for numerical reasons only.
index_threshold = HPZ_Constants.index_threshold;

% this threshold is also needed for numerical reasons only. it is used in
% determining whether two bundles are identical or not.
% for further information, see the documentation in "HPZ_Constants".
identical_choices_threshold = HPZ_Constants.identical_choices_threshold;



% number of observations
[obs_num , ~] = size(data_matrix);

% these are for convenience
basic_flags = [GARP_flags(1) , AFRIAT_flags(1) , VARIAN_flags(1) , HOUTMAN_flags(1) , MPI_flags(1)];
residuals_flags_temp = [GARP_flags(2) , AFRIAT_flags(2) , VARIAN_flags(2) , HOUTMAN_flags(2) , MPI_flags(2)];
residuals_flags = min(basic_flags , residuals_flags_temp);
in_sample_flags = [GARP_flags(3) , AFRIAT_flags(3) , VARIAN_flags(3) , HOUTMAN_flags(3) , MPI_flags(3)];
out_sample_flags = [GARP_flags(4) , AFRIAT_flags(4) , VARIAN_flags(4) , HOUTMAN_flags(4) , MPI_flags(4)];





% define the waitbar
if (active_waitbar)
    indices_strings = {' GARP' , ' Afriat', ' Varian', ' Houtman-Maks', ' Money-Pump'};
    for i=1:length(basic_flags)
        if basic_flags(i) == 0
            indices_strings{i} = '';
        end
    end
    waitbar_name = char(strcat(HPZ_Constants.waitbar_name_calculation, {' '}, '(', HPZ_Constants.current_run_waitbar, {' '}, num2str(current_run), {' '}, HPZ_Constants.total_runs_waitbar, {' '}, num2str(total_runs), ')'));
    waitbar_msg = char(strcat(HPZ_Constants.waitbar_calculation, {' '}, num2str(data_matrix(1,1)), {' '}, HPZ_Constants.waitbar_indices));
    new_bar_val = 0;
    h_wb = wide_waitbar(new_bar_val, {waitbar_msg, ''}, waitbar_name, HPZ_Constants.waitbar_width_multiplier);
end





%% initialization of residuals result Matrix

if sum(residuals_flags) == 0
    % then no residuals are required
    Mat = 0;
else
    % the residuals are required
    
    % determining the size of the result matrix:
    
    % number of columns
    % initialization to 2, since we always print the subject number
    % and the number of the observation in the first and second columns
    num_of_columns = 2;
    % for each of the 4 measures (GARP VIO PAIRS, AFRIAT, VARIAN, HOUTMAN-MAKS)
    % we check whether we need to calculate its residuals and whether to
    % calculate them in sample or out of sample or both
    % (the number of columns will be between 1 and 9)
    for i=1:5
        if (residuals_flags(i) == 1)
            if (in_sample_flags(i) == 1)
                if (i == 3)
                    % there are 3 in sample measures for VARIAN -
                    % one for MIN, one for MEAN and one for AVGSSQ
                    num_of_columns = num_of_columns + 3;
                elseif (i == 1)
                    % we have WARP, GARP and SARP, 
                    % for each of them we have 2 columns 
                    num_of_columns = num_of_columns + 2*sum(GARP_flags(5:7));
                else
                    num_of_columns = num_of_columns + 1;
                end   
            end
            if (out_sample_flags(i) == 1)
                % in out-of-sample we print 3 things:
                %   (1) the full index (when all observation are taken into account)
                %   (2) the partial index (when omitting one observation)
                %   (3) the difference between (1) and (2)
                %   in VARIAN and HOUMAN-MAKS we also have:
                %   (4) the normalized difference between (1) and (2)
                %   also, in VARIAN and HOUMAN-MAKS we also have to print
                %   whether the original index and whether the residual
                %   index were exact or approximated
                if (i == 1)
                    % in violations counting we have WARP, GARP and SARP
                    num_of_columns = num_of_columns + 3*sum(GARP_flags(5:7));
                elseif (i == 3)
                    % there are 3  measures for VARIAN -
                    % one for MIN, one for MEAN and one for AVGSSQ
                    % for each of these 3 we need: full index, out-of-sample
                    % index, difference and normalized difference.
                    % we also need to print is_exact both for the full data
                    % and for the truncated data.
                    num_of_columns = num_of_columns + 3*4+2;
                elseif (i == 4)
                    % in HOUTMAN-MAKS we have the full index and whether it
                    % is exact, the residual index and whether it is exact,
                    % and difference and normalized difference.
                    num_of_columns = num_of_columns + 4+2;
                elseif (i == 5)
                    % in MPI we have 2 aggregators: mean and median
                    num_of_columns = num_of_columns + 3*2;
                else
                    % full index, out-of-sample index, difference
                    % (this is default, but currently implied only for Afriat) 
                    num_of_columns = num_of_columns + 3;
                end
            end
        end
    end
    
    % creating the results matrix for the residuals
    Mat = zeros (obs_num, num_of_columns);
    % entering the subject ID and the observations numbers
    for i=1:obs_num
        Mat(i, 1) = data_matrix(1,1);
        Mat(i, 2) = i;
    end 
    
    % this counter tells us which column we should print to next
    % initialized to 3 (start with the 3rd column)
    col_counter = 3;
    
    % this is a preparation for out of sample residuals - 
    % creating truncated data matrices:
    truncated_datas = zeros(obs_num , (obs_num-1) , 6);
    truncated_datas(1,:,:) = data_matrix(2:obs_num,:);
    truncated_datas(obs_num,:,:) = data_matrix(1:(obs_num-1),:);
    % we assume there are at least 3 observations:
    for i=2:(obs_num-1)
        truncated_datas(i,:,:) = vertcat(data_matrix(1:(i-1),:) , data_matrix((i+1):obs_num,:));
    end
end





% Choices
Choices(:,1:4) = data_matrix(1:obs_num,3:6);   % quantities & prices

% The matrix "expenditure" has at the cell in the i'th row and the j'th
% column, the value of the bundle that was chosen in observation j given the
% prices of observation i
expenditure = (Choices(:,1)*Choices(:,3)' + Choices(:,2)*Choices(:,4)')';

% record the dimensions of expenditure
[rows , cols] = size(expenditure); %#ok<ASGLU>

% the purpose of the matrix identical_choice is to locate identical choices. Value of 
% cell (j,k) is 1 if choices identical, and 0 otherwise. The following loop builds the 
% matrix: 
identical_choice = zeros(rows);

for j=1:rows
    % going through all identical_choice�s cells.
    for k=1:rows
        %if the choices j&k are identical, then set value of cell (j,k) = 1 , otherwise 0. 
        if (abs(Choices(j,1) - Choices(k,1)) <= Choices(j,1)*identical_choices_threshold) && (abs(Choices(j,2) - Choices(k,2)) <= Choices(j,2)*identical_choices_threshold)
              identical_choice(j,k) = 1;
        end
    end
end



% The matrix REF has at the cell in the i'th row and the j'th
% column, the difference between the value of the bundle that was chosen in 
% observation i and the bundle that was chosen in observation j given the 
% prices of observation i

REF = diag(expenditure)*ones(rows,1)' - expenditure;

% DRP is a representation of the relation R^0.
% The matrix DRP has at the cell in the i'th row and the j'th
% column, 1 if and only if the bundle that was chosen in 
% observation i is directly revealed preferred to the bundle that was chosen 
% in observation j. otherwise it equals 0.

% DRP(i,j) = 1 iff REF_DRP(i,j) > 0 iff the bundle that was chosen
% in observation i is directly revealed preferred to the bundle that was chosen 
% in observation j (for a index_threshold small enough)

% we increase REF by a small threshold in order to make sure it is bigger
% than 0, also we add to it "identical_choice" in order to make sure that
% indentical choices will be considered DRP to one another (which is
% crucial for WARP and SARP, otherwise we might have negative numbers in
% WARP and SARP).

DRP = (REF + identical_choice + diag(expenditure)*ones(rows,1)'*index_threshold >= 0) * 1;

% SDRP is a representation of the relation P^0.
% The matrix SDRP has at the cell in the i'th row and the j'th
% column, 1 if and only if the bundle that was chosen in 
% observation i is strictly directly revealed preferred to the bundle that was chosen 
% in observation j, and 0 otherwise.

% we decrease REF by a small threshold in order to make sure that if it is 
% bigger than 0, it is because it should be, and not as a result of a
% computational inprecision.

SDRP = (REF - diag(expenditure)*ones(rows,1)'*index_threshold > 0) * 1;



% statement needed for the graph theory external package to work
% efficiently
set_matlab_bgl_default(struct('full2sparse',1));

% The matrix NS_RP has at the cell in the i'th row and the j'th
% column Inf if and only if the bundle that was chosen in 
% observation i is not revealed preferred to the bundle that was chosen 
% in observation j. Otherwise it includes a positive integer.
NS_RP = all_shortest_paths(DRP);

% The matrix RP has at the cell in the i'th row and the j'th
% column 1 if and only if the bundle that was chosen in 
% observation i is revealed preferred to the bundle that was chosen 
% in observation j. Otherwise, it equals 0. 

% NEW CODE:
RP = ~isinf(NS_RP);

% OLD CODE:
% RP = zeros(rows);
% 
% for j=1:rows
%     % going through all NS_RP�s cells.
%     for k=1:cols
%         % if the length of the path from j to k is 
%         % finite, meaning there is a path from j to k (or in 
%         % other words bundle j is revealed preferred to bundle
%         % k), then RP(j,k)=1, otherwise RP(j,k) stays 0. 
%         if ~isinf(NS_RP(j,k))   
%             RP(j,k) = 1;
%         end
%     end
% end



% The matrix NS_SRP has at the cell in the i'th row and the j'th
% column, zero if and only if the bundle that was chosen in 
% observation i is not strictly revealed preferred to the bundle that was chosen 
% in observation j. Otherwise it includes a positive integer.

% notice that x^iPx^j iff there exist bundles x^s,x^t such that:
% (x^i R x^s) & (x^s p0 x^t) & (x^t R x^j).
% Now, calculating the (i,j) cell of the matrix shows that:   
% NS_SRP(i,j)=sum[t=1 to m] ( sum[s=1 to m] PR(i,s)*SDRP(s,t) )*RP(t,j)
% NS_SRP(i,j)=\ 0 iff there exist bundles x^s,x^t such that:
%(x^i R x^s) & (x^s p0 x^t) & (x^t R x^j). 
NS_SRP = (RP*SDRP)*RP;
% % An alternative, maybe better for our needs SRP (not implemented):
% RP_no_self = RP;
% for i=1:obs_num
%     RP_no_self(i,i) = 0;
% end
% NS_SRP = (RP_no_self*SDRP)*RP_no_self;



% The matrix SRP has at the cell in the i'th row and the j'th
% column, 1 if and only if the bundle that was chosen in 
% observation i is strictly revealed preferred to the bundle that was chosen 
% in observation j. 

% NEW CODE:
SRP = (NS_SRP ~= 0);

% OLD CODE:
% SRP = zeros(rows);
% 
% for j=1:rows
%     % going through all NS_SRP�s cells.
%     for k=1:cols
%         if ~(NS_SRP(j,k) == 0) 
%             % if x^j P x^k then sets SRP(j,k)=1 otherwise, it stays 0.
%             SRP(j,k) = 1;
%         end
%     end
% end



% To test for SARP we will use the definition of SARP1. For every pair of
% choices x and y if xRy and yRx it must be that x=y. We will take RP and
% its transpose and multiply element by element. Every 1 corresponds to 
% xRy and yRx. Then we take off the identity matrix (all the pairs x=y). 
% The final matrix is the zero matrix if and only if SARP is satisfied. 


% The following conforms with Varian (1982) definition of SARP1 and SARP2
% SARP3 would be SARP = RP.*(DRP') - identical_choice;
% As a binary test these definitions are equivalent. However, SARP1 and
% SARP2 may report more violations.

% RP(i,j)*RP'(j,i)=1 iff (x^i R x^j)&(x^j R x^i), otherwise it equals 0. For every i,j identical
% choice SARP(i,j)=0 (including i=j of course) 

SARP = RP.*(RP') - identical_choice;

% A flag that keeps the SARP result. It is 0 if and only if the data
% satisfies SARP.
SARP_FLAG = 1;

SARP_ERRORS = sum(sum(SARP));   % sums the values of all matrix SARP's cells, counting 
                                % the number of violations (notice that every pair of 
                                % bundles x^i,x^j is counted as two violations; one for 
                                % (i,j), and the other for (j,i)). 


% sums the values of all matrix SARP's cells, counting the number of violations 
% pairs (every two pair of bundles x^i,x^j is counted as one violation pair). 
SARP_VIO_PAIRS = sum(sum(triu(SARP|(SARP'))));

% If the SARP matrix is only zeros then the data satisfies SARP
if SARP_ERRORS == 0
  SARP_FLAG = 0;
end 



% To test for GARP we will do the following: for every pair of
% choices x and y if xRy then not yP0x. We will take RP and the transpose of
% SDRP and multiply element by element. Every 1 corresponds to 
% xRy and yP0x. The final matrix is the zero matrix if and only if 
% GARP is satisfied. 
GARP = RP .* (SDRP');

% A flag that keeps the GARP result. It is 0 if and only if the data
% satisfies GARP.
GARP_FLAG = 1;

% GARP_ERRORS sums the values of all matrix GARP's cells. Counting 
% the number of violations (notice that a pair of bundles x^i,x^j might be counted 
% as two violations; one for (i,j), and the other for (j,i)). 
GARP_ERRORS = sum(sum(GARP));

% GARP_VIO_PAIRS sums the values of all matrix
% GARP's cells, counting the number of violations pairs (every two pair of bundles 
% x^i,x^j is counted as one violation pair). 
GARP_VIO_PAIRS = sum(sum(triu(GARP|(GARP'))));

% If the GARP matrix is only zeros then the data satisfies GARP
if GARP_ERRORS == 0
  GARP_FLAG = 0;
end 



% To test for WARP we will do the following: for every pair of
% choices x and y if xR0y then not yR0x. We will take DRP and the transpose of
% DRP and multiply element by element. Every 1 corresponds to 
% xR0y and yR0x. BUT we fix it so identical choices will not be counted, as
% they do not violate WARP.
% The final matrix is the zero matrix if and only if 
% WARP is satisfied. 

WARP = DRP.*(DRP') - identical_choice;

% A flag that keeps the WARP result. It is 0 if and only if the data
% satisfies WARP.
WARP_FLAG = 1;

% WARP_ERRORS sums the values of all matrix WARP's cells. Counting 
% the number of violations (notice that a pair of bundles x^i,x^j might be counted 
% as two violations; one for (i,j), and the other for (j,i)). 
WARP_ERRORS = sum(sum(WARP));

% WARP_VIO_PAIRS sums the values of all matrix
% GARP's cells, counting the number of violations pairs (every two pair of bundles 
% x^i,x^j is counted as one violation pair). 
WARP_VIO_PAIRS = sum(sum(triu(WARP|(WARP'))));

% If the WARP matrix is only zeros then the data satisfies WARP
if WARP_ERRORS == 0
  WARP_FLAG = 0;
end 



% assigning the results to the results vectors
FLAGS(1) = WARP_FLAG;
FLAGS(2) = GARP_FLAG;
FLAGS(3) = SARP_FLAG;

VIO_PAIRS(1) = WARP_VIO_PAIRS;
VIO_PAIRS(2) = GARP_VIO_PAIRS;
VIO_PAIRS(3) = SARP_VIO_PAIRS;

VIOLATIONS(1) = WARP_ERRORS;
VIOLATIONS(2) = GARP_ERRORS;
VIOLATIONS(3) = SARP_ERRORS;





%% residuals for WARP/GARP/SARP VIO PAIRS

% the full number of GARP VIO PAIRS violations
WARP_full = WARP_VIO_PAIRS;
GARP_full = GARP_VIO_PAIRS;
SARP_full = SARP_VIO_PAIRS;

if (sum(residuals_flags .* out_sample_flags) > 0)
    % we need to perform out-of-sample for GARP for 2 possible reasons:
    % (1) the user asked for WARP/GARP/SARP out-of-sample residuals
    % (2) the user asked for AFRIAT/VARIAN/HM/etc. out-of-sample
    %     residuals, and we would like to confirm first whether the 
    %     truncated data is consistenct with GARP, because if it is
    %     then there is no need to calculate the indices themselves.

    % this matrix is a helper matrix for WARP/GARP/SARP out-of-sample residuals. 
    % we perform the calculations now and temporarily store them in
    % this matrix, and later assign them to the main matrix.
    Mat_GARP = zeros (obs_num, 6);

    % THIS CODE (for out-of-sample of WARP/GARP/SARP) IS FASTER, but its advantage is negligible
    % (e.g. running time of calculations is 1 seconds instead of 3 seconds,
    %  but printing the residuals to the file takes 100 seconds).
    % so we use the regular, more naive but simpler code ahead (recursive calls). 
%     % calculation of out-of-sample residuals
%     [subsets_of_observations, num_of_subsets] = HPZ_Indices_Problem_Distribute(RP, RP);
%     for i=1:num_of_subsets
%         
%         current_subset = subsets_of_observations{i};
%         
%         if length(current_subset) == 1
%             
%             Mat_GARP(current_subset(1), 1) = WARP_full;     % WARP partial index
%             Mat_GARP(current_subset(1), 3) = GARP_full;     % GARP partial index
%             Mat_GARP(current_subset(1), 5) = SARP_full;     % SARP partial index
%             Mat_GARP(current_subset(1), 2) = 0;             % WARP difference in %
%             Mat_GARP(current_subset(1), 4) = 0;             % GARP difference in %
%             Mat_GARP(current_subset(1), 6) = 0;             % SARP difference in %
%             
%         else
%             
%             subset_DRP = DRP(current_subset, current_subset);
%             subset_SDRP = SDRP(current_subset, current_subset);
%             subset_identical_choice = identical_choice(current_subset, current_subset);
%             
%             for j=1:length(current_subset)
%                 
%                 remaining_indices = [1:(j-1) , (j+1):length(current_subset)];
%                 truncated_identical_choice = subset_identical_choice(remaining_indices, remaining_indices);
%                 truncated_DRP = subset_DRP(remaining_indices, remaining_indices);
%                 truncated_SDRP = subset_SDRP(remaining_indices, remaining_indices);
%                 
%                 % WARP
%                 truncated_WARP = truncated_DRP.*(truncated_DRP') - truncated_identical_choice;
%                 %truncated_WARP_ERRORS = sum(sum(truncated_WARP));
%                 truncated_WARP_VIO_PAIRS = sum(sum(triu(truncated_WARP|(truncated_WARP'))));
%                 % GARP
%                 [truncated_GARP, truncated_RP] = GARP_based_on_DRP_and_SDRP(truncated_DRP, truncated_SDRP);
%                 %truncated_GARP_ERRORS = sum(sum(truncated_GARP));
%                 truncated_GARP_VIO_PAIRS = sum(sum(triu(truncated_GARP|(truncated_GARP'))));
%                 % SARP
%                 truncated_SARP = truncated_RP.*(truncated_RP') - truncated_identical_choice;
%                 %truncated_SARP_ERRORS = sum(sum(truncated_SARP));
%                 truncated_SARP_VIO_PAIRS = sum(sum(triu(truncated_SARP|(truncated_SARP'))));
%                 
%                 % assignment
%                 Mat_GARP(current_subset(j), 1) = truncated_WARP_VIO_PAIRS;      % WARP partial index
%                 Mat_GARP(current_subset(j), 3) = truncated_GARP_VIO_PAIRS;      % GARP partial index
%                 Mat_GARP(current_subset(j), 5) = truncated_SARP_VIO_PAIRS;      % SARP partial index
%                 Mat_GARP(current_subset(j), 2) = (WARP_full - truncated_WARP_VIO_PAIRS) / WARP_full;    % WARP difference in %
%                 Mat_GARP(current_subset(j), 4) = (GARP_full - truncated_GARP_VIO_PAIRS) / GARP_full;    % GARP difference in %
%                 Mat_GARP(current_subset(j), 6) = (SARP_full - truncated_SARP_VIO_PAIRS) / SARP_full;    % SARP difference in %
%                 
%             end
%         end
%         
%     end

    % calculation of out-of-sample residuals
    for i=1:obs_num
        truncated_data = squeeze(truncated_datas(i,:,:));
        [~,~,~,~,~,~,~,~,~,VIO_PAIRS_temp,VIOLATIONS_temp,~,~,~,~,~,~,~] = HPZ_Subject_Consistency (truncated_data, [1,0,0,0,0,0,0], [0,0,0,0], [0,0,0,0,0], [0,0,0,0], [0,0,0,0], 0); %#ok<ASGLU>
        WARP_partial = VIO_PAIRS_temp(1);
        GARP_partial = VIO_PAIRS_temp(2);
        SARP_partial = VIO_PAIRS_temp(3);   
        % WARP
        Mat_GARP(i, 1) = WARP_partial;                              % partial index
        Mat_GARP(i, 2) = (WARP_full - WARP_partial) / WARP_full;    % difference in %
        % GARP
        Mat_GARP(i, 3) = GARP_partial;                              % partial index
        Mat_GARP(i, 4) = (GARP_full - GARP_partial) / GARP_full;    % difference in %
        % SARP
        Mat_GARP(i, 5) = SARP_partial;                              % partial index
        Mat_GARP(i, 6) = (SARP_full - SARP_partial) / SARP_full;    % difference in %
    end
end

if (residuals_flags(1) == 1)
    
    % WARP
    if (GARP_flags(5) == 1 && WARP_full ~= 0)
        % full value
        Mat(:, col_counter) = WARP_full;
        % update the column
        col_counter = col_counter + 1;
        
        % in sample
        if (in_sample_flags(1) == 1)
            WARP_Pairs = triu(WARP|(WARP'));
            % this loop go through all the observations and finds the residual for each
            for i=1:obs_num
                for j=i:obs_num
                    if WARP_Pairs(i,j) == 1
                        % if there is a GARP violation involving observations
                        % i and j, increase the in-sample residual for i and j
                        Mat(i, col_counter) = Mat(i, col_counter) + 1;  
                        Mat(j, col_counter) = Mat(j, col_counter) + 1;
                    end
                end
            end
            
            % this loop finds the residaul as (%) of the total WARP
            for i=1:obs_num
                Mat(i, col_counter+1) = Mat(i, col_counter) / WARP_VIO_PAIRS;
            end
            
            % update the column
            col_counter = col_counter + 2;
        end
        
        % out of sample
        if (out_sample_flags(1) == 1)
            for i=1:obs_num
                Mat(i, col_counter) = Mat_GARP(i, 1);
                Mat(i, col_counter+1) = Mat_GARP(i, 2);
            end
            
            % update the column
            col_counter = col_counter + 2;
        end
    end
        
    % GARP
    if (GARP_flags(6) == 1 && GARP_full ~= 0)
        % full value
        Mat(:, col_counter) = GARP_full;
        % update the column
        col_counter = col_counter + 1;
        
        % in sample
        if (in_sample_flags(1) == 1)
            GARP_Pairs = triu(GARP|(GARP'));
            % this loop go through all the observations and finds the residual for each
            for i=1:obs_num
                for j=i:obs_num
                    if GARP_Pairs(i,j) == 1
                        % if there is a GARP violation involving observations
                        % i and j, increase the in-sample residual for i and j
                        Mat(i, col_counter) = Mat(i, col_counter) + 1;  %#ok<*AGROW>
                        Mat(j, col_counter) = Mat(j, col_counter) + 1;
                    end
                end
            end
            
            % this loop finds the residaul as (%) of the total GARP
            for i=1:obs_num
                Mat(i, col_counter+1) = Mat(i, col_counter) / GARP_VIO_PAIRS;
            end
            
            % update the column
            col_counter = col_counter + 2;
        end
        
        % out of sample
        if (out_sample_flags(1) == 1)
            for i=1:obs_num
                Mat(i, col_counter) = Mat_GARP(i, 3);
                Mat(i, col_counter+1) = Mat_GARP(i, 4);
            end
            
            % update the column
            col_counter = col_counter + 2;
        end
    end
        
    % SARP
    if (GARP_flags(7) == 1 && SARP_full ~= 0)
        % full value
        Mat(:, col_counter) = SARP_full;
        % update the column
        col_counter = col_counter + 1;
        
        % in sample 
        if (in_sample_flags(1) == 1)
            SARP_Pairs = triu(SARP|(SARP'));
            % this loop go through all the observations and finds the residual for each
            for i=1:obs_num
                for j=i:obs_num
                    if SARP_Pairs(i,j) == 1
                        % if there is a GARP violation involving observations
                        % i and j, increase the in-sample residual for i and j
                        Mat(i, col_counter) = Mat(i, col_counter) + 1;  %#ok<*AGROW>
                        Mat(j, col_counter) = Mat(j, col_counter) + 1;
                    end
                end
            end
            
            % this loop finds the residaul as (%) of the total SARP
            for i=1:obs_num
                Mat(i, col_counter+1) = Mat(i, col_counter) / SARP_VIO_PAIRS;
            end
            
            % update the column
            col_counter = col_counter + 2;
        end
        
        % out of sample
        if (out_sample_flags(1) == 1)
            for i=1:obs_num
                Mat(i, col_counter) = Mat_GARP(i, 5);
                Mat(i, col_counter+1) = Mat_GARP(i, 6);
            end
            
            % update the column
            col_counter = col_counter + 2;
        end
    end
    
end



% updating waitbar after finishing with GARP
if (active_waitbar)
    new_bar_val = HPZ_Constants.waitbar_finish_GARP;
    waitbar(new_bar_val, h_wb, {waitbar_msg, char(strcat({'Completed:'}, indices_strings{1}, {' , Remaining:'}, indices_strings{2:5}))}); %#ok<*NODEF>
end



%% inconsistency indices

% initializations of inconsistency indices results variables

AFRIAT = 0;

VARIAN = [0,0,0];

HM = 0;

HM_exact = 1;

Var_exact = 1;

MPI = [0,0];

% If the data doesn't satisfies GARP, then the program calculates the violation 
% extent according to AFRIAT, Varian and Houtman-Maks indices (which of these
% that the user chose). 
% When it can't calculate an index in a reasonable time, it returns an estimation.
if GARP_FLAG == 1 
    
    %% AFRIAT
    if (AFRIAT_flags(1) == 1)
        % if AFRIAT was chosen
        % (otherwise we don't want to aimlessly waste time in unneeded calculations)
        AFRIAT = HPZ_Afriat_efficiency_index (expenditure, index_threshold);
        
        % AFRIAT residuals
        if AFRIAT_flags(2) == 1
        
            % in sample (not implemented)
            % if AFRIAT_flags(3) == 1
            % 
            %     % we calculate for each observation its residual
            %     for i=1:obs_num
            %         Mat(i, col_counter) = 0;   % here there should be a call to a function that calculates residuals... 
            %     end 
            %     % update the counter
            %     col_counter = col_counter + 1;
            % 
            % end
            
            % out of sample (the only option actually for residuals)
            if AFRIAT_flags(4) == 1
                
                % we calculate for each observation its residual
                for i=1:obs_num
                    
                    if Mat_GARP(i, 3) == 0
                        % perfectly consistent
                        AFRIAT_partial = 0;
                    else
                        % not perfectly consistent
                        truncated_data = squeeze(truncated_datas(i,:,:));
                        truncated_Choices(:,1:2) = truncated_data(:,3:4);   % quantities
                        truncated_Choices(:,3:4) = truncated_data(:,5:6);   % prices
                        truncated_expenditure = (truncated_Choices(:,1)*truncated_Choices(:,3)' + truncated_Choices(:,2)*truncated_Choices(:,4)')';
                        AFRIAT_partial = HPZ_Afriat_efficiency_index (truncated_expenditure, index_threshold);
                    end
                    Mat(i, col_counter) = AFRIAT;                       % full index
                    Mat(i, col_counter + 1) = AFRIAT_partial;           % partial index
                    Mat(i, col_counter + 2) = AFRIAT - AFRIAT_partial;  % difference
                end 
                % update the counter
                col_counter = col_counter + 3;

            end
            
        end
    end
    
    
    
    % updating waitbar after finishing with AFRIAT
    if (active_waitbar)
        new_bar_val = HPZ_Constants.waitbar_finish_AFRIAT;
        waitbar(new_bar_val, h_wb, {waitbar_msg, char(strcat({'Completed:'}, indices_strings{1:2}, {' , Remaining:'}, indices_strings{3:5}))});
    end
    
    
    
    %% VARIAN
    if (VARIAN_flags(1) == 1)
        % if VARIAN was chosen
        % (otherwise we don't want to aimlessly waste time in unneeded calculations)
        %[VARIAN, Var_exact, VAR_in_sample_residuals, ~] = HPZ_Varian_efficiency_index (expenditure, identical_choice, index_threshold);
        [VARIAN, Var_exact, VAR_in_sample_residuals, VAR_out_sample_residuals, ~] = HPZ_Varian_Manager (expenditure, identical_choice, index_threshold, SDRP, VARIAN_flags(2)&&VARIAN_flags(4));
        
        % VARIAN residuals
        if VARIAN_flags(2) == 1

            % in sample
            if VARIAN_flags(3) == 1

                % we assign each observation its residuals (mean and meanssq) 
                for i=1:obs_num
                    Mat(i, col_counter) = VAR_in_sample_residuals(i,1);         % MIN
                    Mat(i, col_counter + 1) = VAR_in_sample_residuals(i,2);     % MEAN 
                    Mat(i, col_counter + 2) = VAR_in_sample_residuals(i,3);     % AVGSSQ
                end 
                % update the counter
                col_counter = col_counter + 3;

            end

            % out of sample
            if VARIAN_flags(4) == 1

                % we assign each observation its residuals (mean and meanssq) 
                for i=1:obs_num
                    % normalized indices (relevant for MEAN and AVGSSQ):
                    VARIAN_Min_norm = VAR_out_sample_residuals(i,1);
                    VARIAN_Mean_norm = VAR_out_sample_residuals(i,2) * (rows-1)/rows;
                    VARIAN_AVGSSQ_norm = sqrt(VAR_out_sample_residuals(i,3)^2 * (rows-1)/rows);
                    % assigning to the matrix
                    Mat(i, col_counter) = VARIAN(1);                            % MIN full index
                    Mat(i, col_counter + 1) = VARIAN(2);                        % MEAN full index
                    Mat(i, col_counter + 2) = VARIAN(3);                        % AVGSSQ full index
                    Mat(i, col_counter + 3) = Var_exact;                        % full indices exact or not
                    Mat(i, col_counter + 4) = VAR_out_sample_residuals(i,1);    % MIN partial index
                    Mat(i, col_counter + 5) = VAR_out_sample_residuals(i,2);    % MEAN partial index
                    Mat(i, col_counter + 6) = VAR_out_sample_residuals(i,3);    % AVGSSQ partial index
                    Mat(i, col_counter + 7) = VAR_out_sample_residuals(i,4);    % partial indices exact or not
                    Mat(i, col_counter + 8) = VARIAN(1) - VAR_out_sample_residuals(i,1);    % MIN difference
                    Mat(i, col_counter + 9) = VARIAN(1) - VARIAN_Min_norm;                  % MIN normalized difference
                    Mat(i, col_counter + 10) = VARIAN(2) - VAR_out_sample_residuals(i,2);   % MEAN difference
                    Mat(i, col_counter + 11) = VARIAN(2) - VARIAN_Mean_norm;                % MEAN normalized difference
                    Mat(i, col_counter + 12) = VARIAN(3) - VAR_out_sample_residuals(i,3);   % AVGSSQ difference
                    Mat(i, col_counter + 13) = VARIAN(3) - VARIAN_AVGSSQ_norm;              % AVGSSQ normalized difference
                end 
                % update the counter
                col_counter = col_counter + 14;
                
            end

        end
    end
    
    
    
    % updating waitbar after finishing with VARIAN
    if (active_waitbar)
        new_bar_val = HPZ_Constants.waitbar_finish_VARIAN;
        waitbar(new_bar_val, h_wb, {waitbar_msg, char(strcat({'Completed:'}, indices_strings{1:3}, {' , Remaining:'}, indices_strings{4:5}))});
    end
    
    
    
    %% HOUTMAN-MAKS
    if (HOUTMAN_flags(1) == 1)
        % if HOUTMAN-MAKS was chosen
        % (otherwise we don't want to aimlessly waste time in unneeded calculations)
        [HM, HM_residuals, ~, ~, HM_exact] = HPZ_Houtman_Maks_Manager (HOUTMAN_flags(2), DRP, SDRP, RP, SRP, is_2_goods);
        
        %% HOUTMAN-MAKS residuals
        if HOUTMAN_flags(2) == 1
            
            % out of sample (the only option actually for residuals)
            if HOUTMAN_flags(4) == 1
                
                for i=1:obs_num
                    Mat(i, col_counter) = HM;                                       % full index
                    Mat(i, col_counter + 1) = HM_exact;                             % full index exact or not
                    Mat(i, col_counter + 2) = HM_residuals(i);                      % partial index
                    Mat(i, col_counter + 3) = HM_exact;                             % partial index exact or not
                    Mat(i, col_counter + 4) = HM - HM_residuals(i);                 % difference
                    % (we added this check, because the whole point of the
                    % normalized difference is that it is either positive
                    % or zero, but due to calculation issues it sometimes
                    % resulted stuff like "2.77555756156289E-17")
                    normalized_difference = HM - HM_residuals(i)*(rows-1)/rows;
                    if abs(normalized_difference) < 10^(-15)
                        normalized_difference = 0;
                    end
                    Mat(i, col_counter + 5) = normalized_difference;                % normalized difference
                end

            end

        end
    end
    
    
    
    % updating waitbar after finishing with HOUTMAN-MAKS
    if (active_waitbar)
        new_bar_val = HPZ_Constants.waitbar_finish_HOUTMAN;
        waitbar(new_bar_val, h_wb, {waitbar_msg, char(strcat({'Completed:'}, indices_strings{1:4}, {' , Remaining:'}, indices_strings{5}))});
    end
    
    
    
    %% MPI
    % (NOT YET IMPLEMENTED - MPI will appear in the next version (version 3) 
%     if (MPI_flags(1) == 1)
%         % if MPI (Money Pump Index) was chosen
%         % (otherwise we don't want to aimlessly waste time in unneeded calculations)
%         MPI = HPZ_Money_Pump_Index(expenditure, DRP);
%         
%         % MPI residuals
%         if MPI_flags(2) == 1
%             
%             % in sample (not implemented)
%             % if MPI_flags(3) == 1
%             % 
%             %     % we calculate for each observation its residual
%             %     for i=1:obs_num
%             %         Mat(i, col_counter) = 0;   % here there should be a call to a function that calculates residuals... 
%             %     end 
%             %     % update the counter
%             %     col_counter = col_counter + 1;
%             % 
%             % end
%             
%             % out of sample (currently the only option actually for residuals)
%             if MPI_flags(4) == 1
%                 
%                 % we calculate for each observation its residual
%                 for i=1:obs_num
%                     
%                     if Mat_GARP(i, 3) == 0
%                         % perfectly consistent
%                         MPI_partial = [0,0];
%                     else
%                         % not perfectly consistent
%                         truncated_data = squeeze(truncated_datas(i,:,:));
%                         truncated_Choices(:,1:2) = truncated_data(:,3:4);   % quantities
%                         truncated_Choices(:,3:4) = truncated_data(:,5:6);   % prices
%                         truncated_expenditure = (truncated_Choices(:,1)*truncated_Choices(:,3)' + truncated_Choices(:,2)*truncated_Choices(:,4)')';
%                         [~ , ~ , truncated_DRP , ~] = GARP_based_on_expenditures(truncated_expenditure, [], index_threshold);
%                         MPI_partial = HPZ_Money_Pump_Index(truncated_expenditure, truncated_DRP);
%                     end
%                     
%                     Mat(i, col_counter) = MPI(1);                       % full index mean
%                     Mat(i, col_counter + 1) = MPI(2);                   % full index median
%                     Mat(i, col_counter + 2) = MPI_partial(1);           % partial index mean
%                     Mat(i, col_counter + 3) = MPI_partial(2);           % partial index median
%                     Mat(i, col_counter + 4) = MPI(1) - MPI_partial(1);  % difference mean
%                     Mat(i, col_counter + 5) = MPI(2) - MPI_partial(2);  % difference median
%                 end 
%                 % update the counter
%                 col_counter = col_counter + 6; %#ok<NASGU>
% 
%             end
%             
%         end
%     end
%     
%     % updating waitbar after finishing with HOUTMAN-MAKS
%     if (active_waitbar)
%         new_bar_val = HPZ_Constants.waitbar_finish_MPI;
%         waitbar(new_bar_val, h_wb, {waitbar_msg, char(strcat({'Completed:'}, indices_strings{1:5}))});
%     end
end



if (active_waitbar)
    % close the waitbar
    close(h_wb, h_wb);
end



end

