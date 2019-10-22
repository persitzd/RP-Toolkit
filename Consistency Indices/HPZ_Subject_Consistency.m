function [expenditure, DRP, SDRP, RP, SRP, WARP, GARP, SARP, FLAGS, VIO_PAIRS, VIOLATIONS, AFRIAT, VARIAN_Bounds, HoutmanMaks, MPI, Mat] = HPZ_Subject_Consistency (data_matrix, Graph_flags, GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags, Mat_num_of_columns, Varian_algorithm_settings, main_folder_for_results, active_waitbar, residuals_waitbar, current_run, total_runs)

% this function calculates all the currently implemented consistency and
% inconsistency indices, and their residuals (and depending on the flags 
% it may calculate only some of them), for a single subject.
% the WARP, GARP and SARP measures are implemented within the function
% itself, while the inconsistency indices (Afriat, Varian, Houtman-Maks) 
% are being calculated by calling other specific functions.

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
%   (1) Max, (2) Mean, (3) Average Sum of Squares.
%   Var_exact and HM_exact tell the user whether the calculation was exact
%   or was it only an approximation, and what type of approximation.
% Mat :
%   a matrix with all the residuals results as should be printed to the
%   results file




% number of observations and calculation of number of goods
[obs_num, num_of_columns] = size(data_matrix);
num_of_goods = (num_of_columns - 2) / 2;



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



% these are for convenience
basic_flags = [GARP_flags(1) , AFRIAT_flags(1) , VARIAN_flags(1) , HOUTMAN_flags(1) , MPI_flags(1)];
residuals_flags_temp = [GARP_flags(2) , AFRIAT_flags(2) , VARIAN_flags(2) , HOUTMAN_flags(2) , MPI_flags(2)];
residuals_flags = min(basic_flags , residuals_flags_temp);
in_sample_flags = [GARP_flags(3) , AFRIAT_flags(3) , VARIAN_flags(3) , HOUTMAN_flags(3) , MPI_flags(3)]; %#ok<NASGU>
out_sample_flags = [GARP_flags(4) , AFRIAT_flags(4) , VARIAN_flags(4) , HOUTMAN_flags(4) , MPI_flags(4)];

% initializations as empty matrices to avoid problems and make the code easier (less if statements) 
Base_Mat = [];
WARP_GARP_SARP_Mat = [];
Afriat_Mat = [];
Varian_Mat = [];
HM_Mat = [];
MPI_Mat = [];





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




%% initialization of the base of the residuals result Matrix
if sum(residuals_flags) ~= 0
    % the residuals are required
    % creating the first part of the results matrix for the residuals -
    % the one that contains the subject number and observation number
    Base_Mat = zeros (obs_num, 2);
    % entering the subject ID and the observations numbers
    for i=1:obs_num
        Base_Mat(i, 1) = data_matrix(1,1);
        Base_Mat(i, 2) = i;
    end 
end





% Choices
%Choices = data_matrix(1:obs_num, 3:6);   % quantities & prices
Choices = data_matrix(1:obs_num, 3:end);   % quantities & prices   %EXTENSION 
% [~ , num_of_goods] = size(Choices)/2;   %EXTENSION 

% The matrix "expenditure" has at the cell in the i'th row and the j'th
% column, the value of the bundle that was chosen in observation j given the
% prices of observation i
%expenditure = (Choices(:,1)*Choices(:,3)' + Choices(:,2)*Choices(:,4)')';
expenditure = (Choices(:, 1:num_of_goods) * Choices(:, (num_of_goods+1):end)')';  %EXTENSION 

% the purpose of the matrix identical_choice is to locate identical choices. Value of 
% cell (j,k) is 1 if choices identical, and 0 otherwise. The following loop builds the 
% matrix: 
identical_choice = zeros(obs_num, obs_num);

for j=1:obs_num
    % going through all identical_choice’s cells.
    for k=1:obs_num
        % if the choices j&k are identical, then set value of cell (j,k) = 1 , otherwise 0. 
        if all( abs(Choices(j, 1:num_of_goods) - Choices(k, 1:num_of_goods)) <= Choices(j, 1:num_of_goods)*identical_choices_threshold )
        % if all( (abs(Choices(j,:) - Choices(k,:)) <= Choices(j,:)*identical_choices_threshold) )   %EXTENSION 
            identical_choice(j,k) = 1;
        end
    end
end





%  create graph/s of revealed preference relations, and save it (if the user specified so) 
if any(Graph_flags(1:4))
    HPZ_draw_RP_graphs (Graph_flags, main_folder_for_results, data_matrix(1,1), expenditure, identical_choice, index_threshold, Choices(:, 1:num_of_goods));
end




% finding GARP, and the relations RP, SRP, DRP and SDRP
[GARP, RP, SRP, DRP, SDRP] = GARP_based_on_expenditures(expenditure, identical_choice, index_threshold);

% calculating WARP GARP and SARP violations, and assigning the results to the results vectors 
[FLAGS, VIO_PAIRS, VIOLATIONS, WARP, SARP] = HPZ_WARP_GARP_SARP(identical_choice, DRP, RP, GARP);
WARP_VIO_PAIRS = VIO_PAIRS(1); %#ok<NASGU>
GARP_VIO_PAIRS = VIO_PAIRS(2);
SARP_VIO_PAIRS = VIO_PAIRS(3); %#ok<NASGU>
WARP_VIOLATIONS = VIOLATIONS(1);
GARP_VIOLATIONS = VIOLATIONS(2);
SARP_VIOLATIONS = VIOLATIONS(3);
GARP_FLAG = FLAGS(2);


%% residuals for WARP/GARP/SARP VIO PAIRS

% the full number of GARP VIO PAIRS violations
WARP_full = WARP_VIOLATIONS; %WARP_VIO_PAIRS;
GARP_full = GARP_VIOLATIONS; %GARP_VIO_PAIRS;
GARP_full_pairs = GARP_VIO_PAIRS;
SARP_full = SARP_VIOLATIONS; %SARP_VIO_PAIRS;

Mat_GARP = [];
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
    Mat_GARP = zeros (obs_num, 8);

    % calculation of out-of-sample residuals
    for i=1:obs_num
        truncated_identical_choice = identical_choice([1:(i-1) , (i+1):end] , [1:(i-1) , (i+1):end]);
        truncated_DRP = DRP([1:(i-1) , (i+1):end] , [1:(i-1) , (i+1):end]);
        truncated_SDRP = SDRP([1:(i-1) , (i+1):end] , [1:(i-1) , (i+1):end]);
        [truncated_GARP, truncated_RP, ~] = GARP_based_on_DRP_and_SDRP(truncated_DRP, truncated_SDRP);
        [~, VIO_PAIRS_temp, VIOLATIONS_temp, ~, ~] = HPZ_WARP_GARP_SARP(truncated_identical_choice, truncated_DRP, truncated_RP, truncated_GARP);
        WARP_partial = VIOLATIONS_temp(1); %VIO_PAIRS_temp(1);
        GARP_partial = VIOLATIONS_temp(2); %VIO_PAIRS_temp(2);
        GARP_partial_pairs = VIO_PAIRS_temp(2);
        SARP_partial = VIOLATIONS_temp(3); %VIO_PAIRS_temp(3);   
        % WARP
        Mat_GARP(i, 1) = WARP_partial;                              % partial index
        Mat_GARP(i, 2) = (WARP_full - WARP_partial) / WARP_full;    % difference in %
        % GARP
        Mat_GARP(i, 3) = GARP_partial;                              % partial index
        Mat_GARP(i, 4) = (GARP_full - GARP_partial) / GARP_full;    % difference in %
        Mat_GARP(i, 5) = GARP_partial_pairs;                                            % partial index
        Mat_GARP(i, 6) = (GARP_full_pairs - GARP_partial_pairs) / GARP_full_pairs;      % difference in %
        % SARP
        Mat_GARP(i, 7) = SARP_partial;                              % partial index
        Mat_GARP(i, 8) = (SARP_full - SARP_partial) / SARP_full;    % difference in %
    end
end

% calculate WARP, GARP and SARP residuals, and assign them to the residuals matrix (and update col_counter)  
if (residuals_flags(1) == 1)
    WARP_GARP_SARP_Mat = HPZ_WARP_GARP_SARP_residuals(GARP_flags, VIO_PAIRS, VIOLATIONS, WARP, GARP, SARP, Mat_GARP);
end



% updating waitbar after finishing with GARP
if (active_waitbar)
    new_bar_val = HPZ_Constants.waitbar_finish_GARP;
    waitbar(new_bar_val, h_wb, {waitbar_msg, char(strcat({'Completed:'}, indices_strings{1}, {' , Remaining:'}, indices_strings{2:5}))}); %#ok<*NODEF>
end



%% inconsistency indices

% initializations of inconsistency indices results variables
AFRIAT = nan;
VARIAN_Bounds = [nan, nan, nan, nan, nan, nan, nan, nan, nan];
HoutmanMaks = nan;
MPI = [nan, nan];

if GARP_FLAG == 0
    % then all the indices are 0 as well
    AFRIAT = 0;
    VARIAN_Bounds = [0, 0, 0, 0, 0, 0, 0, 0, 0];
    HoutmanMaks = 0;
    MPI = [0, 0]; 
end

% If the data doesn't satisfies GARP, then the program calculates the violation 
% extent according to AFRIAT, Varian and Houtman-Maks indices (which of these
% that the user chose). 
% When it can't calculate an index in a reasonable time, it returns an estimation.
if GARP_FLAG == 1 
    
    %% AFRIAT
    if AFRIAT_flags(1)
        % if AFRIAT was chosen
        % (otherwise we don't want to aimlessly waste time in unneeded calculations)
        if AFRIAT_flags(2) && AFRIAT_flags(4)
            GARP_vector = Mat_GARP(:, 3);
        else
            GARP_vector = [];
        end
        % we pass "residuals_waitbar" only when there are a lot of observations, since Afriat is fast enough so normally there is no need in a residuals waitbar  
        residuals_waitbar_AFRIAT = false; %residuals_waitbar & obs_num >= 100;
        [AFRIAT, Afriat_Mat] = HPZ_Afriat_Manager (AFRIAT_flags, expenditure, index_threshold, SDRP, residuals_waitbar_AFRIAT, current_run, total_runs, data_matrix(1,1), GARP_vector);
    end
    
    
    
    % updating waitbar after finishing with AFRIAT
    if (active_waitbar)
        new_bar_val = HPZ_Constants.waitbar_finish_AFRIAT;
        waitbar(new_bar_val, h_wb, {waitbar_msg, char(strcat({'Completed:'}, indices_strings{1:2}, {' , Remaining:'}, indices_strings{3:5}))});
    end
    
    
    
    %% VARIAN
    if VARIAN_flags(1)
        % if VARIAN was chosen
        % (otherwise we don't want to aimlessly waste time in unneeded calculations)
        [VARIAN_Bounds, ~, ~, ~, Varian_Mat] = HPZ_Varian_Manager (VARIAN_flags, expenditure, identical_choice, index_threshold, SDRP, Varian_algorithm_settings, main_folder_for_results);   % , residuals_waitbar, current_run, total_runs, data_matrix(1,1)
    end
    
    
    
    % updating waitbar after finishing with VARIAN
    if (active_waitbar)
        new_bar_val = HPZ_Constants.waitbar_finish_VARIAN;
        waitbar(new_bar_val, h_wb, {waitbar_msg, char(strcat({'Completed:'}, indices_strings{1:3}, {' , Remaining:'}, indices_strings{4:5}))});
    end
    
    
    
    %% HOUTMAN-MAKS
    if HOUTMAN_flags(1)
        % if HOUTMAN-MAKS was chosen
        % (otherwise we don't want to aimlessly waste time in unneeded calculations)
        [HoutmanMaks, ~, ~, ~, HM_Mat] = HPZ_Houtman_Maks_Manager (HOUTMAN_flags, DRP, SDRP, RP, identical_choice);
    end
    
    
    
    % updating waitbar after finishing with HOUTMAN-MAKS
    if (active_waitbar)
        new_bar_val = HPZ_Constants.waitbar_finish_HOUTMAN;
        waitbar(new_bar_val, h_wb, {waitbar_msg, char(strcat({'Completed:'}, indices_strings{1:4}, {' , Remaining:'}, indices_strings{5}))});
    end
    
    
    
    %% MPI
    % (NOT YET IMPLEMENTED - MPI will appear in the next version (version 4) 
%     if (MPI_flags(1) == 1)
%         % if MPI (Money Pump Index) was chosen
%         % (otherwise we don't want to aimlessly waste time in unneeded calculations)
%         MPI = HPZ_Money_Pump_Manager(MPI_flags, expenditure, DRP);
%     
%     % updating waitbar after finishing with HOUTMAN-MAKS
%     if (active_waitbar)
%         new_bar_val = HPZ_Constants.waitbar_finish_MPI;
%         waitbar(new_bar_val, h_wb, {waitbar_msg, char(strcat({'Completed:'}, indices_strings{1:5}))});
%     end

    % joining all residuals matrices of the various indices, to one matrix
    Mat = [Base_Mat , WARP_GARP_SARP_Mat, Afriat_Mat , Varian_Mat , HM_Mat, MPI_Mat];
    
else
    
    % if it satisfies GARP, we fill all requested indices with 0s
    Mat = [Base_Mat, zeros(obs_num, Mat_num_of_columns-2)];
    
end





if (active_waitbar)
    % close the waitbar
    close(h_wb, h_wb);
end



end

