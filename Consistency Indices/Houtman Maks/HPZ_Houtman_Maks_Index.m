function [HM_mean, HM_sum, HM_raw_residuals] = HPZ_Houtman_Maks_Index (DRP, SDRP, identical_choice, out_of_sample_required)

% for this specific algorithm (finding minimal cycles),
% we want DRP to not contain 1 in (i,j) if i and j are identical choices 
DRP = DRP .* (~identical_choice); % (added 20.02.19)

% number of observations
[obs_num , ~] = size(DRP);

% when we use an algorithm to find all minimal cycles, we don't want that 
% there would be edges from a vertex to itself, so we make sure it won't occur
for i=1:obs_num
    DRP(i,i) = 0;
end



% find all minimal cycles (cycles that don't contain sub-cycles)

% --- according to SARP ---
% % NOTE! When using this function to calculate HM according to SARP, the
% % given DRP must NOT include a relation between identical choices, that
% % is: we use DRP & ~identical_choice instead of DRP
% cycles = find_all_minimal_cycles (DRP & ~identical_choice);  

% --- according to GARP ---
cycles = find_all_minimal_GARP_cycles (DRP, SDRP);



% now we make preparations for the use of integer programming (using intlinprog)   

% weight of each observations
weights = ones(1,obs_num);
f = weights';

% all values are integers that receive value of either 0 or 1
intcon = 1:obs_num;
lb = zeros(1, obs_num);
ub = ones(1, obs_num);

% the Ax<=b restriction will make sure that each cycle is covered by at least one observation 
cycles_num = length(cycles);
A = zeros(cycles_num, obs_num);  % initialization
b = - ones(cycles_num, 1);
for c = 1:cycles_num
    currect_cycle = cycles{c};
    A(c, currect_cycle) = -1;
end

% avoid displaying algorithm stages in command window
intlinprog_options = optimoptions('intlinprog', 'Display','off');

% intlinprog function was changed between R2017a and R2017b,
% hence we need to address that. we don't bother addressing version before 2014a, 
% since intlinprog wasn't even introduced before then (there was bintprog instead).  
is_old_version_intlinprog = all(version('-release') == '2017a') || ...
                all(version('-release') == '2016b') || all(version('-release') == '2016a') || ...
                all(version('-release') == '2015b') || all(version('-release') == '2015a') || ...
                all(version('-release') == '2014b') || all(version('-release') == '2014a');

% find the raw HM Index (minimum number of observations that need to be dropped) 
if is_old_version_intlinprog
    [x, fval, exitflag, ~] = intlinprog(f, intcon, A,b, [],[], lb,ub, intlinprog_options);
else
    [x, fval, exitflag, ~] = intlinprog(f, intcon, A,b, [],[], lb,ub, [], intlinprog_options);
end

if exitflag == 1
    % assigning to the matrices
    HM_sum = fval;
    HM_mean = fval / obs_num;
    solution_vector = x;
else
    % assigning nan values and ending the function
    HM_sum = nan;
    HM_mean = nan;
    HM_raw_residuals = nan(1, obs_num);
    return
end





% out-of-sample

% initialization, and also this is what will be returned if residuals are not required   
HM_raw_residuals = nan(1, obs_num); 
% we need to reduce the sizes of these vectors by one
% these vectors (except intcon) contain the same value in all cells so we can truncate them any way we want  
intcon_residuals = intcon(1:(end-1));
lb_residuals = lb(2:end);
ub_residuals = ub(2:end);
f_residuals = f(2:end);

if out_of_sample_required
    
%     if obs_num == 1
%         % a single observation cannot violate GARP, so...
%         HM_raw_residuals(i) = 0;
    
    for i=1:obs_num
        
        if solution_vector(i) == 1
            % then we know that dropping this observations will reduce by 1
            % the number of observations needed to be dropped  
            % (this check is not necessary, but it saves some time)  
            HM_raw_residuals(i) = 1;
        else
            % we need to check if dropping this observation reduces by 1 or
            % not reduces the number of observations needed to be dropped
            
            % cycles that observation i is part of them
            cycles_that_include_observation_i = logical(A(:,i));
            
            % a matrix A without observation i and without all minimal
            % cycles that include observation i
            A_without_obs_i = A(~cycles_that_include_observation_i , [1:(i-1),(i+1):end]);
            b_without_obs_i = b(~cycles_that_include_observation_i);
            
            if isempty(A_without_obs_i)
                % then removing this observations removed all the cycles,
                % and since we only use this function on datasets that
                % violate GARP, we know that dropping this observation
                % reduces the number of observations that need to be dropped (by one)  
                HM_raw_residuals(i) = 1;
            end
            
            % find the raw HM Index without this observations 
            [~, fval, exitflag, ~] = intlinprog(f_residuals, intcon_residuals, A_without_obs_i,b_without_obs_i, [],[], lb_residuals,ub_residuals, [], intlinprog_options);
            if exitflag == 1
                % assigning to the raw residuals vector
                HM_raw_residuals(i) = HM_sum - fval;   % will be either 1 or 0 
            end
        end
        
    end
    
end





end