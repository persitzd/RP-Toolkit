function min_cycles = find_all_minimal_GARP_cycles (DRP, SDRP)

% this function is a Matlab implementation of the FIND ALL MINIMAL NON WHITE CYCLES  
% algorithm from the paper (as of 25.04.19 - working paper): 
% "Code Package for Revealed Preferences Analysis of Choices from Linear Budget Sets"  
% in our context:
%   WHITE = R^0 (=DRP, Directly Revealed Preferred)
%   BLACK = P^0 (=SDRP, Strictly Directly Revealed Preferred)
%   and NON-WHITE Cycles are Cycles that violate GARP
%   (while Cyclces in WHITE without any edge in BLACK only violate SARP)



[n, m] = size(DRP);

% if DRP and SDRP are not formatted correctly - create an error
if n ~= m
    error('DRP must have an equal number of rows and columns')
elseif sum(sum(DRP == 1)) + sum(sum(DRP == 0)) ~= n^2
    error('DRP must contain only values of 0 or 1');
elseif any(diag(DRP))
    error('DRP must have 0 in all its diagonal');
elseif ~all(size(DRP) == size(SDRP))
    error('DRP and SDRP must be of same size');
end

% initialization
min_cycles = {};

% if DRP is empty or there are no edges, then there are no cycles - return empty list 
if isempty(DRP) || sum(sum(DRP)) == 0
   return
end


for i = n:(-1):2
    
    % find all minimal cycles that contain vertex i, 
    % and do not contain vertices i+1,...,n.
    
    % update DRP and SDRP
    DRP = DRP(1:i , 1:i);
    SDRP = SDRP(1:i , 1:i);
    
    % use the recursive function to find all minimal cycles that contain vertex i 
    min_cycles_with_i = find_all_minimal_GARP_cycles_with_path (DRP, SDRP, i);
    
    % update the full minimal cycles list
    min_cycles = [min_cycles , min_cycles_with_i]; %#ok<AGROW>
end



end