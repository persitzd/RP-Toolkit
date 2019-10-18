function min_cycles = find_all_minimal_cycles (G)

% this function is a Matlab implementation of the FIND ALL MINIMAL CYCLES  
% algorithm from the paper (as of 25.04.19 - working paper): 
% "Code Package for Revealed Preferences Analysis of Choices from Linear Budget Sets"  
%
% for calculation of consistency indices, it can be used for two purposes:  
%   (1) to calculate Houtman-Maks Index according to SARP 
%       (we don't do that, we calculate only according to GARP) 
%   (2) to calculate Varian Index after transforming it to a Weighted Houtman-Maks problem. 
%       this is currently the only place in the program where it is used.



[n, m] = size(G);

% if G is not formatted correctly - create an error
if n ~= m
    error('G must have an equal number of rows and columns')
elseif sum(sum(G == 1)) + sum(sum(G == 0)) ~= n^2
    error('G must contain only values of 0 or 1');
elseif any(diag(G))
    error('G must have 0 in all its diagonal');
end

% initialization
min_cycles = {};

% if G is empty or there are no edges, then there are no cycles - return empty list 
if isempty(G) || sum(sum(G)) == 0
   return
end

% logical improves efficiency (in case the input wasn't logical already)
% sparse seem to not affect running time, but it reduces space complexity 
G = sparse(logical(G));



for i = n:(-1):2
    
    % find all minimal cycles that contain vertex i, 
    % and do not contain vertices i+1,...,n.
    
    % update G
    G = G(1:i , 1:i);
    
    % use the recursive function to find all minimal cycles that contain vertex i 
    min_cycles_with_i = find_all_minimal_cycles_with_path (G, i);
    
    % update the full minimal cycles list
    min_cycles = [min_cycles , min_cycles_with_i]; %#ok<AGROW>
    
end



end