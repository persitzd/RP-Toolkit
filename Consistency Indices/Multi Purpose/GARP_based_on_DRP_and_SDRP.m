function [GARP, RP] = GARP_based_on_DRP_and_SDRP(DRP, SDRP)

% statement needed for the graph theory external package to work efficiently
set_matlab_bgl_default(struct('full2sparse',1));

% The matrix NS_RP has at the cell in the i'th row and the j'th column, Inf if and only if the 
% bundle of observation i is not revealed preferred to the bundle of observation j. 
% otherwise it includes a positive integer.
NS_RP = all_shortest_paths(DRP);

% The matrix RP has at the cell in the i'th row and the j'th column, 
% 1 if and only if the bundle i is revealed preferred to the bundle j. 
% Otherwise, it equals 0.      
RP = ~isinf(NS_RP);

% The matrix is the zero matrix if and only if GARP is satisfied.
GARP = RP .* (SDRP');

end
