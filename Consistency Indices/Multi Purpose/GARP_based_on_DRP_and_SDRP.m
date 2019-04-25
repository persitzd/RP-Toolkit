function [GARP, RP, SRP] = GARP_based_on_DRP_and_SDRP(DRP, SDRP)

% number of observations
[obs_num, ~] = size(DRP);

% statement needed for the graph theory external package to work efficiently
set_matlab_bgl_default(struct('full2sparse',1));

% The matrix NS_RP has at the cell in the i'th row and the j'th column, Inf if and only if the 
% bundle of observation i is not revealed preferred to the bundle of observation j. 
% otherwise it includes a positive integer.
NS_RP = all_shortest_paths(DRP);

% The matrix RP has at the cell in the i'th row and the j'th column, 
% 1 if and only if the bundle i is revealed preferred to the bundle j. 
% Otherwise, it equals 0.      
RP = ~isinf(NS_RP)*1;



% The matrix NS_SRP has at the cell in the i'th row and the j'th column, Inf if and only if the 
% bundle of observation i is either not revealed preferred to the bundle of observation j, 
% or it is revealed preferred but the chain of direct revealed preference relations 
% does not contain any strict relation.
% otherwise it includes a positive integer.
NS_SRP = (RP*SDRP)*RP;

% The matrix SRP has at the cell in the i'th row and the j'th column, 
% 1 if and only if the bundle i is strictly revealed preferred to the bundle j. 
% Otherwise, it equals 0.  
SRP = (NS_SRP ~= 0)*1;



% THIS MUST BE *AFTER* NS_SRP IS CALCULATED. OTHERWISE, NS_SRP MAY NOT
% INCLUDE SOME CASES WHERE THERE IS A DIRECT SDRP BETWEEN OBSERVATIONS
% for every practical use, we have no need in the DRP relation between an
% observation to itself. more than that, if the DRP matrix will represent 
% this relation, it may hinder some of the function and calculations.
for i=1:obs_num
    RP(i,i)=0;
end



% To test for GARP we will do the following: for every pair of
% choices x and y if xRy then not yP0x. We will take RP and the transpose of
% SDRP and multiply element by element. Every 1 corresponds to 
% xRy and yP0x. The final matrix is the zero matrix if and only if 
% GARP is satisfied. 
% The matrix is the zero matrix if and only if GARP is satisfied.
GARP = RP .* (SDRP');



end