function [GARP, RP, SRP, DRP, SDRP] = GARP_based_on_expenditures(expenditure, identical_choice, index_threshold)

% NOTE: the DRP relation returned by this function does not include the 
% relation between a bundle to itself (including identical choices in 
% different observations), and the same goes for RP that is calculated on
% the basis of the DRP relation



% number of observations
[obs_num, ~] = size(expenditure);



% if the user didn't specify the identical choice matrix
% (the user should use "[]" as input in this situation),
% we initialize the matrix with no identical choices
if isempty(identical_choice)
    identical_choice = zeros(obs_num, obs_num);
end



% The matrix REF has at the cell in the i'th row and the j'th
% column, the difference between the value of the bundle that was chosen in 
% observation i and the bundle that was chosen in observation j given the 
% prices of observation i
REF = diag(expenditure) * ones(obs_num,1)' - expenditure;



% DRP is a representation of the relation R^0.
% The matrix DRP has at the cell in the i'th row and the j'th
% column, 1 if and only if the bundle that was chosen in 
% observation i is directly revealed preferred to the bundle that was chosen 
% in observation j. otherwise it equals 0.
%
% DRP(i,j) = 1 iff REF_DRP(i,j) > 0 iff the bundle that was chosen
% in observation i is directly revealed preferred to the bundle that was chosen 
% in observation j (for a index_threshold small enough)
%
% We increase REF by a small threshold in order to make sure it is bigger
% than 0, also we add to it "identical_choice" in order to make sure that
% indentical choices will be considered DRP to one another (which is
% crucial for WARP and SARP, otherwise we might have negative numbers in
% WARP and SARP).

%DRP = ceil((REF+index_threshold) / (max(max(abs(REF+index_threshold)))+1));
DRP = (REF + identical_choice + diag(expenditure)*ones(obs_num,1)'*index_threshold >= 0) * 1;

% for every practical use, we have no need in the DRP relation between an
% observation to itself. more than that, if the DRP matrix will represent 
% this relation, it may hinder some of the function and calculations.
for i=1:obs_num
    DRP(i,i)=0;
end



% SDRP - represents the relation P^0. 
% The matrix SDRP has at the cell in the i'th row and the j'th column, 1 if and only if the bundle
% that was chosen in observation i is strictly directly revealed preferred to the bundle that
% was chosen in observation j. otherwise it equals 0.
%
% We decrease REF by a small threshold in order to make sure that if it is 
% bigger than 0, it is because it should be, and not as a result of a
% computational inprecision.

%SDRP = ceil((REF-index_threshold) / (max(max(abs(REF-index_threshold)))+1));
SDRP = (REF - diag(expenditure)*ones(obs_num,1)'*index_threshold > 0) * 1;



% The matrix is the zero matrix if and only if GARP is satisfied.
[GARP, RP, SRP] = GARP_based_on_DRP_and_SDRP(DRP, SDRP);



end