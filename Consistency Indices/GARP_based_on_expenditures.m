function [GARP, RP, DRP, SDRP] = GARP_based_on_expenditures(expenditure, identical_choice, index_threshold)

% number of observations
[obs_num, ~] = size(expenditure);

% if the user didn't specify the identical choice matrix
% (the user should use "[]" as input in this situation),
% we initilize the matrix with no identical choices
if isempty(identical_choice)
    identical_choice = zeros(obs_num, obs_num);
end

% The matrix REF has at the cell in the i'th row and the j'th column, the difference between the value of the bundle that 
% was chosen in observation i and the bundle that was chosen in observation j given the prices of observation i
REF = diag(expenditure) * ones(obs_num,1)' - expenditure;

% DRP - represents the relation R^0. The matrix DRP has at the cell in the i'th row and the j'th
% column 1 if and only if the bundle that was chosen in observation i is directly revealed 
% preferred to the bundle that was chosen in observation j. otherwise it equals 0.
%DRP = ceil((REF+index_threshold) / (max(max(abs(REF+index_threshold)))+1));
DRP = (REF + identical_choice + diag(expenditure)*ones(obs_num,1)'*index_threshold >= 0) * 1;

% SDRP - represents the relation P^0. 
% The matrix SDRP has at the cell in the i'th row and the j'th column, 1 if and only if the bundle
% that was chosen in observation i is strictly directly revealed preferred to the bundle that
% was chosen in observation j. otherwise it equals 0.
%SDRP = ceil((REF-index_threshold) / (max(max(abs(REF-index_threshold)))+1));
SDRP = (REF - diag(expenditure)*ones(obs_num,1)'*index_threshold > 0) * 1;

% The matrix is the zero matrix if and only if GARP is satisfied.
[GARP, RP] = GARP_based_on_DRP_and_SDRP(DRP, SDRP);

end