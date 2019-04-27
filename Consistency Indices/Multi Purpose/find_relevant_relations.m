function [relevant_DRP, relevant_SDRP] = find_relevant_relations(DRP, SDRP)

% this function accepts the R0 (DRP) and P0 (SDRP) relations as matrix in
% which there is a relation between i and j iff mat(i,j) = 1,
% and returns a matrix in which mat(i,j) = 1 only for relations that are a
% part of a GARP violation cycle.

% "DRP" stands for "Directly Revealed Preferred".
% "SDRP" stands for "Strictly Directly Revealed Preferred".
% for our purposes, we take the DRP not including identical choices.

[obs_num , ~] = size(DRP);

% calculate RP and SRP
[~, RP, SRP] = GARP_based_on_DRP_and_SDRP(DRP, SDRP);

% % the above function returns RP without a bundle revealed preferred to itself, but here we need it 
% for i=1:obs_num
%     RP(i,i)=1;
% end
% 
% % calculate SRP
% NS_SRP = (RP*SDRP)*RP;
% SRP = (NS_SRP ~= 0);

% a DRP relation (i,j) belongs to some GARP cycle iff there is SRP relation (j,i) 
% (Note: it could be that it is only a part of a non-simple GARP cycle
% (i.e. some observation appears more than once in the cycle), but I
% haven't found a way to handle it. So it may contain DRP relations that
% are not a part of a simple GARP cycle)
relevant_DRP = DRP .* SRP';

% a SDRP relation (i,j) belongs to some GARP cycle iff there is RP relation (j,i) 
relevant_SDRP = SDRP .* RP';

for k=1:obs_num
    % we need that the process will not break the relation i-DRP-i for
    % every i that this relation existed before
    relevant_DRP(k,k) = DRP(k,k);
end

% a bundle is never directly strictly revealed preferred to itself (just in case) 
for i=1:obs_num
    relevant_SDRP(i,i)=0;
end

% we want to include every relevant relation - be it a DRP or SDRP relation 
%relevant_relations = max(DRP, SDRP);

end