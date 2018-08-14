function [HM, HM_exact] = HPZ_Houtman_Maks_efficiency_index_approx (Choices, index_threshold)

% this function estimates the Houtman-Maks inconsistency index when exact
% calculation has failed.
% the function also returns the residuals of the index.

% If we were not able to directly compute the index, we have to approximate.
% And this is done here using HoutmanMaksCFGK.m that was written by Daniel Martin. 
% The package of Matlab files that calculates the Houtman Maks index was 
% downloaded from Daniel Martin's personal website on November 5th 2011.
% For the implementation I used the readme file that is located next to the
% package on the same website.

% Comment 1: Some adaptation was introduced for the case of two goods

% Comment 2: the approximation written by Daniel Martin works with a single binary relation. However, 
% testing for GARP requires two relations R and P^0. As P^0 and R^0 are often identical, it almost 
% always doesn't matter. However, it is possible of course that xR^0y but not xP^0y, and so 
% we ought to choose a relation for the approximation algorithm. 
% As the approximation algorithm is as an upper bound, we chose "working" with the relation R. 
% The consequence of that, is that we might add cycles that do not exist. 
% While choosing the relation P^0, on the other hand, might cause at missing cycles 
% (those satisfies xRy, yP^0x, but not xPy). So using  R^0 for the algorithm gives 
% us an upper bound of the number of observations one must remove to eliminate all GARP cycles .  

% trials is set to be the number of rows of matrix Choices, meaning m, the 
% number of observations for a single subject.

[trials,~] = size(Choices);

% The matrix "expenditure" has at the cell in the i'th row and the j'th
% column, the value of the bundle that was chosen in observation j given the
% prices of observation i

expenditure = (Choices(:,1)*Choices(:,3)' + Choices(:,2)*Choices(:,4)')';

% The matrix REF has at the cell in the i'th row and the j'th
% column, the difference between the value of the bundle that was chosen in 
% observation i and the bundle that was chosen in observation j given the 
% prices of observation i

REF = diag(expenditure) * ones(trials,1)' - expenditure;

% The matrix DRP has at the cell in the i'th row and the j'th column, 1 if and only if the 
% bundle that was chosen in observation i is directly revealed preferred to the bundle 
% that was chosen in observation j, and 0 otherwise.

DRP = ceil((REF+index_threshold) / (max(max(abs(REF+index_threshold)))+1));


% DRP_num sums the values of DRP cells, meaning the number of cases where 
% a bundle is directly revealed preferred to another bundle.     

DRP_num = sum(sum(DRP));

% input is the matrix of all pairs of bundles number i,j such that x^i is
% directly revealed preferred to x^j. 
% The input matrix has two cells in every row, each contains a number of 
% a bundle, i,j so that x^i is directly revealed preferred to x^j.

input = zeros(DRP_num,2);

% If one wants to insert the p^0 relation use SDRP

%% SDRP = ceil((REF-index_threshold)/(max(max(abs(REF-index_threshold)))+1));

%% SDRP_num = sum(sum(SDRP));

%% input=zeros(SDRP_num,2);

index = 1;

% the purpose of the next loops is to use the DRP matrix in order to insert 
% the matrix "input" every bundles numbers i,j such that x^i is directly 
% revealed preferred to x^j.

for i=1:trials
    for j=1:trials
        
        % if x^i is directly revealed preferred to x^j. insert i,j to as a row 
        % to "input" matrix. Replace DRP with SDRP if P^0 is preferred. 

        if DRP(i,j) == 1
            input(index,1) = j;
            input(index,2) = i;
            index = index+1;
        end
    end
end

% HoutmanMaksCFGK_2G is a function that given the data of R^0, returns an approximation 
% for the minimum number of removals needed for acyclicality in R. 

[SumRemovals , Removals] = HoutmanMaksCFGK_2G(input); %#ok<ASGLU>

% the HM index is calculated as the number of removed relatively to the number of 
% observations.
HM = SumRemovals(1,1) / trials;

% the calculation is not exact, but it is an approximation
HM_exact = 0;

end
