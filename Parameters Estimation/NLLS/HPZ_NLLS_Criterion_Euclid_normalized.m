function [normalized_euclid_criterions] = HPZ_NLLS_Criterion_Euclid_normalized(observations, optimal_bundles)

% The function calculates the normalized-Euclidean NLLS criterion for 
% a given set of observations and predictions.
% The function returns the value of the criterion for the given data. 

% Explanation about the normalized euclidean criterion, and its
% calculation:
% In the normalized euclidean criterion, we don't take the value of the
% distance between the chosen bundle and the optimal bundle as is, but we
% normalize it by dividing it by the full length of the budget line. doing
% so makes sure that observations with small endowments or with mild prices 
% ratios will have the same weight as observations with big endowments or
% extreme prices.
% The formula should be therefore: 
% sqrt( sum((observations(:,1:2) - optimal_bundles).^2, 2)) / sqrt( sum(1 ./ (observations(:,3:4)).^2)) 
% But if you think about it, since the budget line is a straight line, in
% order to find the relative length of a sub-part of the budget line, it is
% enough to use only the X values of the given points (or only the Y values), 
% which leads to the following formula.
% (note the 1./observations(:,3) is the max quantity of the first
% commodity, so this formula is equivalent to dividing by the max quantity)

%normalized_euclid_criterions = sqrt( sum((observations(:,1:2) - optimal_bundles).^2, 2));
normalized_euclid_criterions = (abs(observations(:,1) - optimal_bundles(:,1))) .* observations(:,3);

end
   
    
