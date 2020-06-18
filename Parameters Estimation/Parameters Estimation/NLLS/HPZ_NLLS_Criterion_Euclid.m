function [euclid_criterions] = HPZ_NLLS_Criterion_Euclid(observations, optimal_bundles)

% The function calculates the Euclidean NLLS criterion for 
% a given set of observations and predictions.
% The function returns the value of the criterion for the given data. 

euclid_criterions =  sqrt( sum((observations(:,1:2) - optimal_bundles).^2, 2));

end
   
    
