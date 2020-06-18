function [ldr_criterions] = HPZ_NLLS_Criterion_Ldr(observations, optimal_bundles)

% This function assumes the choice is on the budget
% line, and never beneath it. Therefore, it must not be used for
% experiments where the DM can choose and did choose bundle that are not on
% the budget line, or at least very close to it. If the DM made achose not
% on the budget line, this function will not give any indication for that.

% The function calculates the NLLS criterion (similar to CFGK (2007)) for 
% a given set of observations and predictions.
% The function returns the value of the criterion for the specified
% functional form and the given data. 

ldr_chosen = log(observations(:,2)./observations(:,1));

ldr_predicted = log(optimal_bundles(:,2)./optimal_bundles(:,1));

ldr_criterions = (ldr_chosen - ldr_predicted).^2;

end
   
    
