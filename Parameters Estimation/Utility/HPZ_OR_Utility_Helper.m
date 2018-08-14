function [minus_utility] = HPZ_OR_Utility_Helper(x, param1, param2, function_flag, varargin)

% This function's goal is to adjust the result from HPZ_OR_Utility, so it will
% work well when used in fmincon.
% The adjustments are:
% 1. Turning the utility to be minus instead of plus (because fmincon finds
%   the minimum, while we want to find the maximum utility)
% 2. Doubling the utility by some multiplier in order to make the utility
%   values and/or derivatives of the utility function big enough so that
%   the algorithm won't easily be terminated by reaching the thresholds of fmincon.

% x - a vector of quantities - x(1) is the quantity of x in the given 
%   bundle while x(2) is the quantity of y in the given bundle.
% param1 is alpha - the share parameter.
% param2 is rho - the parameter of the CES function.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% alpha_multiplier is needed in order to avoid a case that the optimization 
% is terminated too soon, as a result of reaching some defined threshold. 
% alpha_multiplier increases the values being tested, hence making it 
% harder for the threshold to be reached.
% (NOTE: there is no connection between alpha the share parameter,
%  and this "alpha_multiplier")
if (length(varargin) == 1)
    alpha_multiplier = cell2mat(varargin(1));
else
    alpha_multiplier = 1;
end



% the real utility
minus_utility_base = (0 - HPZ_OR_Utility(x, param1, param2, function_flag));



% the modified utility
minus_utility = alpha_multiplier * minus_utility_base;


end