function [utility] = CARA (x,A)

% The function calculates the value of the CARA function with parameter A
% and the variable x.

utility = -exp(-(A*x));

end