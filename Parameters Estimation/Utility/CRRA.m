function [utility] = CRRA (x, rho)

% The function calculates the value of the CRRA function with parameter rho
% and the variable x.

% if rho = 1 then the CRRA utility is going to be infinity, in order to
% avoid this problem we use the logaritmic form:
%if abs(rho-1) < 10^(-6)
if (rho == 1)
    utility = log(x);
else
    utility = (x^(1-rho)) / (1-rho);
end

end