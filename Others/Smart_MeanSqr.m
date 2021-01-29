function meansqr = Smart_MeanSqr(vec)

% This meansqr function takes into acount the sign of the values.

meansqr = sum(sign(vec) .* (vec.^2));

end