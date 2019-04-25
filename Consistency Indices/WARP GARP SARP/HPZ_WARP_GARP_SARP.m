function [FLAGS, VIO_PAIRS, VIOLATIONS, WARP, SARP] = HPZ_WARP_GARP_SARP(identical_choice, DRP, RP, GARP)

% To test for SARP we will use the definition of SARP1. For every pair of
% choices x and y if xRy and yRx it must be that x=y. We will take RP and
% its transpose and multiply element by element. Every 1 corresponds to 
% xRy and yRx. Then we take off the identity matrix (all the pairs x=y). 
% The final matrix is the zero matrix if and only if SARP is satisfied. 


% The following conforms with Varian (1982) definition of SARP1 and SARP2
% SARP3 would be SARP = RP.*(DRP') - identical_choice;
% As a binary test these definitions are equivalent. However, SARP1 and
% SARP2 may report more violations.

% RP(i,j)*RP'(j,i)=1 iff (x^i R x^j)&(x^j R x^i), otherwise it equals 0. For every i,j identical
% choice SARP(i,j)=0 (including i=j of course) 

% OLD CODE:
%SARP = RP.*(RP') - identical_choice;
% NEW CODE:
RP_no_identical = min(RP, ~identical_choice);
SARP = RP_no_identical .* (RP_no_identical');

% A flag that keeps the SARP result. It is 0 if and only if the data
% satisfies SARP.
SARP_FLAG = 1;

SARP_ERRORS = sum(sum(SARP));   % sums the values of all matrix SARP's cells, counting 
                                % the number of violations (notice that every pair of 
                                % bundles x^i,x^j is counted as two violations; one for 
                                % (i,j), and the other for (j,i)). 


% sums the values of all matrix SARP's cells, counting the number of violations 
% pairs (every two pair of bundles x^i,x^j is counted as one violation pair). 
SARP_VIO_PAIRS = sum(sum(triu(SARP|(SARP'))));

% If the SARP matrix is only zeros then the data satisfies SARP
if SARP_ERRORS == 0
  SARP_FLAG = 0;
end 



% A flag that keeps the GARP result. It is 0 if and only if the data
% satisfies GARP.
GARP_FLAG = 1;

% GARP_ERRORS sums the values of all matrix GARP's cells. Counting 
% the number of violations (notice that a pair of bundles x^i,x^j might be counted 
% as two violations; one for (i,j), and the other for (j,i)). 
GARP_ERRORS = sum(sum(GARP));

% GARP_VIO_PAIRS sums the values of all matrix
% GARP's cells, counting the number of violations pairs (every two pair of bundles 
% x^i,x^j is counted as one violation pair). 
GARP_VIO_PAIRS = sum(sum(triu(GARP|(GARP'))));


% If the GARP matrix is only zeros then the data satisfies GARP
if GARP_ERRORS == 0
  GARP_FLAG = 0;
end 



% To test for WARP we will do the following: for every pair of
% choices x and y if xR0y then not yR0x. We will take DRP and the transpose of
% DRP and multiply element by element. Every 1 corresponds to 
% xR0y and yR0x. BUT we fix it so identical choices will not be counted, as
% they do not violate WARP.
% The final matrix is the zero matrix if and only if 
% WARP is satisfied. 

% OLD CODE:
%WARP = DRP.*(DRP') - identical_choice;
% NEW CODE:
DRP_no_identical = min(DRP, ~identical_choice);
WARP = DRP_no_identical .* (DRP_no_identical');

% A flag that keeps the WARP result. It is 0 if and only if the data
% satisfies WARP.
WARP_FLAG = 1;

% WARP_ERRORS sums the values of all matrix WARP's cells. Counting 
% the number of violations (notice that a pair of bundles x^i,x^j might be counted 
% as two violations; one for (i,j), and the other for (j,i)). 
WARP_ERRORS = sum(sum(WARP));

% WARP_VIO_PAIRS sums the values of all matrix
% GARP's cells, counting the number of violations pairs (every two pair of bundles 
% x^i,x^j is counted as one violation pair). 
WARP_VIO_PAIRS = sum(sum(triu(WARP|(WARP'))));

% If the WARP matrix is only zeros then the data satisfies WARP
if WARP_ERRORS == 0
  WARP_FLAG = 0;
end 



% assigning the results to the results vectors
FLAGS(1) = WARP_FLAG;
FLAGS(2) = GARP_FLAG;
FLAGS(3) = SARP_FLAG;

VIO_PAIRS(1) = WARP_VIO_PAIRS;
VIO_PAIRS(2) = GARP_VIO_PAIRS;
VIO_PAIRS(3) = SARP_VIO_PAIRS;

VIOLATIONS(1) = WARP_ERRORS;
VIOLATIONS(2) = GARP_ERRORS;
VIOLATIONS(3) = SARP_ERRORS;

end