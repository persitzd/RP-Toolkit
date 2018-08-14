function Removals = InternalSolver(DepMatrix)

% Removals=InternalSolver(DepMatrix)
% Determines a set of clicks of the smallest size needed to remove all
% cycles (i.e., make the graph acyclic).  'Removals' is a vector of ones
% and zeros where a one indicates that the click should be removed.

% Version: 3f
% Date: July 8, 2009

% Lanny altered this file on February 15th, 2016. I defined the number of
% oberservations as n (=DepClicks) and use this in the input for intlinprog
% below

[DepCycles , DepClicks] = size(DepMatrix); % Size of the problem
%disp(['Size of the problem: ' int2str(DepCycles)...
%    ' by ' int2str(DepClicks)]);

% Create constaints for bintprog
f = ones(DepClicks,1); % All clicks have equal wieght
b = ones(DepCycles,1); % All cycles must be removed at least once
% Reverse sign to reverse inequality
DepMatrix = -1 .* DepMatrix;
b = -1 .* b; 

if DepCycles == 1 % If only one removal is needed
    Removals = zeros(DepClicks,1);
    Removals(find(DepMatrix<0, 1 )) = 1; % Remove first click in the cycle
else
    % Optimization Toolbox solver
    % Removals=bintprog(f,DepMatrix,b);
    options = optimoptions('intlinprog','Display','off');
    Removals = intlinprog(f,1:DepClicks,DepMatrix,b,[],[],zeros(DepClicks,1),ones(DepClicks,1),options);    
    
    % GLPKMEX solver
 %   vartype=[];
 %   for i=1:DepClicks
 %       vartype=[vartype 'B'];
 %   end
 %   Removals=glpk(f,DepMatrix,b,[],[],[],vartype);
end

clear DepMatrix f b

end