function [SumRemovals , Removals] = HoutmanMaksCFGK_2G(Data)

% function summary = HoutmanMaksCFGK(Data)
% Input:
% (1) Data set (format described in the Data section)
% Output vector:
% (1) Minimum number of removals needed for acyclicality
% (2) Processing time

% Version: 3g
% Date: May 15, 2010
% Changes since previous version: More error handling added

% tic; % Start timer
% format compact; % Display format

%--------------------------------------------------------------------------
% Data

% Format for 'Data'
% [BundleLess BundleMore]
% BundleLess = Less preferred bundle
% BundleMore = More preferred bundle
%   Note: Each row is called a 'relation'

% Ex. If the selections from the 2nd and 6th budget sets were inside of
% the 4th budget set, but were not selected from that budget set, then
% it would generate this data:
% 2 4
% 6 4

%--------------------------------------------------------------------------
% Prepare subject data 

% Data characteristics needed for subsequent functions:
  
% MaxItem is the biggest number of observation among the 
% observations that are directly revealed 
% preferred to another bundle or that another bundle is 
% directly revealed preferred to them.

MaxItem = max([Data(:,1);Data(:,2)]);

% The number of pairs that satisfy the relation R^0.

NumberRelations = size(Data,1);

% Record preferences in 0-1 matrix
% Note: For each cell, the column number is preferred to the row number
% Example 1: Cs(3,2)=0 if 2 is never directly revealed preferred to 3
% Example 2: Cs(3,2)=1 if 2>3 for R^0
% Note: Cs is a directed graph (digraph)

Cs = zeros(MaxItem,MaxItem); % Initialize matrix (use zeros or sparse)
for m=1:NumberRelations
    Cs(Data(m,1),Data(m,2)) = 1; % Binary switch
end

AreCycles = FloydWarshall(Cs); % Test for cycles of R in the graph

%--------------------------------------------------------------------------
% Next we calculate the set of cycles of length 3 under R^0.
% Then we approximate the number of observations needed to be deleted for
% all the cycles to be broken

% note: Rose (1958) shows that every cycle of observations that violates SARP (xRy & yRx), 
% must contain a pair of observations that violates WARP (wR^0z & zP^0w). 
% Therefore, cycles of length 3 under relation R or under relation R^0 are identical. 

SumRemovals = 0; % Initialize minimum number of removals required

if AreCycles > 0 % Proceed only if cycles are present

    Cycles = SizeTwoCycles(Cs); % Extract/list short cycles (cycles of length 3 according to R^0) 

    % Dependency matrix where every row represents a cycle. 
    % The number of columns is MaxItem. 
    % On row i the values in cell(i,j) and (i,k) equal 1 and the rest are zeros 
    % iff (x^kR^0x^j)&(x^jR^0x^k) (meaning [j,k,j] is a cycle).

    DepMatrix = DependencyCFGK(Cycles,MaxItem); 

    if isempty(DepMatrix)

        % if there are no cycles, no observation is removed(removals is a vector of zeros)

        Removals = zeros(size(DepMatrix,2),1);

    else

        % Determine minimum removals, according to the linear program approximation, 
        % (it is an upper bound since it necessarily eliminates all cycles of size 3)

        Removals = InternalSolver(DepMatrix); 

    end

    SumRemovals = sum(Removals); % Store summary data for this subject

    % -------------------------------------------------------------

    % in the HoutmanMaksCFGK original program, the algorithm doesn't end here, 
    % but performs another check for cycles in R, this time for longer cycles. 
    % However, from Rose (1985) we know it is unnecessary. 

    % Keeps=(1-Removals);
    % KeepsMatrix=Keeps*Keeps';
    % 
    % % Ds is a square matrix of zeros and ones that represents the binary relation 
    % % R^0 after removing the observations of the minimal removals to clear cycles. 
    % 
    % Ds=Cs.*KeepsMatrix;
    % 
    % AreCycles=FloydWarshall(Ds); % Test for cycles in the graph
    % 
    % clear Cycles DepMatrix Removals Keeps KeepsMatrix Ds
    % 
    % if AreCycles>0 % Proceed only if cycles are present
    % 
    % Cycles=LongerCycles(Cs); % Extract/list short cycles  
    % 
    % DepMatrix=DependencyCFGK(Cycles,MaxItem); % Dependency matrix
    % 
    % if isempty(DepMatrix)
    % Removals=zeros(size(DepMatrix,2),1);
    % else
    % Removals=InternalSolver(DepMatrix); % Determine minimum removals
    % SumRemovals=sum(Removals); % Store summary data for this subject
    % end
    % 
    % Keeps=(1-Removals);
    % KeepsMatrix=Keeps*Keeps';
    % Ds=Cs.*KeepsMatrix;
    % 
    % AreCycles=FloydWarshall(Ds); % Test for cycles in the graph
    % 
    % clear Cycles DepMatrix Removals Keeps KeepsMatrix Ds
    % 
    % if AreCycles>0 % Proceed only if cycles are present
    %     
    % Cycles=Johnson(Cs); % Use's Johnson's algorithm to extract/list cycles     
    % 
    % DepMatrix=DependencyCFGK(Cycles,MaxItem); % Create Dependency matrix
    % 
    % clear Cycles
    % 
    % Removals=InternalSolver(DepMatrix); % Determine minimum removals
    % SumRemovals=sum(Removals); % Store summary data for this subject
    % 
    % Keeps=(1-Removals);
    % KeepsMatrix=Keeps*Keeps';
    % Ds=Cs.*KeepsMatrix;
    % 
    % AreCycles=FloydWarshall(Ds); % Test for cycles in the graph
    % 
    % if AreCycles>0
    %     disp('Not a solution');
    % end
    % 
    % clear Cs DepMatrix Removals Keeps KeepsMatrix Ds
    % 
    % end % if AreCycles>0
    % end % if AreCycles>0
    % end % if AreCycles>0
    % 
    % %--------------------------------------------------------------------------
    % Additional

    % ProcessingTime=toc; % Store subject processing time

    % summary=[SumRemovals ProcessingTime];

    % Clear remaining subject variables
    clear Data ans Structure MaxItem NumberRelations m AreCycles
    % clear SumRemovals ProcessingTime 

end

end