function summary=HoutmanMaksCFGK(Data)

% function summary=HoutmanMaksCFGK(Data)
% Input:
% (1) Data set (format described in the Data section)
% Output vector:
% (1) Minimum number of removals needed for acyclicality
% (2) Processing time

% Version: 3g
% Date: May 15, 2010
% Changes since previous version: More error handling added

tic; % Start timer
format compact; % Display format

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

% Data characteristics needed for subsequent functions
MaxItem=max([Data(:,1);Data(:,2)]);
NumberRelations=size(Data,1);

% Record preferences in 0-1 matrix
% Note: For each cell, the column number is preferred to the row number
% Example 1: Cs(3,2)=0 if 2 is never revealed preferred to 3
% Example 2: Cs(3,2)=1 if 2>3 for some relation
% Note: Cs is a directed graph (digraph)
Cs=zeros(MaxItem,MaxItem); % Initialize matrix (use zeros or sparse)
for m=1:NumberRelations
    Cs(Data(m,1),Data(m,2))=1; % Binary switch
end

AreCycles=FloydWarshall(Cs); % Test for cycles in the graph

%--------------------------------------------------------------------------
% Next we calculate the set of cycles of length 3 or less.
% Then we approximate the number of observations needed to be deleted for
% all the cycles to be broken

SumRemovals=0; % Initialize minimum number of removals required

if AreCycles>0 % Proceed only if cycles are present

Cycles=ShortCycles(Cs); % Extract/list short cycles (cycles of length 3 or less) 

DepMatrix=DependencyCFGK(Cycles,MaxItem); % Dependency matrix

if isempty(DepMatrix)
Removals=zeros(size(DepMatrix,2),1);
else
Removals=InternalSolver(DepMatrix); % Determine minimum removals
SumRemovals=sum(Removals); % Store summary data for this subject
end

% -------------------------------------------------------------

Keeps=(1-Removals);
KeepsMatrix=Keeps*Keeps';
Ds=Cs.*KeepsMatrix;

AreCycles=FloydWarshall(Ds); % Test for cycles in the graph

clear Cycles DepMatrix Removals Keeps KeepsMatrix Ds

if AreCycles>0 % Proceed only if cycles are present

Cycles=LongerCycles(Cs); % Extract/list short cycles  

DepMatrix=DependencyCFGK(Cycles,MaxItem); % Dependency matrix

if isempty(DepMatrix)
Removals=zeros(size(DepMatrix,2),1);
else
Removals=InternalSolver(DepMatrix); % Determine minimum removals
SumRemovals=sum(Removals); % Store summary data for this subject
end

Keeps=(1-Removals);
KeepsMatrix=Keeps*Keeps';
Ds=Cs.*KeepsMatrix;

AreCycles=FloydWarshall(Ds); % Test for cycles in the graph

clear Cycles DepMatrix Removals Keeps KeepsMatrix Ds

if AreCycles>0 % Proceed only if cycles are present
    
Cycles=Johnson(Cs); % Use's Johnson's algorithm to extract/list cycles     

DepMatrix=DependencyCFGK(Cycles,MaxItem); % Create Dependency matrix

clear Cycles

Removals=InternalSolver(DepMatrix); % Determine minimum removals
SumRemovals=sum(Removals); % Store summary data for this subject

Keeps=(1-Removals);
KeepsMatrix=Keeps*Keeps';
Ds=Cs.*KeepsMatrix;

AreCycles=FloydWarshall(Ds); % Test for cycles in the graph

if AreCycles>0
    disp('Not a solution');
end

clear Cs DepMatrix Removals Keeps KeepsMatrix Ds

end % if AreCycles>0
end % if AreCycles>0
end % if AreCycles>0

%--------------------------------------------------------------------------
% Additional

ProcessingTime=toc; % Store subject processing time

summary=[SumRemovals ProcessingTime];

% Clear remaining subject variables
clear Data ans Structure MaxItem NumberRelations m AreCycles
clear SumRemovals ProcessingTime 

end