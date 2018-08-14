function DepMatrix=DependencyCFGK(Cycles,MaxItem)

% function DepMatrix=DependencyCFGK(Cycles,MaxItem)
% Create Dependency matrix.

% Version: 3f
% Date: July 8, 2009

n=1;

% Start Dependency matrix (each click is a column)
CyclesLength=length(Cycles); % Minimum number of cycles
DepMatrix=zeros(CyclesLength,MaxItem); % Zeros or sparse
for p=1:CyclesLength % For each cycle
    TempCycle=Cycles{p};
    TempCycleLength=size(TempCycle,2); 
        for r=1:(TempCycleLength-1)
            DepMatrix(n,TempCycle(r))=1;
        end % for r
        n=n+1; % End of cycle, so advance to next row of the matrix
end % for p   

clear Cycles MaxItem
clear n CyclesLength
clear p TempCycle TempCycleLength r

end