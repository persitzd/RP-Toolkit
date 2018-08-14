function AreCycles=FloydWarshall(F)

% function AreCycles=FloydWarshall(F)
% Takes an F matrix of integers (where n>0 in position i,j indicates 
% that jBi) and calculates the transitive closure using the 
% Floyd-Warshall algorithm.

% Version: 3f
% Date: July 8, 2009

n=length(F);

for k=1:n
    for j=1:n
        for i=1:n
            if and(F(i,k)>0,F(k,j)>0)
                F(i,j)=1;
            end;
        end;
    end;
end;

AreCycles=0;
if max(diag(F))>0 % max(diag(F))==0 if there are no cycles
    AreCycles=1; 
end

clear F n k j i

end