function is_there_path = is_there_a_path(G, s, d)

% this function only checks whether there is a path from s to d in graph G
% it return 'true' if there is, 'false' if there isn't.
% it uses BFS (Breadth First Search) approach.

% improves performance, and also necessary sometimes to avoid errors (due to the use of logical indexing) 
G = logical(G);
% 
% % making sure there are no edges from a vertex to itself
% [n, ~] = size(G);
% for i=1:n
%     G(i,i) = false;
% end

while ~G(s,d)
    % s points to these vertices, i.e. there are edges e = (s,v) for each of these vertices 
    s_points_to = G(s,:); %logical(G(s,:));
    % if there is no place left to go, and yet we haven't found
    if ~any(s_points_to)
        is_there_path = false;
        return
    end
    % we don't want to visit these vertices again, so we remove their incoming edges 
    G(:,s_points_to) = false;
    % we add to s all the incoming edges of these vertices
    G(s,:) = any(G(s_points_to,:), 1);
end

% if we got here then G(s,d)==true, then we found a path from s to d
is_there_path = true;

end