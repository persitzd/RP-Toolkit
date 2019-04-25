function min_cycles_with_path = find_all_minimal_cycles_with_path (G, Path)

    % this function is a Matlab implementation of the NEXT STEP
    % algorithm from the paper (as of 25.04.19 - working paper): 
    % "Code Package for Revealed Preferences Analysis of Choices from Linear Budget Sets"  
    %
    % as we mention in one of the footnotes in the article, we also perform
    % BFS (Breadth First Search) in each iteration to check if there is
    % still some path that leads from the current node to the root.
    
    
    
    % initialization
    min_cycles_with_path = {};
    
    % Initial values: d is the last element in P; r is the first element in P.  
    %
    % we have a path starting from 'Root' and ending with 'Current',
    % and we want to complete the path to a cycle by finding a path
    % from current to root
    Root = Path(1);
    Current = Path(end);

    % 1. Use G to calculate OutG(d).  
    OutCurrent_Indexes = G(Current,:);
    OutCurrent = find(OutCurrent_Indexes);   % OutCurrent contains all the vertices that Current points to
    
    % 2. If OutG(d) is empty, terminate.
    if isempty(OutCurrent)
       return 
    end

    % 3. Else if r in OutG(d), Add P to MIN CYCLES and terminate.
    %
    % if Current points directly to Root, then we managed to close a
    % cycle, also there is no need to continue to the other vertices in
    % OutCurrent as any cycle they may lead to, this cycle we now found
    % will be its subcycle
    if G(Current, Root)
        min_cycles_with_path = {Path};
        return
    end

    
    
    % ===== NOT PART OF THE ORIGINAL ALGORITHM - BUT IMPROVES RUNNING TIME ===== 
    % for improved time complexity, we check if the current subgraph of
    % the original G still contains some path from Current to Root, if
    % it does not - we immediately return an empty list
    if Current ~= Root
        is_there_path = is_there_a_path(G, Current, Root);
        if ~is_there_path
           return 
        end
        % we prefer to keep Current in the matrix and not delete it,
        % but in order to "block" it, we turn all its rows and columns to 0. 
        G(Current,:) = 0;
        G(:,Current) = 0;
    end
    
    
    
    % Now we go through each of Current's neighbors, to see if any of
    % them can lead the path to a minimal cycle
    for Next = OutCurrent

        % 'Next' stands for one of the vertices in OutCurrent
        
        % 4(c). set Others to OutG(d)\{d'}.
        %
        % OutCurrentOthers contains all OutCurrent except for the current Next 
        OutCurrentOthers_Indexes = OutCurrent_Indexes;
        OutCurrentOthers_Indexes(Next) = 0;

        % 4(d). set Pointers to InG(d')\{r}
        %
        % InNext contains all the vertices that point to Next, 
        % except for Current, that was handled earlier in a different manner  
        InNext_Indexes = G(:,Next)';
        InNext_Indexes(Current) = 0;

        % 4(e). set G' to G without the columns and rows corresponding to Others and pointers  
        %
        % for this Next, we create G_Next as a subgraph of G by removing
        % from G all vertices in OutCurrentOthers and in InNext_Indexes  
        Verices_To_Drop_Indexes = OutCurrentOthers_Indexes | InNext_Indexes;
        Verices_To_Remain_Indexes = ~Verices_To_Drop_Indexes;
        Verices_To_Remain = find(Verices_To_Remain_Indexes);
        G_Next = G(Verices_To_Remain , Verices_To_Remain);

        % We need to translate the indexes of Path and of Next to fit
        % the indexes of G_Next which is a truncated version of G
        New_Indexes_Of_Verices_To_Remain = cumsum(Verices_To_Remain_Indexes);
        Path_Next = New_Indexes_Of_Verices_To_Remain([Path,Next]);

        % 4(f). call NEXT STEP with parameters G' and [P,d'].
        %
        % Now we perform a recursive call with G_sub and with new-path = [Path,Next] 
        min_cycles_with_path_next = find_all_minimal_cycles_with_path (G_Next, Path_Next);

        for c = 1:length(min_cycles_with_path_next)
            % the indexes that were returned are different, since we
            % dropped vertices from G to create G_sub, so now we need
            % to translate the indexes back to their original form
            min_cycles_with_path_next{c} = Verices_To_Remain(min_cycles_with_path_next{c});
        end

        % update the minimal cycles list
        min_cycles_with_path = [min_cycles_with_path , min_cycles_with_path_next]; %#ok<AGROW>
    end
    
    
end