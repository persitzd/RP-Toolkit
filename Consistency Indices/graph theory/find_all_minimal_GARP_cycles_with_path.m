function min_cycles_with_path = find_all_minimal_GARP_cycles_with_path (DRP, SDRP, Path, varargin)
    
    % when calling this function from another function, please leave 
    % varargin empty, otherwise the function will epically fail.
    
    % this function is a Matlab implementation of the ADAPTED NEXT STEP
    % algorithm from the paper (as of 25.04.19 - working paper): 
    % "Code Package for Revealed Preferences Analysis of Choices from Linear Budget Sets"  
    %
    % as we mention in one of the footnotes in the article, we also perform
    % BFS (Breadth First Search) in each iteration to check if there is
    % still some path that leads from the current node to the root.
    
    
    
    % initialization
    min_cycles_with_path = {};

    % we have a path starting from 'Root' and ending with 'Current',
    % and we want to complete the path to a cycle by finding a path
    % from current to root
    Root = Path(1);
    Current = Path(end);
    
    
    
    if isempty(varargin)
        
        % then it is the initial call to this function, we need to initialize some variables:  
        
        SDRP_State = 0;   % BLACK STATE = 0 in the initial call
        
        
        % in the initial call, Leftover Others and Potential Pointers are empty:  
        
        % list of Out(v) (for some v in the path) indexes according to DRP 
        % but not according to SDRP, hence we didn't drop them. Initially set to none.  
        Leftover_Others_Indexes = false(1, length(DRP));
        % list of In(v) (for some v in the path) indexes according to DRP 
        % but not according to SDRP, hence we didn't drop them. Initially set to none.  
        Potential_Pointers_Indexes = false(1, length(DRP));
        
    else
        
        % if it is not the initial call but it is a recursive call, then it
        % must have passed some important arguments:
        
        SDRP_State = varargin{1};   % BLACK STATE is the NEXT BLACK STATE that was passed from the previous recursive call  
        
        
        % Leftover Others and Potential Pointers are the NEXT Leftover Others 
        % and NEXT Potential Pointers were passed from the previous recursive call: 
        
        % list of Out(v) (for some v in the path) indexes according to DRP 
        % but not according to SDRP, hence we didn't drop them.
        Leftover_Others_Indexes = varargin{2};
        % list of In(v) (for some v in the path) indexes according to DRP 
        % but not according to SDRP, hence we didn't drop them. Initially set to none.  
        Potential_Pointers_Indexes = varargin{3};
        
    end
    
    
    
    % notations:
    % r  <=>  Root
    % d  <=>  Current
    % d' <=>  Next
    % SDRP_State <=> BLACK STATE
    
    
    
    % 2. If any of the following conditions is satisfied, terminate:
    %   (a) OutWhite(d) is empty.
    %   (b) Black State = 2 and r in OutBlack(d).
    %	(c) Black State = 2 and OutDiff(d) is empty.
    if ~any(DRP(Current,:)) || ...
        ( SDRP_State == 2 && SDRP(Current, Root) ) || ...
        ( SDRP_State == 2 && ~any( DRP(Current,:) & ~SDRP(Current,:) ) )
        
        return
    end
    
    % 3. If any of the following conditions is satisfied, 
    %    add P to MIN CYCLES and terminate:
    %   (a) r in OutBlack(d) and Black State < 2
    %   (b) r in OutDiff(d) and Black State > 0
    if ( SDRP(Current, Root) && SDRP_State < 2 ) || ...
       ( DRP(Current, Root) && ~SDRP(Current, Root) && SDRP_State > 0 )
        
        min_cycles_with_path = [min_cycles_with_path , Path];
        return
    end
    
    
    
    % ===== NOT PART OF THE ORIGINAL ALGORITHM - BUT IMPROVES RUNNING TIME ===== 
    % for improved time complexity, we check if the current subgraph of
    % the original DRP still contains some path from Current to Root, if
    % it does not - we immediately return an empty list
    if Current ~= Root
        % if BLACK STATE is not two, we check if there is a path in WHITE,
        % whereas if it is two, we check if there is a path in WHITE \ BLACK  
        if SDRP_State < 2
            Graph_to_check = DRP;
        else
            Graph_to_check = DRP & ~SDRP;
        end
        % checking if there is a path using BFS
        is_there_path = is_there_a_path(Graph_to_check, Current, Root);
        if ~is_there_path
            % if there isn't a path, there is no point in continuing 
            return 
        end
    end
    
    
    
    % 4. loop over members of OutWhite(d)
    OutCurrent = find(DRP(Current,:));
    for Next = OutCurrent
        
        % 4(b)ii. If any of the following conditions is satisfied, Go to 4a:  
        %   A. d' = r
        %   B. Black State = 0 and d' not-in OutBlack(d) and d' in Leftover Others   
        %   C. Black State = 2 and d' in OutBlack(d)
        %   D. d' in OutBlack(d) and d' in Potential Pointers
        if ( Next == Root ) || ...
           ( SDRP_State == 0 && ~SDRP(Current, Next) && Leftover_Others_Indexes(Next) ) || ...
           ( SDRP_State == 2 && SDRP(Current, Next) ) || ...
           ( SDRP(Current, Next) && Potential_Pointers_Indexes(Next) )
            
            continue
        end
        
        % 4(b)iii. 
        % If d' in Leftover Others, 
        %       set Next Black State = 2.
        % Elseif d' in OutBlack(d), 
        %       set Next Black State = 1.  
        % Else, 
        %       set Next Black State = Black State.
        if Leftover_Others_Indexes(Next)
            Next_SDRP_State = 2;
        elseif SDRP(Current, Next)
            Next_SDRP_State = 1;
        else
            Next_SDRP_State = SDRP_State;
        end
        
        % 4(c). 
        % If Black State > 0, 
        %       set Others to OutWhite(d)\{d'} and set
        %       Next Leftover Others to Leftover Others.
        % Else, 
        %       set Others to OutBlack(d)\{d'} and set Next Leftover Others 
        %       to the union of Leftover Others and OutDiff(d).
        if SDRP_State > 0
            Others_Indexes = DRP(Current,:);
            Others_Indexes(Next) = false;
            Next_Leftover_Others_Indexes = Leftover_Others_Indexes;
        else
            Others_Indexes = SDRP(Current,:);
            Others_Indexes(Next) = false;
            Next_Leftover_Others_Indexes = Leftover_Others_Indexes | (DRP(Current,:) & ~SDRP(Current,:));
        end
        
        % 4(d)i. Set Pointers to OutBlack(d')/{r}.
        Pointers_Indexes = SDRP(:,Next)';
        Pointers_Indexes(Root) = 0;
        
        % 4(d)ii. 
        % If d' in OutBlack(d), 
        %       add Potential Pointers to Pointers
        %       and set Next Potential Pointers to InDiff(d')\{r}.
        % Else, 
        %       set Next Potential Pointers to the union of 
        %       Potential Pointers and InDiff(d')\{r}.
        if SDRP(Current, Next)
            Pointers_Indexes = Pointers_Indexes | Potential_Pointers_Indexes;
            Next_Potential_Pointers_Indexes = (DRP(:,Next)' & ~SDRP(:,Next)');
            Next_Potential_Pointers_Indexes(Root) = false;
        else
            Next_Potential_Pointers_Indexes = Potential_Pointers_Indexes | (DRP(:,Next)' & ~SDRP(:,Next)');
            Next_Potential_Pointers_Indexes(Root) = false;
        end
        
        % 4(e). Set Next White and Next Black to White and Black, respectively, 
        % without the columns and rows corresponding to Others, Pointers and to d (if d ~= r). 
        Next_DRP = DRP;
        Next_SDRP = SDRP;
        Obs_to_Remove = Others_Indexes | Pointers_Indexes;
        Obs_to_Remove(Current) = true;
        Obs_to_Remove(Root) = false;
        Next_DRP(Obs_to_Remove , :) = false;
        Next_DRP(: , Obs_to_Remove) = false;
        Next_SDRP(Obs_to_Remove , :) = false;
        Next_SDRP(: , Obs_to_Remove) = false;
        
        % 4(f). Call ADAPTED NEXT STEP with parameters Next White, Next
        % Black, [P,d'], Next Black State, Next Leftover Others, Next Potential Pointers,
        % (perform a recursive call with the extended path)
        min_cycles_with_path_next = find_all_minimal_GARP_cycles_with_path (Next_DRP, Next_SDRP, [Path, Next], ...
                                                            Next_SDRP_State, Next_Leftover_Others_Indexes, Next_Potential_Pointers_Indexes);
        
        % (update the minimal cycles list)
        min_cycles_with_path = [min_cycles_with_path , min_cycles_with_path_next]; %#ok<AGROW>
    end
    

end