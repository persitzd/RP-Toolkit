function [HM, HM_raw, HM_raw_pseudo_residuals] = HPZ_Houtman_Maks_Graph_Approach (GARP_couples, approach_type)

% this function is, in fact, an implementation of an algorithm of minimum vertex cover 
% to see that, you may rename the variables as follows:
%   HM_raw -> size_mVC (size of minimum vertex cover)
%   HM -> percent_mVC (size of minimum vertex cover / number of vertices)
%   HM_raw_residuals -> is_in_mVC (a vector that says for each vertex if it is in any mVC, or not) 
%   HM_raw_pseudo_residuals -> if it is 1, we can say for sure that the vertex 
%                              is in some mVC, but if it is 0, we are not sure 
%   GARP_couples -> G (the graph as a matrix; (i,j)=1 iff there is an edge between vertices i and j) 

% WHEN FIRST CALLED (not inside the recursion), 
% approach_type should be equal 0
% (unless you are sure the graph contains no leaves)



% --- Explanation about the Algorithm ---
% The algorithm only looks for SOME minimum vertex cover, and not necessarily  
% all such covers. This aspect is what allows us to use Algorithm 2.
% --- Algorithm 1 ---
%   Divide the graph to its connected components,
%   and find the minimum vertex cover for each (then sum them up)
% --- Algorithm 2 ---
%   If there is a leaf in the graph, add the vertex connected to the leaf
%   to the on-build vertex cover, then drop the leaf and its connected
%   vertex and all their edges from the graph. Continue this until there
%   are no leaves, then pass the remaining graph to another algorithm that
%   will solve the rest (Algorithm 3).
% --- Algorithm 3 ---
%   Find a vertex with maximum degree in the graph, denote it by v.
%   Solve the minimum vertex cover for 2 options:
%   option I  : v is in the cover - then we drop v and its edges and solve
%              the remaining graph.
%   option II : v is not in the cover - then all the vertices connected to
%               it must be in the cover, so we drop v and all its connected
%               vertices and solve the remaining graph.
%   Finally, compare the results of the 2 options, and take the one that
%   achieved a smaller vertex cover.
% --- Combining the three Algorithms ---
%   in each iteration, we first use Algorithm 1 to divide the problem to a
%   few smaller problems (if possible), then we call either Algorithm 2 (if
%   in the previous iteration we used Algorithm 3) or Algorithm 3 (if
%   in the previous iteration we used Algorithm 2)
% --- The Importance of each of the three Algorithms ---
%   Algorithm 3 can solve the problem on its own, but it is not always
%   efficient enough, and it is exponential at worst.
%   Algorithm 2 gives a huge boost if there are leaves.
%   Algorithm 1 allows faster computation since the complexity inside 
%   each iteration is polynomial (of degree 2 or 3) in the number of observations.



% number of observations
[obs_num , ~] = size(GARP_couples);

% if there is only 1 observation, we don't need to perform any calculations
if obs_num == 1
    HM = 0;
    HM_raw = 0;
    HM_raw_pseudo_residuals = 0;
    return
end



% dividing the observations to subsets
[subsets_of_observations, num_of_subsets] = HPZ_Indices_Problem_Distribute([], [], GARP_couples);


% initializations
HM_raw = 0;
HM_raw_pseudo_residuals = zeros(obs_num, 1);


for i=1:num_of_subsets

    subset_obs_num = length(subsets_of_observations{i});
    
    if subset_obs_num > 1

%         % create new DRP & SDRP matrices with only the observations of this subset 
%         subset_DRP = DRP(subsets_of_observations{i} , subsets_of_observations{i});
%         subset_SDRP = SDRP(subsets_of_observations{i} , subsets_of_observations{i});
% 
%         subset_GARP = subset_DRP .* subset_SDRP';
%         subset_GARP_couples = max(GARP, GARP');

        % create a new GARP_couples matrix with only the observations of this subset 
        subset = subsets_of_observations{i};
        subset_GARP_couples = GARP_couples(subset , subset);
        
        
        if (approach_type == 0)

            % The 1st approach to ease the problem - breaking the graph apart using 
            % the fact that in vertex cover, it is always better to choose a vertex 
            % that is linked to a leaf, rather than the leaf itself

            % initializations
            is_handled = zeros(subset_obs_num, 1);
            obs_to_drop = zeros(subset_obs_num, 1);
            obs_to_keep = zeros(subset_obs_num, 1);
            previous_handled = -1;
            handled = 0;
            while (previous_handled ~= handled) && (handled ~= subset_obs_num)
                for j=1:subset_obs_num
                    %fprintf('obs: %d , relations = %d\n', biggest_subset(i), sum(new_GARP_couples(i,:)));
                    if sum(subset_GARP_couples(j,:)) == 1  && is_handled(j) == 0

                        connected_to = find(subset_GARP_couples(j,:));

                        is_handled(j) = 1;
                        is_handled(connected_to) = 1;

                        obs_to_drop(j) = 1;
                        obs_to_keep(connected_to) = 1;

                    elseif sum(subset_GARP_couples(j,:)) == 0

                        is_handled(j) = 1;
                        if obs_to_keep(j) ~= 1
                            obs_to_drop(j) = 1;
                        end
                    end
                end

                % delete edges from the graph that are connected to the vertices that were handled 
                for j=find(obs_to_drop)
                    subset_GARP_couples(j,:) = 0;
                    subset_GARP_couples(:,j) = 0;
                end
                for j=find(obs_to_keep)
                    subset_GARP_couples(j,:) = 0;
                    subset_GARP_couples(:,j) = 0;
                end

                % keep the previous handled before updating it
                previous_handled = handled;
                % update how many were handled
                handled = sum(is_handled);
            end

            raw_HM_counter = sum(obs_to_keep);
            reverse_raw_HM_counter = sum(obs_to_drop);
            if (raw_HM_counter + reverse_raw_HM_counter ~= handled)
                % this should never happen, in fact an error might be more fitting here than a warning 
                warning('counter = %d , reverse-counter = %d , handled = %d', raw_HM_counter, reverse_raw_HM_counter, handled);
            end


            % the remaining observations are all those that were not handled
            remaining_obs = find(~is_handled);

            % initializations
            %remaining_HM_raw = 0;
            %remaining_HM_raw_residuals = zeros(remaining_obs, 1);

            if ~isempty(remaining_obs)
                % relevant GARP is without transitivity when there are only 
                % 2 goods (thanks to the Rose (1958) Theorem)
                % calculate GARP without transitivity
                remaining_GARP_couples = subset_GARP_couples(remaining_obs , remaining_obs);
                %remaining_DRP = DRP(remaining_obs , remaining_obs);
                %remaining_SDRP = SDRP(remaining_obs , remaining_obs);

                % solving the remaining observations
                % (NOTE! it is curcial to use "~approach_type" in order
                % for this to work and not start an infinite loop)
                [~, remaining_HM_raw, remaining_HM_raw_residuals] = HPZ_Houtman_Maks_Graph_Approach (remaining_GARP_couples, ~approach_type);
            else
                remaining_HM_raw = 0;
                remaining_HM_raw_residuals = [];
            end

            % updating the full data values
            HM_raw = HM_raw + raw_HM_counter + remaining_HM_raw;
            HM_raw_pseudo_residuals(subset(find(obs_to_keep))) = 1; %#ok<*FNDSB>
            HM_raw_pseudo_residuals(subset(remaining_obs)) = remaining_HM_raw_residuals;

        else   % (approach_type == 1)

            % The 2nd approach to ease the problem - choosing the vertex with the
            % biggest degree, then referring to 2 different options:
            % either to discard it (meaning it was handled) or to discard all the
            % vertices that are connected to it (meaning both it and those vertices 
            % were handled)
            
            % the degree of each vertex (= observation)
            obs_degrees = sum(subset_GARP_couples);

            % the index of the vertex with biggest degree
            max_degree_obs = find((obs_degrees == max(obs_degrees)), 1);

            % indices of observation linked to this vertex
            linked_obs = find(subset_GARP_couples(max_degree_obs, :));

            % create a new matrix without this vertex / observation
            GARP_couples_1 = subset_GARP_couples;
            GARP_couples_1(max_degree_obs, :) = 0;
            GARP_couples_1(:, max_degree_obs) = 0;
            % run again on the new matrix
            [~, HM_raw_1, HM_raw_residuals_1] = HPZ_Houtman_Maks_Graph_Approach (GARP_couples_1, ~approach_type);
            % we also drop all this observations
            HM_raw_1 = HM_raw_1 + 1;
            
            % create a new matrix without these vertices / observations
            GARP_couples_2 = subset_GARP_couples;
            GARP_couples_2([linked_obs, max_degree_obs], :) = 0;
            GARP_couples_2(:, [linked_obs, max_degree_obs]) = 0;
            % run again on the new matrix
            [~, HM_raw_2, HM_raw_residuals_2] = HPZ_Houtman_Maks_Graph_Approach (GARP_couples_2, ~approach_type);
            % we also drop all the linked observations
            HM_raw_2 = HM_raw_2 + length(linked_obs);
            
            
            % we want to find the minimum number of observations to drop / minimum vertex cover 
            if HM_raw_1 < HM_raw_2
                HM_raw = HM_raw + HM_raw_1;
                %HM_exact = min(HM_exact, HM_exact_1);
                HM_raw_pseudo_residuals(subset(find(HM_raw_residuals_1))) = 1;
            elseif HM_raw_1 > HM_raw_2
                HM_raw = HM_raw + HM_raw_2;
                %HM_exact = min(HM_exact, HM_exact_2);
                HM_raw_pseudo_residuals(subset(find(HM_raw_residuals_2))) = 1;
            else % HM_raw_1 == HM_raw_2
                HM_raw = HM_raw + HM_raw_1;
                %HM_exact = min(HM_exact, max(HM_exact_1, HM_exact_2));
                HM_raw_pseudo_residuals(subset(find(HM_raw_residuals_1))) = 1;
                HM_raw_pseudo_residuals(subset(find(HM_raw_residuals_2))) = 1;
            end
        end

        
        % update the accumulated HM
        %HM_raw = HM_raw + subset_HM_raw;

        % update the residuals for these observations
        % if the subset residual is bigger than the subset HM, it means
        % that the raw residual was 0, and if it is smaller, it means that
        % the raw residual was 1
        %HM_raw_residuals(subset) = subset_HM_raw_pseudo_residuals;
    end

end



% calculating the HM index
HM = HM_raw / obs_num;



end
