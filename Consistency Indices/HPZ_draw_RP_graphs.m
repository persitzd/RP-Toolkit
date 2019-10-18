function HPZ_draw_RP_graphs (Graph_flags, main_folder_for_results, subject_num, expenditure, identical_choice, index_threshold, bundles_chosen)

% This function creates and saves (as PNG) graphs that represent the
% revealed preference relation of a given subject.
% Since drawing all the relations results graphs full of edges and
% unreadable, we reduce the number of edges in one of 3 ways:
%   1. (if subject is consistent) - we delete any direct relation, that can
%      be derived by transitivity from other relations. By doing so, we
%      create a hierarchy graph of the revealed preference relation.
%   2. (if subject is isconsistent) - we delete all edges (relations),
%      except for those involved in any GARP-violating cycles.
%      this graph allows to investigate the violations from consistnecy.
%   3. (if subject is isconsistent) - we delete all edges (relations),
%      except for those involved in any GARP-violating *minimal* cycles.
%      this graph allows to investigate the basic violations from consistnecy. 

obs_num = length(expenditure);

% finding GARP, and the relations DRP and SDRP
[GARP, ~, ~, DRP, SDRP] = GARP_based_on_expenditures(expenditure, identical_choice, index_threshold);

% a more "relaxed" identical choice matrix (treats bundles as identical even if they are slightly off)
identical_choice_relaxed = zeros(obs_num, obs_num);   % Value of cell (j,k) is 1 if choices identical, and 0 otherwise.
% the threshold we will use here is:
identical_choices_relaxed_threshold = 0.001;   % a difference of no more than 0.1% in each of the goods
for j=1:obs_num
    % going through all identical_choice’s cells.
    for k=1:obs_num
        % if the choices j&k are identical, then set value of cell (j,k) = 1 , otherwise 0. 
        if all( abs(bundles_chosen(j, :) - bundles_chosen(k, :)) <= bundles_chosen(j, :)*identical_choices_relaxed_threshold )
        % if all( (abs(Choices(j,:) - Choices(k,:)) <= Choices(j,:)*identical_choices_threshold) )   %EXTENSION 
            identical_choice_relaxed(j,k) = 1;
        end
    end
end


%  create a graph of revealed preference relations, and print it (if the user specified so) 
RP_graph_file_name_base = strcat(main_folder_for_results, '/', HPZ_Constants.results_files_dir, '/', 'Subject-', num2str(subject_num), '-', num2str(Graph_flags(end)));
RP_graph_file_name_extra = {'Hierarchy Bundles', 'Hierarchy Observations', 'All Cycles', 'Minimal Cycles'};
for graph_type=1:4
    % whether to draw - depends on whether the user specified, and if consistent or not 
    % (layout can be: 'auto' | 'circle' | 'force' | 'layered' | 'subspace' | 'force3' | 'subspace3')   
    if graph_type==1
        to_draw = (Graph_flags(graph_type) && sum(sum(GARP)) == 0);
        layout = 'layered';
    elseif graph_type==2
        to_draw = (Graph_flags(graph_type) && sum(sum(GARP)) == 0);
        layout = 'layered';
    elseif graph_type==3 || graph_type==4
        to_draw = (Graph_flags(graph_type) && sum(sum(GARP)) > 0);
        layout = 'force';
        % when printing weights, they are unreadable unless we use circle layout  
        if Graph_flags(5)==1
            layout = 'circle';
        end
    end
    
    % now to the real thing:
    if to_draw
        
        % create the matrix that represents the graph (depends on graph type) 
        if graph_type==1
            % Hierarchy Graph with names = bundles (for consistent subjects)
            
            Graph_Matrix = DRP + SDRP;
            
            % this part is needed in order to turn observations numbers to bundles  
            if true
                % which obs to drop from Graph_Matrix
                obs_to_drop = false(1, obs_num);
                % sets of identical bundles
                sets_of_identical_bundles = cell(0);
                % this vector keeps track of those observations that were already dealt with  
                was_found_identical = false(1, obs_num);
                for obs = 1:obs_num
                    if was_found_identical(obs)
                        % we already dealt with this observation,
                        % which means it was identical to a previous one.  
                        obs_to_drop(obs) = true;
                        continue
                    end
                    % find which observations are identical to the current observation  
                    identical_to_obs = find(identical_choice_relaxed(obs, :));
                    % save this set of observations with identical chosen bundle    
                    sets_of_identical_bundles{end+1} = identical_to_obs; %#ok<AGROW>
                    % remember that we dealt with these observations
                    was_found_identical(identical_to_obs) = true;
                end
                % number of non-identical chosen bundles
                num_of_non_identical_bundles = length(sets_of_identical_bundles);
                % initialization of matrix that will hold all different (non-identical) bundles 
                bundles_chosen_merged_identicals = nan(num_of_non_identical_bundles, size(bundles_chosen,2));
                % now we assign values and update Graph_Matrix  
                for identical_obs_set = 1:num_of_non_identical_bundles
                    current_set = sets_of_identical_bundles{identical_obs_set};
                    first_obs = current_set(1);   % only this obs will remain from this set
                    % combine the direct revealed preference relation of all observations with this identical chosen bundle 
                    Graph_Matrix(first_obs, :) = max(Graph_Matrix(current_set, :), [], 1); %max(Graph_Matrix(first_obs, :) , max(Graph_Matrix(current_set, :), 1));
                    Graph_Matrix(:, first_obs) = max(Graph_Matrix(:, current_set), [], 2); %max(Graph_Matrix(:, first_obs) , max(Graph_Matrix(:, current_set), 2));
                    % the bundle that was chosen in each observation in this set   
                    bundles_chosen_merged_identicals(identical_obs_set , :) = bundles_chosen(first_obs , :);
                end
                % remove all unnecessary observations
                Graph_Matrix = Graph_Matrix(~obs_to_drop , ~obs_to_drop);
                % % update the number of observations that will be in the graph 
                % obs_num = num_of_non_identical_bundles;

            end   
            
            % delete edges that are not necessary to display the hierarchy 
            for i=1:num_of_non_identical_bundles
                for j=[1:(i-1) , (i+1):num_of_non_identical_bundles]
                    if Graph_Matrix(i,j)
                        Graph_Matrix_temp = Graph_Matrix;
                        Graph_Matrix_temp(i,j) = 0;
                        is_there_path = is_there_a_path(Graph_Matrix_temp, i, j);
                        %[cost, ~] = dijkstra_edge_count(Graph_Matrix_temp, i, j);
                        if is_there_path %~isinf(cost)
                            Graph_Matrix(i,j) = 0;
                        end
                    end
                end
            end
            
        elseif graph_type==2
            % Hierarchy Graph names = observations indexes (for consistent subjects)
            
            if Graph_flags(5)==0
                Graph_Matrix = DRP + SDRP;  
            elseif Graph_flags(5)==1
                Graph_Matrix = DRP .* expenditure;
            end  
            % delete edges that are not necessary to display the hierarchy 
            for i=1:obs_num
                for j=[1:(i-1) , (i+1):obs_num]
                    if Graph_Matrix(i,j)
                        Graph_Matrix_temp = Graph_Matrix;
                        Graph_Matrix_temp(i,j) = 0;
                        is_there_path = is_there_a_path(Graph_Matrix_temp, i, j);
                        %[cost, ~] = dijkstra_edge_count(Graph_Matrix_temp, i, j);
                        if is_there_path %~isinf(cost)
                            Graph_Matrix(i,j) = 0;
                        end
                    end
                end
            end
            
        elseif graph_type==3
            % All Cycles Graph
            [relevant_DRP, relevant_SDRP] = find_relevant_relations(DRP, SDRP);
            if Graph_flags(5)==0
                Graph_Matrix = relevant_DRP + relevant_SDRP;
            elseif Graph_flags(5)==1
                Graph_Matrix = relevant_DRP .* expenditure;
            end   
            
        elseif graph_type==4
            % Minimal Cycles Graph
            minimal_cycles = find_all_minimal_GARP_cycles (DRP, SDRP);
            Graph_Matrix = zeros(obs_num,obs_num);
            for c=1:length(minimal_cycles)
                current_cycle = minimal_cycles{c};
                for v=1:(length(current_cycle)-1)
                    Graph_Matrix(current_cycle(v) , current_cycle(v+1)) = 1;
                end
                Graph_Matrix(current_cycle(end) , current_cycle(1)) = 1;
            end
            if Graph_flags(5)==0
                Graph_Matrix = Graph_Matrix .* (DRP + SDRP);
            elseif Graph_flags(5)==1
                Graph_Matrix = Graph_Matrix .* expenditure;
            end 
        end
        
        
        % making sure there are no edges from a vertex to itself
        for i=1:length(Graph_Matrix)
            Graph_Matrix(i,i) = 0;
        end
        % creating the graph
        G = digraph(Graph_Matrix);
        % (Graph_flags(5)==1): we want that 1 will represent SDRP, and 0 will symbolize DRP 
        % (Graph_flags(5)==2): we want that 1 will represent expenditure  
        G.Edges.Weight = G.Edges.Weight - 1;
        % here we want that weak relations that are in fact not a relation
        % but identical bundles, will be represented by dashed lines
        is_identical = zeros(1, length(G.Edges.Weight));
        if graph_type==1
            % do nothing; no identicals left
        else
            % find identical bundles and "mark" them as such
            for edge_index=1:length(G.Edges.Weight)
                if identical_choice_relaxed(G.Edges.EndNodes(edge_index,1) , G.Edges.EndNodes(edge_index,2))
                    is_identical(edge_index) = 1;
                end
            end
        end
        line_style = cell(1,length(G.Edges.Weight));
        line_style(:) = {'-'};
        line_style(is_identical~=0) = {'--'};
        
        % names of vertices/observations - {'(1)' '(2)' '(3)' ... }'
        names = cell(length(Graph_Matrix),1);
        
        if graph_type==1
            for i=1:length(Graph_Matrix)
                current_bundle = bundles_chosen_merged_identicals(i,:);
                current_bundle_str = '(';
                for j=1:length(current_bundle)
                    if j < length(current_bundle)
                        current_bundle_str = [current_bundle_str , num2str(current_bundle(j), '%10.4g'), ' , ']; %#ok<AGROW>
                    else
                        current_bundle_str = [current_bundle_str , num2str(current_bundle(j), '%10.4g'), ')']; %#ok<AGROW>
                    end
                end
                names{i} = current_bundle_str;
            end
        else   % graph_type==2 || graph_type==3 || graph_type==4           
            for i=1:length(Graph_Matrix)
                names{i} = strcat('(', num2str(i), ')');
            end
        end
        G.Nodes.Name = names;
        
        % font size of node label
        if obs_num <= 15
            node_font_size = 10;
            edge_font_size = 7.0;
        elseif obs_num <= 25
            node_font_size = 9;
            edge_font_size = 6.0;
        elseif obs_num <= 50
            node_font_size = 8;
            edge_font_size = 5.0;
        elseif obs_num <= 100
            node_font_size = 7;
            edge_font_size = 4.5;
        elseif obs_num <= 200
            node_font_size = 6;
            edge_font_size = 4.0;
        else   % obs_num > 200
            node_font_size = 5;
            edge_font_size = 3.5;
        end
        
        % print the revealed preference graph to .png file :
        f = figure('visible', 'off'); 
        if Graph_flags(5)==0 || graph_type==1
            RP_graph = plot(G, 'Layout',layout, 'LineWidth',(2*G.Edges.Weight+1)/3, 'LineStyle',line_style, 'NodeLabel',names, 'NodeFontSize',node_font_size, 'EdgeFontSize',edge_font_size);
        else   % elseif Graph_flags(5)==1 && graph_type~=1
            %RP_graph = plot(G, 'Layout',layout, 'EdgeLabel',G.Edges.Weight*(-1));
            weights = G.Edges.Weight*(-1);
            weights_labels = cell(length(weights),1);
            for w=1:length(weights)
                weights_labels{w} = num2str(weights(w), '%10.3f');
            end
            RP_graph = plot(G, 'Layout',layout, 'EdgeLabel',weights_labels, 'LineStyle',line_style, 'NodeLabel',names, 'NodeFontSize',node_font_size, 'EdgeFontSize',edge_font_size);
        end
        %RP.XData
        %RP.YData
        RP_graph_file_name = strcat(RP_graph_file_name_base, '-', RP_graph_file_name_extra{graph_type});
        saveas(RP_graph, strcat(RP_graph_file_name, '.png'));
        close(f);
    end
end

end