function [subsets_of_observations, num_of_subsets] = HPZ_Indices_Problem_Distribute(relevant_RP, relevant_SRP, varargin)

% NOTE: This function is used to separate the observations (nodes) to 
% "strongly connected components". But in our case it might be a bit
% comlicated, as we don't have just one graph but two graphs
% (one for DRP and one for SDRP).

% When calculating various Inconsistency Indices (e.g. Varian, Houtman-Maks), 
% the calculations can be performed separately for disjoint subsets of
% observations, then in the end aggregated together.
%
% Such distribution of the calculation problem to smaller problems, assists
% and improves running time in 2 aspects:
%   (1) if the complexity of number of iterations of the algorithm is
%       O(f(n)), and normally in complicated indices f is either polynomial
%       with a high degree or even exponential, reducing the size of n by
%       this distribution significantly reduces running time.
%   (2) the inner complexity of each iteration is smaller, since the sizes
%       of all matrices are smaller.
%   * Also any subset that has only 1 observation is a subset that does not
%     affect the index, therefore it requires no further calculations.
%
% In order to achieve all these great fits, this function finds the best
% possible distribution (distribution to smallest subsets), and returns
% the subsets.
%
% The distribution to subsets differs depending on the index and the
% calculation problem related to it, for example: 
% * in the case of the Varian Index, we would like to have disjoint subsets 
%   such that there is no cycle of strict-revealed-preferred relations that 
%   involves observations from different subsets.
% * in the case of the Houtman-Maks Index with 3 or more goods, we would 
%   like to have disjoint subsets such that there is no GARP-violating 
%   cycle that involves observations from different subsets.
% * in the case of the Houtman-Maks Index with 2 goods, we would like to 
%   have disjoint subsets such that there is no GARP-violating cycle of
%   length 2 that involves observations from different subsets.
% The parameters "relevant_RP" & "relevant_SRP" allows the distribution
% according to different types of cycles and relations.



if ~isempty(varargin)
    % the varargin allows the user to directly pass this function
    % the relevant_GARP_couples matrix, instead of the relevant_RP & relevant_SRP 
    relevant_GARP_couples = varargin{1};
    
    % number of observations
    [obs_num , obs_num2] = size(relevant_GARP_couples);
    
    if obs_num ~= obs_num2
        error('number of rows and number of column of input matrix should be equal');
    end
else
    % calculate the relevant GARP relation
    % (relevant GARP is without transitivity in Houtman-Maks when there are 
    % only 2 goods (thanks to the Rose (1958) Theorem), and it is with 
    % transitivity in Houtman-Maks with more than 2 goods, also it is with 
    % transitivity when we calculate Varian)
    relevant_GARP = relevant_RP .* relevant_SRP';

    % number of observations
    [obs_num , ~] = size(relevant_GARP);
    
    % making sure an observation is not related to itself (just in case)
    for i=1:obs_num
        relevant_GARP(i,i) = 0;
    end

    % we want a matrix in which (i,j) is 1 iff there is a cycle
    % which contains observations i and j (regardless of which of them is
    % reveald preferred and which is strictly revealed preferred)
    relevant_GARP_couples = max(relevant_GARP, relevant_GARP');
    %copy_GARP_couples = relevant_GARP_couples;
end





% we want to divide the observations to subsets, so that none of the
% observations in subset X is in a cycle with any observation
% in subset Y, for every two different subsets X and Y
subsets_of_observations = cell(0);

% observations that were yet assigned to a subset
obs_left = 1:obs_num;

% the current subset we assign observations to
current_subset = 1;

while sum(obs_left) > 0
    
    % the current observation we assign to a subset
    current_obs = find(obs_left, 1, 'first');
    
    % initialization of this subset
    temp_subset = current_obs;
    
    % all observations that have a cycle with this observation
    related_obs = find(relevant_GARP_couples(current_obs,:));
    % make this row zero, we won't be needing it again
    relevant_GARP_couples(current_obs,:) = 0;
    
    % now we want to add the related observations to the subset, but also
    % the observations that are related to them, and so on
    while ~isempty(related_obs)
        
        % add this related observation to the subset
        temp_subset = [temp_subset, related_obs(1)]; %#ok<AGROW>
        
        % add its related observations to the related_obs vector
        related_obs_temp = find(relevant_GARP_couples(related_obs(1),:));
        % make this row zero, we won't be needing it again
        relevant_GARP_couples(related_obs(1),:) = 0;
        
        % deleting observations that are already in temp_subset
        j = 1;
        while j <= length(related_obs_temp)
            if any(related_obs_temp(j) == temp_subset)   % ismember(related_obs_temp(j),temp_subset)
                related_obs_temp = related_obs_temp([1:(j-1) , (j+1):end]);
            else
                j = j + 1;
            end
        end
        
        % add the new related observations
        related_obs = [related_obs, related_obs_temp]; %#ok<AGROW>
        % delete the current related obs that was handled
        related_obs = related_obs(2:end);
        
        % deleting duplicates (use 'stable' to keep original order, but there's no need for this here) 
        related_obs = unique([related_obs, related_obs_temp]);
    end
    
    % update the cell of subsets and its counter
    subsets_of_observations{current_subset} = temp_subset;
    current_subset = current_subset + 1;
    
    % deleting the observations that were assigned to this subset
    obs_left(temp_subset) = 0;
    
%     if length(temp_subset) > 1
%         temp_subset
%         temp_GARP_couples = copy_GARP_couples(temp_subset , temp_subset)
%     end
end



% final number of subsets
num_of_subsets = current_subset - 1;


end