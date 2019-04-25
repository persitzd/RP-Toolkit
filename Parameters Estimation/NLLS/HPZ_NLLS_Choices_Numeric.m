function [choice_matrix, param] = HPZ_NLLS_Choices_Numeric (param, prices, endowments, treatment, function_flag, asymmetric_flag, pref_class, debugger_mode)

% This function finds and returns the optimal choice (a bundle) of the DM,
% given his/hers preferences, prices, and other things.
% In case of multiple optimal choices - it returns the optimal bundle that
% is closest to the observed bundle.

% The output is a matrix. The number of rows is the number of observations, 
% while the number of columns is the number of goods (currently 2).
% Each row is a vector of optimal choices (one quantity for each good).

% "param" is returned because it may change (be rounded) during the
% calculations

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



if pref_class == HPZ_Constants.risk_pref   % risk preferences
    
    % we need these threshold in order to avoid precision of print problems. 
    % for more detailed, see the documentation in "HPZ_Constants".
    if abs(param(1) - (-1)) < HPZ_Constants.print_threshold
        param(1) = -1;
    end
    if abs(param(2) - 1) < HPZ_Constants.print_threshold
        param(2) = 1;
    end
    
    % prob_x is the probability of the account x
    prob_x = 1/2;
    
    if asymmetric_flag == 1
        if treatment == 2
            prob_x = 1/3;
        elseif treatment == 3
            prob_x = 2/3;
        end
    end
    
    [observations,~] = size(prices);
    
    choice_matrix = zeros(observations, 2);
    
    if function_flag == HPZ_Constants.CRRA_func   % CRRA
        
        % this loop goes one set (of prices and endowment) after another,
        % and calculates the optimal choice for each of them
        parfor i=1:observations       
            
            % Threshold initialization
            % In the first stage we run two optimizations with two different starting points.
            % If the results are close enough, we interpret them as equal, and we use it 
            % as the prediction. 
            % The "close enough" criterion is based on the following thresholds   
            ratio_diff = 0.05;   % Threshold in terms of the ratio between the results
            Diff_threshod = 10^(-2);   % Threshold in terms of the difference between the results
            
            %% Choose two points A and B near the boundaries from the budget line:
            
            % compute the maximum quantity of x1 available for the given budget limit
            x_1_init = endowments(i) / prices(i,1);
            
            % The (2x2) matrix of two initial points: X_init
            X_init = zeros(2,2);
            
            % Point A - The 1st row of the X_init matrix:
            % Divide the budget line to 1000 pieces and take the point
            % that is 102 pieces away from the y axis, and 898 pieces away 
            % from the x axis:
            X_init(1,1) = (102/(1000+1)) * x_1_init; % (x1,_)
            X_init(1,2) = (endowments(i) - (X_init(1,1)*prices(i,1)))/prices(i,2); % (_,x2)
            
            % Point B - The 2nd row of the X_init matrix:
            % Divide the budget line to 1000 pieces and and take the point
            % that is 898 pieces away from the y axis, and 102 pieces away 
            % from the x axis:
            X_init(2,1) = (898/(1000+1)) * x_1_init;  % (x1,_)
            X_init(2,2) = (endowments(i) - (X_init(2,1)*prices(i,1)))/prices(i,2); % (_,x2)
            
            % setting the algoritm that will be used
            opts_phase1 = optimset('Algorithm','sqp','Display','off');        
            
            % Optimal solution matrix (2x3) for points A and B: [x1,x2,fval]
            Optimal_AB = zeros(2,3);
            
            % Run fmincon using:
            % 1. Default HPZ_Risk_Utility_Helper function as objective function
            % 2. Subject to EQUALITY constraint: Ax = b
            % 3. lower bound set to : [0 0]
            % 4. TolX = 0
            % 5. Initial points set to A = [X_init(1,:)] AND B=[X_init(2,:)]
            [Optimal_AB(1, 1:2), Optimal_AB(1, 3)] = fmincon (@(x) HPZ_Risk_Utility_Helper(prob_x, x, param(1), param(2), 0, function_flag), X_init(1,:), [], [], prices(i,:), endowments(i), [0 0], [], [], opts_phase1);
            [Optimal_AB(2, 1:2), Optimal_AB(2, 3)] = fmincon (@(x) HPZ_Risk_Utility_Helper(prob_x, x, param(1), param(2), 0, function_flag), X_init(2,:), [], [], prices(i,:), endowments(i), [0 0], [], [], opts_phase1);
            
            % calculating difference between the 2 solutions,
            % in terms of difference and in terms of difference in ratio
            Diff_Xopts_AB = abs(Optimal_AB(1,1:2) - Optimal_AB(2,1:2));
            ratio_Xopts_AB = abs((Optimal_AB(1,1:2) ./ Optimal_AB(2,1:2)) - 1);
            
            % if the distance between two optimal solutions is less than the
            % threshold (less than at least one of the thresholds):
            if (((Diff_Xopts_AB(1) <= Diff_threshod) && (Diff_Xopts_AB(2) <= Diff_threshod)) || ...
                (ratio_Xopts_AB(1) <= ratio_diff)    && (ratio_Xopts_AB(2) <= ratio_diff))
            else
                % compute the minimum between these two optimal values
                Fval_Diff_phase1 = abs(Optimal_AB(1,3) - Optimal_AB(2,3));
                
                if (Fval_Diff_phase1 > 0)
                    % Alpha_Factor:
                    % Fval Multiplier to remove the effect of penalizer factor
                    % inside the score function:
                    % utility multiplier (used inside GlobalSearch objective
                    % function)
                    if (Fval_Diff_phase1 >= 1)
                        Alpha_Factor = 10^(6);
                    else
                        % compute the number of digits after the "0."
                        digit_num = ceil(abs(log10(Fval_Diff_phase1)));
                        % Fval Multiplier to remove the effect of penalizer 
                        % factor inside the score function:
                        Alpha_Factor = 10^(digit_num + 5);
                    end
                    
                    [Optimal_AB(1, 1:2), Optimal_AB(1, 3)] = fmincon (@(x) HPZ_Risk_Utility_Helper(prob_x, x, param(1), param(2), 0, function_flag, Alpha_Factor), X_init(1,:), [], [], prices(i,:), endowments(i), [0 0], [], [], opts_phase1);
                    [Optimal_AB(2, 1:2), Optimal_AB(2, 3)] = fmincon (@(x) HPZ_Risk_Utility_Helper(prob_x, x, param(1), param(2), 0, function_flag, Alpha_Factor), X_init(2,:), [], [], prices(i,:), endowments(i), [0 0], [], [], opts_phase1);
                end
            end
            
            % find the best choice of the ones that were found
            [~,index_tt] = min(Optimal_AB(:,3));
            choice_matrix(i,:) = Optimal_AB(index_tt,1:2);
            
        end
        
    elseif function_flag == HPZ_Constants.CARA_func    % CARA
        
        % this loop goes one set (of prices and endowment) after another,
        % and calculates the optimal choice for each of them
        parfor i=1:observations       
            
            % Threshold initialization
            % In the first stage we run two optimizations with two different starting points.
            % If the results are close enough, we interpret them as equal, and we use it 
            % as the prediction. 
            % The "close enough" criterion is based on the following thresholds   
            ratio_diff = 0.05;   % Threshold in terms of the ratio between the results
            Diff_threshod = 10^(-2);   % Threshold in terms of the difference between the results
            
            %% Choose two points A and B near the boundaries from the
            %% budget line:
            
            % compute the maximum quantity of x1 available for the given budget limit
            x_1_init = endowments(i)/prices(i,1);
            
            % The (2x2) matrix of two initial points: X_init
            X_init = zeros(2,2);
            
            % Point A - The 1st row of the X_init matrix:
            % Divide the budget line to 1000 pieces and take the point
            % that is 102 pieces away from the y axis, and 898 pieces away 
            % from the x axis:
            X_init(1,1) = (102/(1000+1)) * x_1_init; % (x1,_)
            X_init(1,2) = (endowments(i)-(X_init(1,1)*prices(i,1)))/prices(i,2); % (_,x2)
            
            % Point B - The 2nd row of the X_init matrix:
            % Divide the budget line to 1000 pieces and and take the point
            % that is 898 pieces away from the y axis, and 102 pieces away 
            % from the x axis:
            X_init(2,1) = (898/(1000+1)) * x_1_init;  % (x1,_)
            X_init(2,2) = (endowments(i)-(X_init(2,1)*prices(i,1)))/prices(i,2); % (_,x2)
            
            % setting the algoritm that will be used
            opts_phase1 = optimset('Algorithm','sqp','Display','off');        
            
            % Optimal solution matrix (2x3) for points A and B: [x1,x2,fval]
            Optimal_AB = zeros(2,3);
            
            % Run fmincon using:
            % 1. Default HPZ_Risk_Utility_Helper function as objective function
            % 2. Subject to EQUALITY constraint: Ax = b
            % 3. lower bound set to : [0 0]
            % 4. TolX = 0
            % 5. Iinitial points set to A=[X_init(1,:)] AND B=[X_init(2,:)]
            [Optimal_AB(1, 1:2), Optimal_AB(1, 3)] = fmincon (@(x) HPZ_Risk_Utility_Helper(prob_x, x, param(1), 0, param(2), function_flag), X_init(1,:), [], [], prices(i,:), endowments(i), [0 0], [], [], opts_phase1);
            [Optimal_AB(2, 1:2), Optimal_AB(2, 3)] = fmincon (@(x) HPZ_Risk_Utility_Helper(prob_x, x, param(1), 0, param(2), function_flag), X_init(2,:), [], [], prices(i,:), endowments(i), [0 0], [], [], opts_phase1);
            
            % calculating difference between the 2 solutions,
            % in terms of difference and in terms of difference in ratio
            Diff_Xopts_AB = abs(Optimal_AB(1,1:2) - Optimal_AB(2,1:2));
            ratio_Xopts_AB = abs((Optimal_AB(1,1:2) ./ Optimal_AB(2,1:2)) - 1);
            
            % if the distance between two optimal solutions is less than the
            % threshold (less than at least one of the thresholds):
            if (((Diff_Xopts_AB(1) <= Diff_threshod) && (Diff_Xopts_AB(2) <= Diff_threshod)) || ...
                (ratio_Xopts_AB(1) <= ratio_diff)    && (ratio_Xopts_AB(2) <= ratio_diff))
            else
                % compute the minimum between these two optimal values
                Fval_Diff_phase1 = abs(Optimal_AB(1,3) - Optimal_AB(2,3));
                
                if (Fval_Diff_phase1 > 0)
                    % Alpha_Factor:
                    % Fval Multiplier to remove the effect of penalizer factor
                    % inside the score function:
                    % utility multiplier (used inside GlobalSearch objective
                    % function)
                    if (Fval_Diff_phase1 >= 1)
                        Alpha_Factor = 10^(6);
                    else
                        % compute the number of digits after the "0."
                        digit_num = ceil(abs(log10(Fval_Diff_phase1)));
                        % Fval Multiplier to remove the effect of penalizer 
                        % factor inside the score function:
                        Alpha_Factor = 10^(digit_num + 5);
                    end
                    
                    [Optimal_AB(1, 1:2), Optimal_AB(1, 3)] = fmincon (@(x) HPZ_Risk_Utility_Helper(prob_x, x, param(1), 0, param(2), function_flag, Alpha_Factor), X_init(1,:), [], [], prices(i,:), endowments(i), [0 0], [], [], opts_phase1);
                    [Optimal_AB(2, 1:2), Optimal_AB(2, 3)] = fmincon (@(x) HPZ_Risk_Utility_Helper(prob_x, x, param(1), 0, param(2), function_flag, Alpha_Factor), X_init(2,:), [], [], prices(i,:), endowments(i), [0 0], [], [], opts_phase1);
                    
                end
            end
            
            % find the best choice of the ones that were found
            [~,index_tt] = min(Optimal_AB(:,3));
            choice_matrix(i,:) = Optimal_AB(index_tt,1:2);
            
        end   % end of loop
        
    end
    
elseif pref_class == HPZ_Constants.OR_pref   % other regarding preferences
    
    % we need these threshold in order to avoid precision of print problems. 
    % for more detailed, see the documentation in "HPZ_Constants".
    if abs(param(1) - 1) < HPZ_Constants.print_threshold
        param(1) = -1;
    end
    if abs(param(2) - 1) < HPZ_Constants.print_threshold
        param(2) = 1;
    end
    
    [observations,~] = size(prices);
    
    choice_matrix = zeros(observations,2);
    
    if function_flag == HPZ_Constants.CES_func   % CES
        
        % this loop goes one set (of prices and endowment) after another,
        % and calculates the optimal choice for each of them
        parfor i=1:observations  
            
            % Threshold initialization
            % In the first stage we run two optimizations with two different starting points.
            % If the results are close enough, we interpret them as equal, and we use it 
            % as the prediction. 
            % The "close enough" criterion is based on the following thresholds   
            ratio_diff = 0.05;   % Threshold in terms of the ratio between the results
            Diff_threshod = 10^(-2);   % Threshold in terms of the difference between the results
            
            %% Choose two points A and B near the boundaries from the
            %% budget line:
            
            % compute the maximum quantity of x1 available for the given budget limit
            x_1_init = endowments(i)/prices(i,1);
            % The (2x2) matrix of two initial points: X_init
            X_init = zeros(2,2);
            
            % Point A - The 1st row of the X_init matrix:
            % Divide the budget line to 1001 pieces and take the 102nd point:       
            X_init(1,1) = (102/(1000+1)) * x_1_init; % (x1,_)
            X_init(1,2) = (endowments(i)-(X_init(1,1)*prices(i,1)))/prices(i,2); % (_,x2)
            
            % Point B - The 2nd row of the X_init matrix:
            % Divide the budget line to 1001 pieces and take the 898th point:
            X_init(2,1) = (898/(1000+1)) * x_1_init;  % (x1,_)
            X_init(2,2) = (endowments(i)-(X_init(2,1)*prices(i,1)))/prices(i,2); % (_,x2)
            opts_phase1 = optimset('Algorithm','sqp','Display','off');   
            
            % Optimal solution matrix (2x3) for point A and B: [x1,x2,fval]
            Optimal_AB = zeros(2,3);
            
            % Run fmincon using:
            % 1. Default HPZ_OR_Utility_Helper function as objective function
            % 2. Subject to EQUALITY constraint: Ax = b
            % 3. lower bound set to : [0 0]
            % 4. TolX = 0
            % 5. Iinitial points set to A=[X_init(1,:)] AND B=[X_init(2,:)]
            [Optimal_AB(1, 1:2), Optimal_AB(1, 3)] = fmincon (@(x) HPZ_OR_Utility_Helper(x, param(1), param(2), function_flag), X_init(1,:), [], [], prices(i,:), endowments(i), [0 0], [], [], opts_phase1);
            [Optimal_AB(2, 1:2), Optimal_AB(2, 3)] = fmincon (@(x) HPZ_OR_Utility_Helper(x, param(1), param(2), function_flag), X_init(2,:), [], [], prices(i,:), endowments(i), [0 0], [], [], opts_phase1);
            
            % calculating difference between the 2 solutions,
            % in terms of difference and in terms of difference in ratio
            Diff_Xopts_AB = abs(Optimal_AB(1,1:2) - Optimal_AB(2,1:2));
            ratio_Xopts_AB = abs((Optimal_AB(1,1:2) ./ Optimal_AB(2,1:2)) - 1);
            
            % if the distance between two optimal solutions is less than the
            % threshold (less than at least one of the thresholds):
            if (((Diff_Xopts_AB(1) <= Diff_threshod) && (Diff_Xopts_AB(2) <= Diff_threshod)) || ...
                (ratio_Xopts_AB(1) <= ratio_diff)    && (ratio_Xopts_AB(1) <= ratio_diff))
            else
                % compute the minimum between these two optimal values
                Fval_Diff_phase1 = abs(Optimal_AB(1,3) - Optimal_AB(2,3));
                if (Fval_Diff_phase1 > 0)
                    % Alpha_Factor:
                    % Fval Multiplier to remove the effect of penalizer factor
                    % inside the score function:
                    % utility multiplier (used inside GlobalSearch objective
                    % function)
                    if (Fval_Diff_phase1 >= 1)
                        Alpha_Factor = 10^(6);
                    else
                        % compute the number of digits after the "0."
                        digit_num = ceil(abs(log10(Fval_Diff_phase1)));
                        % Fval Multiplier to remove the effect of penalizer 
                        % factor inside the score function:
                        Alpha_Factor = 10^(digit_num + 5);
                    end
                    [Optimal_AB(1, 1:2), Optimal_AB(1, 3)] = fmincon (@(x) HPZ_OR_Utility_Helper(x, param(1), param(2), function_flag, Alpha_Factor), X_init(1,:), [], [], prices(i,:), endowments(i), [0 0], [], [], opts_phase1);
                    [Optimal_AB(2, 1:2), Optimal_AB(2, 3)] = fmincon (@(x) HPZ_OR_Utility_Helper(x, param(1), param(2), function_flag, Alpha_Factor), X_init(2,:), [], [], prices(i,:), endowments(i), [0 0], [], [], opts_phase1); %#ok<*PFBNS>
                
                end
            end
            
            % find the best choice of the ones that were found
            [~,index_tt] = min(Optimal_AB(:,3));
            choice_matrix(i,:) = Optimal_AB(index_tt, 1:2);
            
        end   % end of loop
        
    end
    
end



% warnings to consol for debugging purposes (only when debugger mode is activated) 
% (it is in a seperate loop, in order not to disturb the "parfor" parallel computing)  
for i=1:observations
    if debugger_mode
        if (isnan(choice_matrix(i,1)) || isnan(choice_matrix(i,2)))
            warning('At least one of the commodities in the optimal bundle is NaN when Param(2) is %.20g, Param(1) is %.20g, p1 is %.20g, p2 is %.20g, x1 is %.20g, x2 is %.20g and i is %d.', param(2), param(1), prices(i,1), prices(i,2), observations(i,1), observations(i,2), i);
        elseif (isinf(choice_matrix(i,1)) || isinf(choice_matrix(i,2)))
            warning('At least one of the commodities in the optimal bundle is +Inf or -Inf when Param(2) is %.20g, Param(1) is %.20g, p1 is %.20g, p2 is %.20g, x1 is %.20g, x2 is %.20g and i is %d.', param(2), param(1), prices(i,1), prices(i,2), observations(i,1), observations(i,2), i);
        end
    end
end



end