function [choice_matrix, param] = HPZ_NLLS_Choices_Analytic(param, observations, function_flag, pref_class)

% This function finds and returns the optimal choice (a bundle) of the DM,
% given his/hers preferences, prices, and other things.
% In case of multiple optimal choices - it returns the optimal bundle that
% is closest (in euclidean metric) to the observed bundle.

% The output is a matrix. The number of rows is the number of observations, 
% while the number of columns is the number of goods (currently 2).
% Each row is a vector of optimal choices (one quantity for each good).

% "param" is returned because it may change (be rounded) during the
% calculations

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% number of observations
[num_obs,~] = size(observations);

% a matrix of intersection points of budget lines,
% which also represents the max quantity that can be
% purchased under budget restrictions
max_vals = 1 ./ observations(:,3:4);

% extract data of observed bundles (quantity of each good - (x1,x2)) 
% from the observation matrix
observed_bundle = observations(:,1:2);

% set prices from observation matrix
%prices = observations(:,3:4);

% NOTE: these code lines are currently redundant
% the percentage (number between 0 to 1) of the endowment that
% the DM used in each of his choices. for example, if the DM used 60$ to
% buy the quantities of x1 and x2, but he had 80$ as endowment, then this
% variable will be 60/80 = 0.75
%endowment_used = observations(:,1) .* observations(:,3) + observations(:,2) .* observations(:,4);

% the output matrix initialization
choice_matrix = zeros(num_obs,2);

%% risk preferences
if pref_class == HPZ_Constants.risk_pref
    
    % we need these threshold in order to avoid precision of print problems. 
    % for more detailed, see the documentation in "HPZ_Constants".
    if abs(param(1) - (-1)) < HPZ_Constants.print_threshold
        param(1) = -1;
    end
    if abs(param(2) - 1) < HPZ_Constants.print_threshold
        param(2) = 1;
    end
    
    %% CRRA
    if function_flag == HPZ_Constants.CRRA_func
        
        % it is more convenient to use these notations
        beta = param(1);
        rho  = param(2);
        
        % update "param" to be returned (currently not needed in CRRA)
        %param = [beta , rho];
        
        % this loop goes one set (of prices and endowment) after another,
        % and calculates the optimal choice for each of them
        % cannot use parfor any more since x2* computation depends on x1*
        for i=1:num_obs 
            

%             % THIS PART IS CURRENTLY NOT IN USE
%             % we now calculate the "adjusted beta", that is, the theoretical
%             % beta that would lead, with prob_x = alpha = 0.5, to the same
%             % results for the utility function
%             
%             prob_x = 0.8;
%             
%             % By Gul (1991), in the case of two prizes, alpha is the probability of the
%             % higher prize.
%             if max_vals(i, 1) > max_vals(i, 2)
%                 alpha = prob_x;
%             else
%                 alpha = 1 - prob_x;
%             end
%             
%             % By Gul (1991), the weights are calculated using gamma(alpha)
%             gamma = alpha / ( 1 + ((1-alpha)*beta) );
%             
%             param(1) = (1 / gamma) - 2;
            
            
            % p = p1 / p2 = Max_Y / Max_X
            p = max_vals(i, 2) / max_vals(i,1);
            
            
            %% if Rho = 0 and Beta >= 0 CASE #4
            if beta >= 0 && rho == 0

                % In this section we did something slightly different
                % from the original analytical solution. in the original 
                % analytical solution, there's a special treatment when p
                % equals to (1+beta) or 1/(1+beta), which makes every bundle
                % on the budget line from the cheap bundle's corner to the
                % bundle with equal amounts, to be an optimal bundle.
                % This results that beta that equals *precisely* to some
                % specific price ratio minus 1, gains unfair advantage (what
                % are the chances that we "scored" to a price ratio which is
                % *exactly* (1+beta)?), therefore we did not account for
                % this case, but rather we joined p = 1+beta with p > 1+beta 
                % to make a one case of p >= 1+beta, and equivalently we
                % made a case of p <= 1/(1+beta) instead of 2 cases.
                % An Exception: when beta=0 (CASE 5), we decided to accept 
                % the idea that all the budget line is an optimal choice.
                
                % if p = 1, beta = 0, rho = 0, then every choice on the
                % budget line is optimal (CASE 5)
                if (p == 1) && ( beta == 0 || (beta > 0 && beta <= 10^(-6)) )
                    
                    % every choice on the budget line is an optimal choice
                    choice_matrix(i,:) = observed_bundle(i,:); % Error is 0
                    % if (1 + Beta) < p   OR   (1 + Beta) = p
                    
                elseif p >= (1 + beta)
                    
                    choice_matrix(i,:) = [0, max_vals(i,2)]; % (0, max_y)
                    % if p < 1 / (1 + Beta)   OR   p = 1 / (1 + Beta)
                    
                elseif p <= 1 / (1 + beta)
                    
                    choice_matrix(i,:) = [max_vals(i,1), 0]; % (max_x, 0)
                    % if (1 / (1 + Beta)) < p < (1 + Beta)
                    
                elseif p > 1 / (1 + beta) && p < (1 + beta)
                    
                    choice_matrix(i,:) = [max_vals(i,2) / (p + 1), max_vals(i,2) / (p + 1)]; % (max_y/(p+1), max_y/(p+1))
                
                end
                
            %% if  (rho > 0 and beta = -1) OR (Rho = 0 and -1 <= Beta < 0)    CASE #3
            elseif (beta < 0 && rho == 0) || (rho > 0 && beta == -1)
                
                if p < 1
                    
                    choice_matrix(i,:) = [max_vals(i,1), 0]; % (max_x, 0)
                    
                elseif p == 1
                    
                    % each of the 2 corner choices is an
                    % optimal choice. Therefore the error is the distance
                    % between the observed bundle to the *nearest* corner.
                    if observed_bundle(i, 1) >= observed_bundle(i,2)
                        
                        %(max_x, 0) is the closer corner
                        choice_matrix(i,:) = [max_vals(i,1), 0];
                        
                    else
                        
                        % (0, max_y) is the closer corner
                        choice_matrix(i,:) = [0, max_vals(i,2)];
                        
                    end
                    
                elseif p > 1
                    
                    choice_matrix(i,:) = [0, max_vals(i,2)]; % (0, max_y)
                    
                end

            %% if Rho > 0 and Beta >= 0 CASE #1
            elseif beta >= 0 && rho > 0

                % if p < 1 / (1 + Beta)
                if p < 1 / (1 + beta)
                    
                    % x1* = Max_x / [1 + ( (p*(1+beta))^(1/rho) ) / p]
                    choice_matrix(i,1) = max_vals(i,1) / ( 1 + ( (p*(1+beta))^(1/rho) ) / p );
                    % x2* = x1* *((1+Beta)*p)^(1/Rho)
                    %choice_matrix(i,2) = choice_matrix(i,1)*((1+beta)*p)^(1/rho);
                    % x2* = Max_y / [1 + p / ( (p(1+beta))^(1/rho) )]
                    choice_matrix(i,2) = max_vals(i,2) / (1 + p / ( (p*(1+beta))^(1/rho) ));
                % if 1 / (1 + Beta) <= p <= 1 + Beta
                
                elseif p >= 1 / (1 + beta) && p <= 1 + beta
                    
                    % x1* = (max_y)/(p+1)
                    choice_matrix(i,1) = (max_vals(i,2) ) / ( p + 1 ) ;
                    % x2* = (max_y)/(p+1)
                    choice_matrix(i,2) = (max_vals(i,2) ) / ( p + 1 ) ;
                % if p > 1 / (1 + Beta)
                
                elseif p > 1 + beta
                    
                    % x1* = (max_x)/(1 + ((p/(1+Beta))^(1/Rho))/p )
                    choice_matrix(i,1) = (max_vals(i,1) ) / ( 1 + ( ( p/(1+beta) )^(1/rho) )/ p ) ;
                    % x2* = x1* *(p/(1+Beta))^(1/Rho)
                    %choice_matrix(i,2) = choice_matrix(i,1) * ( ( p/(1+param(1)) )^(1/param(2)) );
                    % x2* = Max_Y / (1 + p / ( (p/(1+beta))^(1/rho) ) )
                    choice_matrix(i,2) = max_vals(i,2) / (1 + p / ( (p/(1+beta))^(1/rho) ) );
                    
                end
                
            %% if Rho > 0 and -1 < Beta < 0 CASE #2
            elseif rho > 0 && ( beta > -1 && beta < 0 )
                
                % if p < 1
                if p < 1
                    
                    % x1* = Max_x / [1 + (p*(1+Beta)^(1/rho))/p]
                    choice_matrix(i,1) = max_vals(i,1) / ( 1 + ( (p*(1+beta)) ^(1/rho) ) / p );
                    % x2* = x1* *((1+Beta)*p)^(1/Rho)
                    %choice_matrix(i,2) = choice_matrix(i,1)*((1+beta)*p)^(1/rho);
                    % x2* = Max_y / [1 + p / ( p(1+beta)^(1/rho) )]
                    choice_matrix(i,2) = max_vals(i,2) / (1 + p / ( (p*(1+beta))^(1/rho) ));
                % if p = 1
                
                elseif p == 1
                    
                    % there are 2 optimal choices, one with
                    % x1* >= x2* and one with x2* >= x1*.
                    % Therefore the error is the distance between the
                    % observed bundle to the *nearest* of these 2.
                    % Therefore, if the observed bundle has x1 >= x2, he
                    % must be closer to the optimal bundle with x1* >= x2*,
                    % and vice versa.
                    if observed_bundle(i, 1) >= observed_bundle(i, 2)
                        
                        % x1 >= x2, therefore we take the
                        % optimal bundle such that x1* >= x2*
                       
                        % x1* = Max / (1 + (1+beta)^(1/rho))
                        choice_matrix(i,1) = max_vals(i,1) / (1 + ( (1 + beta)^(1/rho) ) );
                        % x2* = Max / (1 + (1+beta)^(-1/rho))
                        choice_matrix(i,2) = max_vals(i,1) / (1 + ( (1 + beta)^(-1/rho) ) );
                        
                    else
                        
                        % x2 >= x1, therefore we take the
                        % optimal bundle such that x2* >= x1*
                        
                        % x1* = Max / (1 + (1+beta)^(-1/rho))
                        choice_matrix(i,1) = max_vals(i,1) / (1 + ( (1 + beta)^(-1/rho) ) );
                        % x2* = Max / (1 + (1+beta)^(1/rho))
                        choice_matrix(i,2) = max_vals(i,1) / (1 + ( (1 + beta)^(1/rho) ) );
                        
                    end
                    
                % if p > 1
                elseif p > 1
                    
                    % x1* = Max_X / (1 + ( (p/(1+beta))^(1/rho) ) / p )
                    choice_matrix(i,1) =  max_vals(i,1) / ( 1 + ( ( p/(1+beta) )^(1/rho) ) / p);
                    % x2* = x1* *(p/(1+beta))^(1/rho)
                    %choice_matrix(i,2) =  choice_matrix(i,1) * ( ( p/(1+beta) )^(1/rho) );
                    % x2* = Max_Y / (1 + p / ( (p/(1+beta))^(1/rho) ) )
                    choice_matrix(i,2) = max_vals(i,2) / (1 + p / ( (p/(1+beta))^(1/rho) ) );
                    
                end
                
            end
            
            % warnings to consol for debugging purposes (only when debugger mode is activated) 
            if (HPZ_Constants.debugger_mode)
                if (isnan(choice_matrix(i,1)) || isnan(choice_matrix(i,2)))
                    warning('At least one of the commodities in the optimal bundle is NaN when Rho is %.20g, Beta is %.20g, P is %.20g and i is %d.', rho, beta, p, i);
                elseif (isinf(choice_matrix(i,1)) || isinf(choice_matrix(i,2)))
                    warning('At least one of the commodities in the optimal bundle is +Inf or -Inf when Rho is %.20g, Beta is %.20g, P is %.20g and i is %d.', rho, beta, p, i);
                end
            end
            
        end   % end of loop over all observations

    %% CARA
    elseif function_flag == HPZ_Constants.CARA_func
        
        % it is more convenient to use these notations
        beta = param(1);
        A  = param(2);
        
        % update "param" to be returned (currently not needed in CARA)
        %param = [beta , A];
        
        % this loop goes one set (of prices and endowment) after another,
        % and calculates the optimal choice for each of them
        for i=1:num_obs
            
            max_1 = max_vals(i,1);
            max_2 = max_vals(i,2);
            
            % p = p1 / p2 = Max_Y / Max_X
            p = max_2 / max_1;
            
            % threshold for the CARA parameter (A). 
            % Cases 17, 18 and 19 in the DA-2 Document refer to the case where A approaches 0. 
            % However, the distinction between A values that are close to zero
            % and those that are not, depend on the value of x.
            % Therefore, the threshold in the code is on the multiplication of
            % A and x.
            % As a value for x we take min(x) and not max(x) because they are both in a
            % negative power and therefore min(x) is more dominant than max(x).
            % If A*min(x)>0 is less than this threshold, we use the cases as stated 
            % in the document. 
            % However if A*min(x)=0 we do the same using A*max(x).
            A_threshold = eps * 2^22;   % 2^(-30)
            
            %% A -> 0 - CASEs 17, 18, 19
            if A >= 0 && ( ( min(max_1, max_2) > 0 && A*min(max_1, max_2) < A_threshold ) || ( min(max_1, max_2) == 0 && A*max(max_1, max_2) < A_threshold ) )
                
                sum_prices = 1 / max_1 + 1 / max_2;
                
                if beta > 0   % CASE 17
                    
                    if p < (1 / (1+beta))
                        choice_matrix(i,:) = [max_1, 0];
                    elseif p >= (1 / (1+beta)) && p <= (1+beta) 
                        choice_matrix(i,:) = [1 / sum_prices, 1 / sum_prices];
                    else  % p > (1+beta) 
                        choice_matrix(i,:) = [0, max_2];
                    end
                    
                elseif beta == 0   % CASE 18
                    
                    if p < 1
                        choice_matrix(i,:) = [max_1, 0];
                    elseif p == 1 
                        % every choice on the budget line is an optimal choice
                        choice_matrix(i,:) = observed_bundle(i,:); % Error is 0
                    else  % p > 1
                        choice_matrix(i,:) = [0, max_2];
                    end
                    
                else % beta < 0 , CASE 19
                    
                    if p < 1
                        choice_matrix(i,:) = [max_1, 0];
                    elseif p == 1
                        % there are 2 optimal choices, one with
                        % x1* >= x2* and one with x2* >= x1*.
                        % Therefore the error is the distance between the
                        % observed bundle to the *nearest* of these 2.
                        % Therefore, if the observed bundle has x1 >= x2, he
                        % must be closer to the optimal bundle with x1* >= x2*,
                        % and vice versa.
                        if observed_bundle(i, 1) >= observed_bundle(i, 2)
                            % x1 >= x2, therefore we take the
                            % optimal bundle such that x1* >= x2*
                            choice_matrix(i,:) = [max_1, 0];
                        else
                            % x2 >= x1, therefore we take the
                            % optimal bundle such that x2* >= x1*
                            choice_matrix(i,:) = [0, max_2];
                        end
                    else  % p > 1
                        choice_matrix(i,:) = [0, max_2];
                    end
                    
                end
                
            %% non-negative beta - CASE 14 (CASE 17 and CASE 18 are included)
            elseif beta >= 0
                
                if p < (exp(-A*max_1) / (1+beta))
                    % x1* = Max_X
                    % x2* = 0
                    choice_matrix(i,:) = [max_1, 0];

                elseif p >= (exp(-A*max_1) / (1+beta)) && p < (1 / (1+beta))
                    % x1* = 1/(p+1) * (Max_Y - (1/a)*ln(p*(1+beta)))
                    choice_matrix(i,1) = (max_2 - (log(p*(1+beta))) / A ) / (1+p); 
                    % x2* = 1/(p+1) * (Max_Y + (p/a)*ln(p*(1+beta)))
                    choice_matrix(i,2) = (max_2 + (p*log(p*(1+beta))) / A ) / (1+p); 

                elseif p >= (1 / (1+beta)) && p <= (1+beta)
                    % x1* = Max_Y / (p+1)
                    % x2* = Max_Y / (p+1)
                    choice_matrix(i,:) = [max_2/(p+1), max_2/(p+1)];

                elseif p > (1+beta) && p <= (1+beta)*exp(A*max_vals(i,2))
                    % x1* = 1/(p+1) * (Max_Y - (1/a)*ln(p/(1+beta)))
                    choice_matrix(i,1) = (max_2 - (log(p/(1+beta))) / A ) / (1+p); 
                    % x2* = 1/(p+1) * (Max_Y + (p/a)*ln(p/(1+beta)))
                    choice_matrix(i,2) = (max_2 + (p*log(p/(1+beta))) / A ) / (1+p); 

                elseif p > (exp(A*max_2)*(1+beta))
                    % x1* = 0
                    % x2* = Max_Y
                    choice_matrix(i,:) = [0, max_2];
                end

            %% negative beta - CASE 15 (some of CASE 19 is included)
            elseif beta > -1 && beta < 0
                
                p1 = 1 / max_1;
                p2 = 1 / max_2;
                % how much can we buy if we take equal quantities from both goods 
                max_both = 1 / (p1 + p2);
                
                % this is the solution, assuming x1 >= x2, and without
                % the restrictions of x1,x2 >= 0
                x1_bigger_solution = [0,0];
                % x1* = 1/(p+1) * (Max_Y - (1/a)*ln(p*(1+beta)))
                x1_bigger_solution(1) = (max_2 - (log(p*(1+beta))) / A ) / (1+p);
                % x2* = 1/(p+1) * (Max_Y + (p/a)*ln(p*(1+beta)))
                x1_bigger_solution(2) = (max_2 + (p*log(p*(1+beta))) / A ) / (1+p);
                % if the solution is such that x2 < 0, then we "compromise"
                % on a solution with x2 = 0
                if x1_bigger_solution(2) < 0
                    x1_bigger_solution(1) = max_1;
                    x1_bigger_solution(2) = 0;
                % if the solution is such that x2 > x1, then we "compromise"
                % on a solution with x2 = x1
                elseif x1_bigger_solution(2) > x1_bigger_solution(1)
                    x1_bigger_solution(1) = max_both;
                    x1_bigger_solution(2) = max_both;
                end
                
                % this is the solution, assuming x2 >= x1, and without
                % the restrictions of x1,x2 >= 0
                x2_bigger_solution = [0,0];
                % x1* = 1/(p+1) * (Max_Y - (1/a)*ln(p/(1+beta)))
                x2_bigger_solution(1) = (max_2 - (log(p/(1+beta))) / A ) / (1+p); 
                % x2* = 1/(p+1) * (Max_Y + (p/a)*ln(p/(1+beta)))
                x2_bigger_solution(2) = (max_2 + (p*log(p/(1+beta))) / A ) / (1+p); 
                % if the solution is such that x1 < 0, then we "compromise"
                % on a solution with x1 = 0
                if x2_bigger_solution(1) < 0
                    x2_bigger_solution(1) = 0;
                    x2_bigger_solution(2) = max_2;
                % if the solution is such that x1 > x2, then we "compromise"
                % on a solution with x1 = x2
                elseif x2_bigger_solution(1) > x2_bigger_solution(2)
                    x2_bigger_solution(1) = max_both;
                    x2_bigger_solution(2) = max_both;
                end
                
                
                if p == 1
                    
                    % there are 2 optimal choices, one with 
                   % x1* >= x2* and one with x2* >= x1*.  
                   % Therefore the error is the distance between the 
                   % observed bundle to the *nearest* of these 2.
                   % Therefore, if the observed bundle has x1 >= x2, he
                   % must be closer to the optimal bundle with x1* >= x2*,
                   % and vice versa.
                   if observed_bundle(i, 1) >= observed_bundle(i, 2)
                       % x1 >= x2, therefore we take the
                       % optimal bundle such that x1* >= x2*
                       choice_matrix(i,1:2) = x1_bigger_solution(1:2);
                   else
                       % x2 >= x1, therefore we take the
                       % optimal bundle such that x2* >= x1*
                       choice_matrix(i,1:2) = x2_bigger_solution(1:2);
                   end
                    
                elseif p < 1
                    
                    % we take the optimal bundle such that x1* >= x2*,
                    % since x1 is cheaper than x2
                    choice_matrix(i,1:2) = x1_bigger_solution(1:2);
                    
                else % p > 1
                    
                    % we take the optimal bundle such that x2* >= x1*,
                    % since x2 is cheaper than x1
                    choice_matrix(i,1:2) = x2_bigger_solution(1:2);
                    
                end
                
                
%                 if p < (exp(-A*max_1) / (1+beta))
%                     % x1* = Max_X
%                     % x2* = 0
%                     choice_matrix(i,:) = [max_1, 0];
% 
%                 elseif p >= (exp(-A*max_1) / (1+beta)) && p < 1
%                     % x1* = 1/(p+1) * (Max_Y - (1/a)*ln(p*(1+beta)))
%                     choice_matrix(i,1) = (max_2 - (log(p*(1+beta))) / A ) / (1+p);
%                     % x2* = 1/(p+1) * (Max_Y + (p/a)*ln(p*(1+beta)))
%                     choice_matrix(i,2) = (max_2 + (p*log(p*(1+beta))) / A ) / (1+p); 
%                     
%                 elseif p == 1
%                    % there are 2 optimal choices, one with 
%                    % x1* >= x2* and one with x2* >= x1*.  
%                    % Therefore the error is the distance between the 
%                    % observed bundle to the *nearest* of these 2.
%                    % Therefore, if the observed bundle has x1 >= x2, he
%                    % must be closer to the optimal bundle with x1* >= x2*,
%                    % and vice versa.
%                    if observed_bundle(i, 1) >= observed_bundle(i, 2)
%                        % x1 >= x2, therefore we take the
%                        % optimal bundle such that x1* >= x2*
%                        
%                        % x1* = 1/2 * (Max_Y - (1/a)*ln(1+beta))
%                        choice_matrix(i,1) = (max_2 - log(1+beta) / A ) / 2;
%                        % x2* = 1/2 * (Max_Y + (1/a)*ln(1+beta))
%                        choice_matrix(i,2) = (max_2 + log(1+beta) / A ) / 2; 
%                    else
%                        % x2 >= x1, therefore we take the
%                        % optimal bundle such that x2* >= x1*
%                        
%                        % x1* = 1/2 * (Max_Y + (1/a)*ln(1+beta))
%                        choice_matrix(i,1) = (max_2 + log(1+beta) / A ) / 2;
%                        % x2* = 1/2 * (Max_Y - (1/a)*ln(1+beta))
%                        choice_matrix(i,2) = (max_2 - log(1+beta) / A ) / 2; 
%                    end
% 
%                 elseif p > 1 && p <= (1+beta)*exp(A*max_2)
%                     % x1* = 1/(p+1) * (Max_Y - (1/a)*ln(p/(1+beta)))
%                     choice_matrix(i,1) = (max_2 - (log(p/(1+beta))) / A ) / (1+p); 
%                     % x2* = 1/(p+1) * (Max_Y + (p/a)*ln(p/(1+beta)))
%                     choice_matrix(i,2) = (max_2 + (p*log(p/(1+beta))) / A ) / (1+p); 
% 
%                 elseif p > (exp(A*max_2) * (1+beta))
%                     % x1* = 0
%                     % x2* = Max_Y
%                     choice_matrix(i,:) = [0, max_2];
%                 end

            %% beta is equal to -1 - CASE 16 (some of CASE 19 is included)
            elseif beta == -1
                if p < 1
                    % x1* = Max_X
                    % x2* = 0
                    choice_matrix(i,:) = [max_1, 0];
                elseif p == 1
                    if observed_bundle(i, 1) >= observed_bundle(i, 2)
                        % x1* = Max_X
                        % x2* = 0
                        choice_matrix(i,:) = [max_1, 0];
                    else
                        % x1* = 0
                        % x2* = Max_Y
                        choice_matrix(i,:) = [0, max_2];
                    end
                elseif p > 1
                    % x1* = 0
                    % x2* = Max_Y
                    choice_matrix(i,:) = [0, max_2];
                end

            end   % end of conditions on beta

            % warnings to consol for debugging purposes (only when debugger mode is activated) 
            if (HPZ_Constants.debugger_mode)
                if (isnan(choice_matrix(i,1)) || isnan(choice_matrix(i,2)))
                    warning('At least one of the commodities in the optimal bundle is NaN when A is %.20g, Beta is %.20g, P is %.20g and i is %d.', A, beta, p, i);
                elseif (isinf(choice_matrix(i,1)) || isinf(choice_matrix(i,2)))
                    warning('At least one of the commodities in the optimal bundle is +Inf or -Inf when A is %.20g, Beta is %.20g, P is %.20g and i is %d.', A, beta, p, i);
                end
            end

        end   % end of loop over all observations

    end   % end of CARA

elseif pref_class == HPZ_Constants.OR_pref   % other regarding preferences
    
    % we need these threshold in order to avoid precision of print problems. 
    % for more detailed, see the documentation in "HPZ_Constants".
    if abs(param(1) - 1) < HPZ_Constants.print_threshold
        param(1) = 1;
    end
    if abs(param(2) - 1) < HPZ_Constants.print_threshold
        param(2) = 1;
    elseif abs(param(2) - (-1)) < HPZ_Constants.print_threshold
        param(2) = -1;
    end
    
    %% CES
    if function_flag == HPZ_Constants.CES_func

        % it is more convenient to use these notations
        alpha = param(1);
        rho = param(2);
        
        % this epsilon is the threshold for deciding when a value of alpha
        % or rho is close to enough to 0 or to 1, to be considered 0 or 1,
        % respectively
        epsilon = eps * 2^12;   % 2^(-40)
        
        % rounding rho to 0 or to 1 when needed
        if (abs(rho) < epsilon)
            rho = 0;
        elseif (abs(rho-1) < epsilon)
            rho = 1;
        end
        
        % rounding alpha to 0 or to 1 when needed
        if (alpha < epsilon)
            alpha = 0;
        elseif (1-alpha < epsilon)
            alpha = 1;
        end
        
        % update "param" to be returned
        param = [alpha , rho];
        
        % alpha / (1 - alpha)
        alpha_ratio = alpha / (1-alpha);
        % Elasticity of substitution (when rho ~= 0,1)
        EOS = 1 / (1 - rho); 
        
        % this loop goes one set (of prices and endowment) after another,
        % and calculates the optimal choice for each of them
        % cannot use parfor any more since x2* computation depends on x1*
        for i=1:num_obs 
            
            % Maximum quantities availables for self (X) and for the other (Y) 
            Max_X = max_vals(i,1);
            Max_Y = max_vals(i,2);
            
            % p = p1 / p2 = Max_Y / Max_X
            p = Max_Y / Max_X;
            
            
            if alpha == 0   %alpha <= 10^(-6)
                
                % Case 1  -  alpha is 0 
                
                choice_matrix(i,:) = [0, Max_Y];
                
            elseif alpha == 1   %alpha >= (1-10^(-6))
                
                % Case 2  -  alpha is 1 
                
                choice_matrix(i,:) = [Max_X, 0];
                
            elseif rho == 0   %(abs(rho) <= 10^(-6))
                
                % Case 4  -  alpha is not 0 or 1, and rho is very close to 0  
                
                choice_matrix(i,:) = [alpha * Max_X , (1 - alpha) * Max_Y]; 
                
            elseif rho == 1   %((rho <= (1 + 10^(-6))) && (rho >= (1 - 10^(-6))))
                
                % Case 5  -  alpha is not to 0 or 1, and rho is 1 
                
                if p > alpha_ratio
                    
                    choice_matrix(i,:) = [0, Max_Y];
                    
                elseif p < alpha_ratio
                    
                    choice_matrix(i,:) = [Max_X, 0];
                    
                else
                    
                    % if p equals exactly to alpha ratio, then
                    % every choice on the budget line is an optimal choice
                    % (despite the fact that allowing the subject to have
                    % some exact prices ratio in which all bundles are
                    % optimal, is "cheating" in some way, because it is
                    % unlikely that the subject's parameters exatly fit to
                    % one of the prices ratios he was exposed to, we
                    % decided not to change this, as to be consistent with
                    % our decision not to change it in CASE 5 in CRRA)
                    choice_matrix(i,:) = observed_bundle(i,:);
                    
                end
                
            elseif rho < 1   %rho < (1 - 10^(-6))
                
                % Case 6  -  alpha is not 0 or 1, and rho is
                %            not 0, and rho is strictly-smaller than 1
                % Case 3 is included in Case 6; when rho -> -Inf, Case 6
                % also provides the Leontief function (min(x,y)) as the
                % result.
                
                base = (p / alpha_ratio) ^ EOS; 
                
                choice_matrix(i,:) = [Max_Y / (p + base), Max_Y / (1 + (p / base))];  
                
            elseif rho > 1   %rho > (1 + 10^(-6))
                
                % Case 7  -  alpha is not 0 or 1, and rho is 
                %            strictly-greater than 1 
                
                if p > alpha_ratio^(1/rho)
                    
                    choice_matrix(i,:) = [0, Max_Y];
                    
                elseif p < alpha_ratio^(1/rho)
                    
                    choice_matrix(i,:) = [Max_X, 0];
                    
                else
                    
                    % when p == alpha_ratio^(1/rho),
                    % both (max_x,0) and (0,max_y) are optimal,
                    % therefore we take the one of these two that is
                    % closer to the chosen bundle
                    
                    % assuming the choice is on the budget line, if its X
                    % value is greater than max_x/2 then it is closer to
                    % the (max_x,0) corner, otherwise it is closer to the
                    % (0,max_y) corner)
                    if observed_bundle(i,1) > (Max_X / 2)
                        
                        choice_matrix(i,:) = [Max_X, 0];
                        
                    else
                        
                        choice_matrix(i,:) = [0, Max_Y];
                        
                    end
                    
                end            
                
            end
            
            % warnings to consol for debugging purposes (only when debugger mode is activated) 
            if (HPZ_Constants.debugger_mode)
                if (isnan(choice_matrix(i,1)) || isnan(choice_matrix(i,2)))
                    warning('At least one of the commodities in the optimal bundle is NaN when Rho is %.20g, Alpha is %.20g, P is %.20g and i is %d.', rho, alpha, p, i);
                elseif (isinf(choice_matrix(i,1)) || isinf(choice_matrix(i,2)))
                    warning('At least one of the commodities in the optimal bundle is +Inf or -Inf when Rho is %.20g, Alpha is %.20g, P is %.20g and i is %d.', rho, alpha, p, i);
                end
            end
            
        end   % end of loop
        
    end   % end of CES
    
end   % end of other regarding



end   % end of function

