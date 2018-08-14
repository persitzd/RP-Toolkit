function [criterions, param] = HPZ_MMI_Analytic(param, observations, function_flag, pref_class)

% The function calculates  and returns the MMI criterion per observation. 
% Given a specific functional form and prices, we look for the lowest 
% expenditure level that yields at least the same level of utility as does 
% the observed choice.

% The function returns the value of the level of expenditure for 
% each observaion, therefore criterion is a vector with length equal to the 
% number of observations.

% "param" is returned because it may change (be rounded) during the
% calculations

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% number of observations
[obs_num,~] = size(observations);

% set prices from observation matrix
prices = observations(:,3:4);

% set observed bundle (x1,x2) from observation matrix
observed_bundle = observations(:,1:2);

% initialization of criterion matrix
MMI_values = zeros(obs_num,1);

%% risk preferences
if pref_class == HPZ_Constants.risk_pref
    
    %% CRRA
    if function_flag == HPZ_Constants.CRRA_func
        
        % When beta's value or rho's value is too close to -1 or to 1, problems 
        % may arise due to Matlab's rounding scheme (the case of rounding values 
        % that are close to zero is treated by the order of calculations). 
        % A simple example for such a problem is the following:
        % Suppose beta = (-1 + 2^-53) then beta will be assigned this number, 
        % and Matlab will perform no rounding (while beta = (-1 + 2^54)
        % will be rounded to -1). But notice that (2+beta) = (1 + 2^-53)
        % requires one more significant binary digit and therefore will be 
        % rounded to 1 while beta itself was not rounded to -1. 
        % Since in many cases, beta and 2+beta are raised to high
        % powers (e.g. when rho is high) this may cause significant
        % miscalculations.
        % These inconsistencies in Matlab's rounding scheme require a
        % threshold that should be at least: eps*64 = 2^-44, in order to
        % avoid significant miscalculations.
        % To be sure, we take an even slightly bigger threshold.
        %   (NOTE: we don't need here HPZ_Constants.print_threshold 
        %    because these threshold are anyway a lot bigger than it)
        threshold = eps * 2^12;   % 2^(-40)
        if abs(param(1) - (-1)) < threshold
            param(1) = -1;
        end
        if abs(param(2) - 1) < threshold
            param(2) = 1;
        end
        
        % it is more convenient to use these notations
        beta = param(1);
        rho  = param(2);
        
        % loop over all observations, to determine the criterion separately
        % for each observation, assuming the given values for the 
        % parameters (beta and rho) 
        for i=1:obs_num

            p1 = prices(i, 1);
            p2 = prices(i, 2);
            
            % p = Max_Y / Max_X => p = p1 / p2
            p = p1 / p2;
            
            % minimal and maximal observed quantities in this observation 
            min_x  = min (observed_bundle(i, 1), observed_bundle(i, 2));
            max_x  = max (observed_bundle(i, 1), observed_bundle(i, 2));
            
            % u is u(x_0,y_0)
            u = HPZ_Risk_Utility (0.5, [max_x, min_x], beta, rho, 0, function_flag);
            
            % u_0 is (1-rho)*(2+beta)*u(x_0,y_0)
            % (u_0 should be used only when rho ~= 1)
            u_0 = (1-rho) * (2+beta) * u;
            
            % initialization to avoid Error
            x1_opt = NaN;
            x2_opt = NaN; %#ok<*NASGU>
            
            
            % CASE 7
            if (rho == 1) && (beta >= 0)
                if  p < 1/(1+beta)
                    % compute the optimal expenditure
                    MMI_values(i,1) = (2+beta) * ( (p1^(1/(2+beta))) * (p2/(1+beta))^((1+beta)/(2+beta)) ) * exp(u);
                elseif 1/(1+beta) <= p && p <= (1+beta)
                    % compute the optimal expenditure
                    MMI_values(i,1) = (p1 + p2) * exp(u);     
                elseif p > (1+beta)
                    % compute the optimal expenditure
                    MMI_values(i,1) = (2+beta) * ( (p2^(1/(2+beta))) * (p1/(1+beta))^((1+beta)/(2+beta)) ) * exp(u);
                end
            
            % CASE 6   (note that (Rho ~= 1))
            elseif rho > 0 && beta >= 0
                
                % here in CASE 6, the order of calculation matters, because
                % an unappropriate order may result in rounding to 0 or to
                % Inf in intermediate calculations.
                % it is therefore important to have all arguments that
                % are raised by the same power (power of 1/rho, and the
                % arguments are p and (1+beta)) together before raising
                % to the power.
                % that is the reason for the formulas having a different
                % order of calculation than in the DA2 document.
                if  p < 1/(1+beta)
                    % x1* optimal:
                    x1_opt = ( u_0 / ( 1 + ( (p*(1+beta))^(1/rho) ) / p ) ) ^ (1/(1-rho));
                    % x2* optimal:
                    x2_opt = ( u_0 / ( (1+beta) + (p*(1+beta))^(-(1-rho)/rho) ) ) ^ (1/(1-rho));
                    % compute the optimal expenditure, using optimal bundle
                    % (x1*, x2*)
                    MMI_values(i,1) = (p1 * x1_opt) + (p2 * x2_opt);
                elseif 1/(1+beta) <= p && p <= (1+beta)
                    % x1* optimal:
                    x1_opt = ( u_0  / ( 2+beta ) ) ^ (1/(1-rho));
                    % x2* optimal: x1*
                    x2_opt = x1_opt;
                    % compute the optimal expenditure, using optimal bundle
                    % (x1*, x2*)
                    MMI_values(i,1) = (p1 * x1_opt) + (p2 * x2_opt);
                elseif p > (1+beta)
                    % x1* optimal:
                    x1_opt = ( u_0 / ( (1+beta) + ( p / (1+beta) )^( (1-rho)/rho ) ) ) ^ (1/(1-rho));
                    % x2* optimal:
                    x2_opt = ( u_0 / ( 1 + ((1+beta)/p)^(1/rho) * p ) ) ^ (1/(1-rho));
                    % compute the optimal expenditure, using optimal bundle
                    % (x1*, x2*)
                    MMI_values(i,1) = (p1 * x1_opt) + (p2 * x2_opt);
                end
                
                
            % CASE 9
            elseif (rho == 1) && (beta > -1 && beta < 0)
                
                % here in CASE 9, the order of calculation matters, because
                % an unappropriate order may result in rounding to 0 or to
                % Inf in intermediate calculations.
                % it is therefore important to have all arguments that
                % are raised by the same power (power of [(1+beta)/(2+beta)],
                % and the arguments are (1+beta) and either p1 or p2) 
                % together before raising to the power.
                % that is the reason for the formulas having a different
                % order of calculation than in the DA2 document.
                if  p < 1
                    % compute the optimal expenditure
                    MMI_values(i,1) = (2+beta) * ( (p1^(1/(2+beta))) * (p2/(1+beta))^((1+beta)/(2+beta)) ) * exp(u);
                elseif p == 1
                    % compute the optimal expenditure
                    MMI_values(i,1) = (2+beta) * p1 * ((1 / (1 + beta)) ^ ((1+beta)/(2+beta))) * exp(u);
                elseif p > 1
                    % compute the optimal expenditure
                    MMI_values(i,1) = (2+beta) * ( (p2^(1/(2+beta))) * (p1/(1+beta))^((1+beta)/(2+beta)) ) * exp(u);
                end
                
            % CASE 8   (note that (Rho ~= 1))
            elseif rho > 0 && (beta > -1 && beta < 0)
                
                % here in CASE 8, the order of calculation matters, because
                % an unappropriate order may result in rounding to 0 or to
                % Inf in intermediate calculations.
                % it is therefore important to have all arguments that
                % are raised by the same power (power of 1/rho, and the
                % arguments are p and (1+beta)) together before raising
                % to the power.
                % that is the reason for the formulas having a different
                % order of calculation than in the DA2 document.
                if  p < 1
                    % x1* optimal:
                    x1_opt = ( u_0 / ( 1 + ( (p*(1+beta))^(1/rho) ) / p ) ) ^ (1/(1-rho));
                    % x2* optimal:
                    x2_opt = ( u_0 / ( (1+beta) + (p*(1+beta))^(-(1-rho)/rho) ) ) ^ (1/(1-rho));
                    % compute the optimal expenditure, using optimal bundle
                    % (x1*, x2*)
                    MMI_values(i,1) = (p1 * x1_opt) + (p2 * x2_opt);
                elseif p == 1
                    % x1* optimal:
                    x1_opt = ( u_0 / ( 1 + ( (1+beta)^(1/rho) ) ) ) ^ (1/(1-rho));
                    % x2* optimal:
                    x2_opt = ( u_0 / ( (1+beta) + (1+beta)^(-(1-rho)/rho) ) ) ^ (1/(1-rho));
                    % compute the optimal expenditure, using optimal bundle
                    % (x1*, x2*)
                    MMI_values(i,1) = (p1 * x1_opt) + (p2 * x2_opt);        
                elseif p > 1
                    % x1* optimal:
                    x1_opt = ( u_0 / ( (1+beta) + ( p / (1+beta) )^( (1-rho)/rho ) ) ) ^ (1/(1-rho));
                    % x2* optimal:
                    x2_opt = ( u_0 / ( 1 + ((1+beta)/p)^(1/rho) * p ) ) ^ (1/(1-rho));
                    % compute the optimal expenditure, using optimal bundle
                    % (x1*, x2*)
                    MMI_values(i,1) = (p1 * x1_opt) + (p2 * x2_opt);
                end
           
            % CASE 11
            elseif (rho == 1) && (beta == -1)
                MMI_values(i,1) = min (p1, p2) * exp(u);
                
            % CASE 10   (note that (Rho ~= 1))
            elseif (rho > 0) && (beta == -1)
                % note that 2+beta = 1, hence u_0 = (1-rho)*u 
                %MMI_values(i,1) = min (p1, p2) * u_0 ^ (1/(1-rho));
                % note that when beta = -1, u = 1/(1-rho) * max(x0,y0)^(1-rho),  
                % hence: u_0 ^ (1/(1-rho)) = [(1-rho)*u] ^ (1/(1-rho)) =
                % = [(1-rho)* 1/(1-rho) * max(x0,y0)^(1-rho)] ^ (1/(1-rho)) = 
                % = [max(x0,y0)^(1-rho)] ^ (1/(1-rho)) = max(x0,y0)
                MMI_values(i,1) = min (p1, p2) * max_x;
           
            % CASE 12
            elseif (rho == 0) && (beta >= 0)
                if  p < 1/(1+beta)
                    % note that rho = 0, hence u = (max{x,y} +(1+beta)min{x,y})/(2+beta),
                    % therefore (2+beta)*u = (max{x,y} +(1+beta)min{x,y})
                    x1_opt = max_x + (1+beta)*min_x;
                    x2_opt = 0;
                    % compute the optimal expenditure, using optimal bundle
                    % (x1*, x2*)
                    MMI_values(i,1) = (p1 * x1_opt) + (p2 * x2_opt);
                elseif p > (1+beta)
                    x1_opt = 0;
                    % note that rho = 0, hence u = (max{x,y} +(1+beta)min{x,y})/(2+beta),
                    % therefore (2+beta)*u = (max{x,y} +(1+beta)min{x,y})
                    x2_opt = max_x + (1+beta)*min_x;
                    % compute the optimal expenditure, using optimal bundle
                    % (x1*, x2*)
                    MMI_values(i,1) = (p1 * x1_opt) + (p2 * x2_opt);
                elseif 1/(1+beta) <= p && p <= (1+beta)
                    % note that rho = 0, hence u = (max{x,y} +(1+beta)min{x,y})/(2+beta)
                    MMI_values(i,1) = ( ( p1 + p2 ) * ( max_x + (1+beta)*min_x ) ) / (2+beta);
                end
            
            % CASE 13
            elseif rho == 0 && (beta >= -1) && (beta < 0)
                % note that 1-rho = 1, hence u_0 = (2+beta)*u 
                MMI_values(i,1) = min (p1, p2) * u_0 ;
            end
            
            % warnings to consol for debugging purposes (only when debugger mode is activated) 
            if (HPZ_Constants.debugger_mode)
                if isnan(MMI_values(i,1))
                    warning('Criterion is NaN when Rho is %.20g, Beta is %.20g, 1+Beta is %.20g, P is %.20g, u is %.20g, u_0 is %.20g, max_x is %.20g, min_x is %.20g and i is %d.', rho, beta, 1+beta, p, u, u_0, max_x, min_x, i);
                elseif isinf(MMI_values(i,1))
                    % note that the MMI Value is the opposite of the
                    % criterion: Criterion = 1 - MMI_Value
                    if (MMI_values(i,1) > 0)
                        inf_str = '-';
                    else
                        inf_str = '+';
                    end
                    warning('Criterion is %sInf when Rho is %.20g, Beta is %.20g, 1+Beta is %.20g, P is %.20g, u is %.20g, u_0 is %.20g, max_x is %.20g, min_x is %.20g and i is %d.', inf_str, rho, beta, 1+beta, p, u, u_0, max_x, min_x, i);
                elseif (MMI_values(i,1) < -HPZ_Constants.MMI_threshold) || (MMI_values(i,1) > 1+HPZ_Constants.MMI_threshold)
                    warning('Criterion is not in the range [0,1]: Criterion is %.20g when Rho is %.20g, Beta is %.20g, 1+Beta is %.20g, P is %.20g, u is %.20g, u_0 is %.20g, max_x is %.20g, min_x is %.20g and i is %d.', 1-MMI_values(i,1), rho, beta, 1+beta, p, u, u_0, max_x, min_x, i);
                end
            end
            
        end   % end of loop over all observations of the subject 
        
    %% CARA
    elseif function_flag == HPZ_Constants.CARA_func
        
        % if beta's value is too close to -1, it may cause problems, 
        % because the accuracy of matlab is lost in values that have too 
        % many meaningful digits
        % (that is also why it is not needed when rho is close to 0 - there
        % aren't so many meaningful digits because matlab uses E-100 and such)
        threshold = eps * 2^12;   % 2^(-40)
        if abs(param(1) - (-1)) < threshold
            param(1) = -1;
        end
        
        % threshold for the CARA parameter (A). 
        % Cases 27 and 28 in the DA-2 Document refer to the case where A approaches 0. 
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
        
        % it is more convenient to use these notations
        beta = param(1);
        A = param(2);
        
        % loop over all observations, to determine the criterion separately
        % for each observation, assuming the given values for the 
        % parameters (beta and A) 
        for i=1:obs_num
            
            min_x  = min (observed_bundle(i, 1), observed_bundle(i, 2));
            max_x  = max (observed_bundle(i, 1), observed_bundle(i, 2));
            
            p1 = prices(i,1);   % price of x
            p2 = prices(i,2);   % price of y
            p = p1 / p2;
            
            % u is u(x_0,y_0)
            u = HPZ_Risk_Utility (0.5, [max_x, min_x], beta, 0, A, function_flag);
            
            % log_u is -ln(-u(x_0,y_0))            
            % we calculate -ln(-u) = ln(-1/u) directly, because otherwise the utility 
            % might be rounded to -0 (if e^(A*min(x)) is too small), then -1/u 
            % will be calculated to be Inf, and then ln(-1/u) will be calculated as Inf, 
            % (when 0<gamma<1 it can be -ln(1-gamma)+A*min(x) at most)
            % HPZ_Risk_Log_Utility provides the direct calculation of -ln(-u) 
            % to avoid such computational calculation problems.
            log_u = HPZ_Risk_Log_Utility (0.5, [max_x, min_x], beta, 0, A, function_flag);
            
            % u_0 is (2+beta)*u(x_0,y_0)
            %u_0 = -exp(-A*max_x) - (1+beta)*exp(-A*min_x);   % Old Code
            u_0 = (2+beta) * u;
            
            
            %% non-negative beta
            
            % CASE 29
            if  A == 0 
                % this case is needed because the optimization process is
                % limited by a closed interval, therefore the interval for
                % A includes 0. To avoid this case we assign it the worst
                % possible criterion (0), so A will never converge to this value.
                MMI_values(i,1) = 0;
                
            % CASE 27
            % See comment where A_threshold is defined
            elseif beta >= 0 && A > 0 && ( (min_x > 0 && A*min_x < A_threshold) || (min_x == 0 && A*max_x < A_threshold) )
                % w(x) = -e^(-ax)
                % when a->0 : w(x)=ax-1. 
                % therefore we can simplify the formula in the DA-2 document as follows:
                % (1+u) = 1 + [gamma*(a*max(x)-1)+(1-gamma)*(a*min(x)-1)] =
                %       = 1 + [gamma*a*max(x)+(1-gamma)*a*min(x) -1]
                %       = gamma*a*max(x)+(1-gamma)*a*min(x)
                % therefore (since gamma = 1/(2+beta)):
                % (1/a)*(2+beta)*(1+u) = max(x) + (1+beta)*min(x)
                
                % u_0_0 = (1/a)*(2+beta)*(1+u)
                u_0_0 = max_x + (1+beta)*min_x;
                
                if p <= 1/(1+beta)
                    MMI_values(i,1) = p1 * u_0_0;
                elseif p <= 1+beta && p >= 1/(1+beta)
                    MMI_values(i,1) = ((p1 + p2) / (2+beta)) * u_0_0;
                elseif p >= 1+beta
                    MMI_values(i,1) = p2 * u_0_0;
                end
                
            % CASE 21
            elseif beta >= 0 && -u_0 > 1+beta
                if p < -((u_0 + (1+beta)) / (1+beta))
                    MMI_values(i,1) = - ( p1 * log(-(u_0 + (1+beta))) ) / A;
                elseif p >= -((u_0 + (1+beta)) / (1+beta)) && p < 1/(1+beta)
                    % criterion(i,1) =
                    % = p1 * ( log(-( (1+p)/(p*u_0) )) / A )  +  p2 * ( log(- ( ((1+p)*(1+beta))/u_0) ) / A ) = 
                    % = p1 * ( log(-( (1+p)/(p*(2+beta)*u) )) / A )  +  p2 * ( log(- (((1+p)*(1+beta))/((2+beta)*u) ) ) / A ) =
                    % = p1 * ( ( log( (1+p)/(p*(2+beta)) ) + log(-1/u) ) / A )  +  p2 * ( ( log( ((1+p)*(1+beta))/(2+beta) ) + log(-1/u) ) / A ) =
                    % = p1 * ( ( log( (1+p)/(p*(2+beta)) ) + log_u ) / A )  +  p2 * ( ( log( ((1+p)*(1+beta))/(2+beta) ) + log_u ) / A )
                    MMI_values(i,1) = p1 * ( ( log( (1+p)/(p*(2+beta)) ) + log_u ) / A ) + ...
                                     p2 * ( ( log( ((1+p)*(1+beta))/(2+beta) ) + log_u ) / A );  
                elseif p >= 1/(1+beta) && p <= (1+beta)
                    % criterion(i,1) = 
                    % = (p1 + p2) * (log(-(2+beta)/u_0) / A) =
                    % = (p1 + p2) * (log(1/(-u)) / A) =
                    % = (p1 + p2) * ((-1)*(log(-u) / A) =
                    % = (p1 + p2) * log_u / A
                    MMI_values(i,1) = (p1 + p2) * (log_u / A);
                elseif p > 1+beta && p <= - ((1+beta)/(u_0 + (1+beta)))
                    % criterion(i,1) = 
                    % = p1 * ( log(- ( ((1+p)*(1+beta))/(p*u_0) ) ) / A )  +  p2 * ( log(- ( (1+p)/(u_0) )) / A ) = 
                    % = p1 * ( log(- ( ((1+p)*(1+beta))/(p*(2+beta)*u) ) ) / A )  +  p2 * ( log(- ( (1+p)/((2+beta)*u) )) / A ) = 
                    % = p1 * ( ( log(- ( ((1+p)*(1+beta))/(p*(2+beta)) ) ) + log(-1/u) ) / A )  +  p2 * ( ( log(- ( (1+p)/((2+beta)) )) + log(-1/u) ) / A ) =
                    % = p1 * ( ( log( ((1+p)*(1+beta))/(p*(2+beta)) ) + log_u ) / A )  +  p2 * ( ( log( (p+1)/((2+beta)) ) + log_u ) / A )
                    MMI_values(i,1) = p1 * ( ( log( ((1+p)*(1+beta))/(p*(2+beta)) ) + log_u ) / A ) + ...
                                     p2 * ( ( log( (p+1)/((2+beta)) ) + log_u ) / A );   
                elseif p > - ((1+beta)/(u_0 + (1+beta)))
                    MMI_values(i,1) = - ( p2 * log(-(u_0 + (1+beta))) ) / A;
                end
                
            % CASE 22    
            elseif beta >=0 && -u_0 <= 1+beta
                if p < 1/ (1+beta)
                    % criterion(i,1) = 
                    % = p1 * ( log(-( (1+p)/(p*u_0) )) / A )  +  p2 * (log(- ( ((1+p)*(1+beta))/u_0 ) ) / A ) = 
                    % = p1 * ( log(-( (1+p)/(p*(2+beta)*u) )) / A )  +  p2 * ( log(- ( ((1+p)*(1+beta))/((2+beta)*u) ) ) / A ) = 
                    % = p1 * ( ( log( (1+p)/(p*(2+beta)) ) + log(-1/u) ) / A )  +  p2 * ( ( log( ((1+p)*(1+beta))/(2+beta) ) + log(-1/u) ) / A ) = 
                    % = p1 * ( ( log( (1+p)/(p*(2+beta)) ) + log_u ) / A )  +  p2 * ( ( log( ((1+p)*(1+beta))/(2+beta) ) + log_u ) / A ) 
                    MMI_values(i,1) = p1 * ( ( log( (1+p)/(p*(2+beta)) ) + log_u ) / A ) + ...
                                     p2 * ( ( log( ((1+p)*(1+beta))/(2+beta) ) + log_u ) / A );
                elseif p >= 1/(1+beta) && p <= 1+beta
                    % criterion(i,1) = 
                    % = (p1 + p2) * (log(-(2+beta)/u_0) / A) =
                    % = (p1 + p2) * (log(1/(-u)) / A) =
                    % = (p1 + p2) * (-1)*(log(-u) / A) =
                    % = (p1 + p2) * log_u / A
                    MMI_values(i,1) = (p1 + p2) * log_u / A;
                elseif p > 1 + beta
                    % criterion(i,1) = 
                    % = p1 * ( log(- ( ((1+p)*(1+beta))/(p*u_0) ) ) / A )  +  p2 * ( log(- ( (1+p)/(u_0) )) / A ) = 
                    % = p1 * ( log(- ( ((1+p)*(1+beta))/(p*(2+beta)*u) ) ) / A )  +  p2 * ( log(- ( (1+p)/((2+beta)*u) )) / A ) = 
                    % = p1 * ( ( log( ((1+p)*(1+beta))/(p*(2+beta)) ) + log(-1/u) ) / A )  +  p2 * ( ( log( (1+p)/((2+beta)) ) + log(-1/u) ) / A ) = 
                    % = p1 * ( ( log( ((1+p)*(1+beta))/(p*(2+beta)) ) + log_u ) / A )  +  p2 * ( ( log( (1+p)/((2+beta)) ) + log_u ) / A )
                    MMI_values(i,1) = p1 * ( ( log( ((1+p)*(1+beta))/(p*(2+beta)) ) + log_u ) / A ) + ...
                                     p2 * ( ( log( (1+p)/((2+beta)) ) + log_u ) / A );
                end
                
            %% negative beta
            
            % beta equals to -1 : CASE 26 and the extreme of 28
            elseif beta == -1
                % We can simplify the formula in the DA-2 document as follows:                
                % when beta=-1 then u_0=-e^-(A*max(x)),
                % and then ln(-1/u_0)=A*max(x)
                % and 1/A times ln(-1/u_0) is exactly max(x).
                if p <= 1
                    MMI_values(i,1) = p1 * max_x;
                else
                    MMI_values(i,1) = p2 * max_x;
                end
                
            % CASE 28 
            % See comment where A_threshold is defined
            elseif beta >= -1 && beta < 0 && A > 0 && ( (min_x > 0 && A*min_x < A_threshold) || (min_x == 0 && A*max_x < A_threshold) )
                % w(x) = -e^(-ax)
                % when a->0 : w(x)=ax-1. 
                % therefore we can simplify the formula in the DA-2 document as follows:
                % (1+u) = 1 + [gamma*(a*max(x)-1)+(1-gamma)*(a*min(x)-1)] =
                %       = 1 + [gamma*a*max(x)+(1-gamma)*a*min(x) -1]
                %       = gamma*a*max(x)+(1-gamma)*a*min(x)
                % therefore (since gamma = 1/(2+beta)):
                % (1/a)*(2+beta)*(1+u) = max(x) + (1+beta)*min(x)
                
                % u_0_0 = (1/a)*(2+beta)*(1+u)
                u_0_0 = max_x + (1+beta)*min_x; 
                MMI_values(i,1) = min(p1, p2) * u_0_0;
            
            % CASES 23-24
            elseif beta > -1 && beta < 0  && -u_0 > 1+beta
                
                % CASE 23
                if - ((u_0 + (1+beta))/(1+beta)) > 1
                    if p <= 1   % p_x <= p_y
                        MMI_values(i,1) = - (p1 * log(-(u_0 + (1+beta)) )) / A;
                    else        % p_x > p_y
                        MMI_values(i,1) = - (p2 * log(-(u_0 + (1+beta)) )) / A;
                    end
                    
                % CASE 24    
                else
                    if p < - ((u_0 + (1+beta))/(1+beta))
                        MMI_values(i,1) = - ( p1 * log(-(u_0 + (1+beta))) ) / A;
                    elseif p >= - ((u_0 + (1+beta))/(1+beta)) && p <= 1
                        % criterion(i,1) = 
                        % = p1 * ( log(-( (1+p)/(p*u_0) )) / A )  +  p2 * ( log(- ( ((1+p)*(1+beta))/u_0 ) ) / A ) = 
                        % = p1 * ( log(-( (1+p)/(p*(2+beta)*u) )) / A )  +  p2 * ( log(- ( ((1+p)*(1+beta))/((2+beta)*u) ) ) / A ) = 
                        % = p1 * ( ( log( (1+p)/(p*(2+beta)) ) + log(-1/u) ) / A )  +  p2 * ( ( log( ((1+p)*(1+beta))/(2+beta) ) + log(-1/u) ) / A ) = 
                        % = p1 * ( ( log( (1+p)/(p*(2+beta)) ) + log_u ) / A )  +  p2 * ( ( log( ((1+p)*(1+beta))/(2+beta) ) + log_u ) / A )
                        MMI_values(i,1) = p1 * ( ( log( (1+p)/(p*(2+beta)) ) + log_u ) / A ) + ...
                                         p2 * ( ( log( ((1+p)*(1+beta))/(2+beta) ) + log_u ) / A );    
                    elseif p > 1 && p <= - ((1+beta)/(u_0 + (1+beta)))
                        % criterion(i,1) =
                        % = p1 * ( log(- ( ((1+p)*(1+beta))/(p*u_0) ) ) / A )  +  p2 * ( log(- ( (1+p)/(u_0) )) / A ) = 
                        % = p1 * ( log(- ( ((1+p)*(1+beta))/(p*(2+beta)*u) ) ) / A )  +  p2 * ( log(- ( (1+p)/((2+beta)*u) )) / A ) = 
                        % = p1 * ( ( log( ((1+p)*(1+beta))/(p*(2+beta)) ) + log(-1/u) ) / A )  +  p2 * ( ( log( (1+p)/((2+beta)) ) + log(-1/u) ) / A ) = 
                        % = p1 * ( ( log( ((1+p)*(1+beta))/(p*(2+beta)) ) + log_u ) / A )  +  p2 * ( ( log( (1+p)/((2+beta)) ) + log_u ) / A ) 
                        MMI_values(i,1) = p1 * ( ( log( ((1+p)*(1+beta))/(p*(2+beta)) ) + log_u ) / A ) + ...
                                         p2 * ( ( log( (1+p)/((2+beta)) ) + log_u ) / A );    
                    elseif p > - ((1+beta)/(u_0 + (1+beta)))
                        MMI_values(i,1) = - ( p2 * log(-(u_0 + (1+beta))) ) / A;
                    end 
                end
            
            % CASE 25 
            elseif beta > -1 && beta < 0  && -u_0 <= 1+beta
                if p <= 1
                    % criterion(i,1) = 
                    % = p1 * ( log(-( (1+p)/(p*u_0) )) / A )  +  p2 * ( log(- ( ((1+p)*(1+beta))/u_0 ) ) / A ) = 
                    % = p1 * ( log(-( (1+p)/(p*(2+beta)*u) )) / A )  +  p2 * ( log(- ( ((1+p)*(1+beta))/((2+beta)*u) ) ) / A ) = 
                    % = p1 * ( ( log( (1+p)/(p*(2+beta)) ) + log(-1/u) ) / A )  +  p2 * ( ( log( ((1+p)*(1+beta))/(2+beta) ) + log(-1/u) ) / A ) = 
                    % = p1 * ( ( log( (1+p)/(p*(2+beta)) ) + log_u ) / A )  +  p2 * ( ( log( ((1+p)*(1+beta))/(2+beta) ) + log_u ) / A )
                    MMI_values(i,1) = p1 * ( ( log((1+p)/(p*(2+beta)) ) + log_u ) / A ) + ...
                                     p2 * ( ( log(((1+p)*(1+beta))/(2+beta) ) + log_u ) / A );
                else    
                    % criterion(i,1) =
                    % = p1 * ( log(- ( ((1+p)*(1+beta))/(p*u_0) ) ) / A )  +  p2 * ( log(- ( (1+p)/(u_0) )) / A ) = 
                    % = p1 * ( log(- ( ((1+p)*(1+beta))/(p*(2+beta)*u) ) ) / A )  +  p2 * ( log(- ( (1+p)/((2+beta)*u) )) / A ) = 
                    % = p1 * ( ( log( ((1+p)*(1+beta))/(p*(2+beta)) ) + log(-1/u) ) / A )  +  p2 * ( ( log( (1+p)/((2+beta)) ) + log(-1/u) ) / A ) = 
                    % = p1 * ( ( log( ((1+p)*(1+beta))/(p*(2+beta)) ) + log_u ) / A )  +  p2 * ( ( log( (1+p)/((2+beta)) ) + log_u ) / A )   
                    MMI_values(i,1) = p1 * ( ( log( ((1+p)*(1+beta))/(p*(2+beta)) ) + log_u ) / A ) + ...
                                     p2 * ( ( log( (1+p)/((2+beta)) ) + log_u ) / A );
                end
            end
            
            % warnings to consol for debugging purposes (only when debugger mode is activated) 
            if (HPZ_Constants.debugger_mode)
                if isnan(MMI_values(i,1))
                    warning('Criterion is NaN when A is %.20g, Beta is %.20g, 1+Beta is %.20g, P is %.20g, u is %.20g, u_0 is %.20g, log_u is %.20g, max_x is %.20g, min_x is %.20g and i is %d.', A, beta, 1+beta, p, u, u_0, log_u, max_x, min_x, i);
                elseif isinf(MMI_values(i,1))
                    % note that the MMI Value is the opposite of the
                    % criterion: Criterion = 1 - MMI_Value
                    if (MMI_values(i,1) > 0)
                        inf_str = '-';
                    else
                        inf_str = '+';
                    end
                    warning('Criterion is %sInf when A is %.20g, Beta is %.20g, 1+Beta is %.20g, P is %.20g, u is %.20g, u_0 is %.20g, log_u is %.20g, max_x is %.20g, min_x is %.20g and i is %d.', inf_str, A, beta, 1+beta, p, u, u_0, log_u, max_x, min_x, i);
                elseif (MMI_values(i,1) < -HPZ_Constants.MMI_threshold) || (MMI_values(i,1) > 1+HPZ_Constants.MMI_threshold)
                    warning('Criterion is not in the range [0,1]: Criterion is %.20g when A is %.20g, Beta is %.20g, 1+Beta is %.20g, P is %.20g, u is %.20g, u_0 is %.20g, log_u is %.20g, max_x is %.20g, min_x is %.20g and i is %d.', 1-MMI_values(i,1), A, beta, 1+beta, p, u, u_0, log_u, max_x, min_x, i);
                end
            end
            
        end   % end of loop over all observations of the subject 
        
    end   % end of CRRA and CARA
    
%% other-regarding preferences
elseif pref_class == HPZ_Constants.OR_pref
    
    %% CES
    if function_flag == HPZ_Constants.CES_func
        
        % it is more convenient to use these notations
        alpha = param(1);
        rho = param(2);
        
        % this epsilon is the threshold for deciding when a value of alpha
        % or rho is close enough to 0 or to 1, to be considered as 0 or 1,
        % respectively.
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
        
        
        % loop over all observations, to determine the criterion separately
        % for each observation, assuming the given values for the 
        % parameters (alpha and rho) 
        for i = 1 : obs_num
            
            p_x = prices(i,1);
            p_y = prices(i,2);
            p = p_x / p_y;
            x_0 = observed_bundle(i, 1);
            y_0 = observed_bundle(i, 2);
            
            % CASE 8
            if alpha == 0              
                % when alpha == 0, the subject cares only about the other
                MMI_values(i,1) = p_y * y_0 ;
                
            % CASE 9    
            elseif alpha == 1        
                % when alpha == 1, the subject cares only about himself  
                MMI_values(i,1) = p_x * x_0 ;
                
            % CASE 11    
            elseif rho == 0
                MMI_values(i,1) = (((p_x * x_0)/alpha)^alpha) * (((p_y * y_0)/(1-alpha))^(1-alpha)); 
                
            % CASE 12    
            elseif rho == 1
                if alpha_ratio < p
                    MMI_values(i,1) = p_y * (y_0 + (alpha_ratio * x_0));                
                elseif alpha_ratio == p
                    % If the choice is on the budget line (as we assume), 
                    % (p_x * x_0) + (p_y * y_0) will always be equal to 1.
                    MMI_values(i,1) = 1;                                
                elseif alpha_ratio > p
                    MMI_values(i,1) = p_x * (x_0 + ((1 / alpha_ratio) * y_0));                
                end    
                
            % CASE 13 (includes CASE 10 as well)
            elseif (rho < 1)
                % the utility of the subject given the being-tested
                % parameters and the subject's choice
                u = HPZ_OR_Utility ([x_0,y_0], alpha, rho, function_flag);
                if (rho > -1)
                    % if rho is in (-1,0)U(0,1), we need to worry about rho
                    % being too close to 1 since then related expressions
                    % may go either to 0 or to Inf.
                    % To avoid this problem we multiply and divide all elements 
                    % before raising them to the power of 1/(rho-1).
                    den_Part_1 = ( alpha + (((alpha/p)^rho)/(1-alpha))^(1/(rho-1)) )^(1/rho);
                    den_Part_2 = ( (1-alpha) + ((((1-alpha)*p)^rho)/alpha)^(1/(rho-1)) )^(1/rho);
                else
                    % if rho is in (-Inf,-1), we need to worry about rho
                    % being too close to -INF since then related expressions
                    % may go either to 0 or to Inf.
                    % To avoid this problem we multiply and divide all elements 
                    % before raising them to the power of rho/(rho-1).
                    rho_ratio = rho / (rho - 1);
                    den_Part_1 = ( alpha + ((alpha/p)^rho_ratio)/((1-alpha)^(1/(rho-1))) )^(1/rho);
                    den_Part_2 = ( (1-alpha) + (((1-alpha)*p)^rho_ratio)/(alpha^(1/(rho-1))) )^(1/rho);
                end
                MMI_values(i,1) = u * ((p_x / den_Part_1) + (p_y / den_Part_2)); 
                
            % CASE 14
            elseif (rho > 1)
                if alpha_ratio^(1/rho) <= p
                    if (y_0 > x_0)
                        % we divide and multiply by y_0
                        MMI_values(i,1) = p_y * y_0 * ( ((alpha/(1 - alpha) * ((x_0/y_0) ^ rho)) + 1) ) ^ (1 / rho);
                    else
                        % we divide and multiply by x_0
                        MMI_values(i,1) = p_y * x_0 * ( alpha/(1 - alpha) + ((y_0/x_0) ^ rho) ) ^ (1 / rho);
                    end
                else
                    if (y_0 > x_0)
                        % we divide and multiply by y_0
                        MMI_values(i,1) = p_x * y_0 * ( (((x_0/y_0) ^ rho)) + (1 - alpha)/alpha ) ^ (1 / rho);
                    else
                        % we divide and multiply by x_0
                        MMI_values(i,1) = p_x * x_0 * ( ((1 + (1 - alpha)/alpha * ((y_0/x_0) ^ rho))) ) ^ (1 / rho);
                    end
                end
            end
            
            % warnings to consol for debugging purposes (only when debugger mode is activated) 
            if (HPZ_Constants.debugger_mode)
                if isnan(MMI_values(i,1))
                    warning('Criterion is NaN when Rho is %.20g, Alpha is %.20g, P is %.20g, u is %.20g, x_0 is %.20g, y_0 is %.20g and i is %d.', rho, alpha, p, u, x_0, y_0, i);
                elseif isinf(MMI_values(i,1))
                    % note that the MMI Value is the opposite of the
                    % criterion: Criterion = 1 - MMI_Value
                    if (MMI_values(i,1) > 0)
                        inf_str = '-';
                    else
                        inf_str = '+';
                    end
                    warning('Criterion is %sInf when Rho is %.20g, Alpha is %.20g, P is %.20g, u is %.20g, x_0 is %.20g, y_0 is %.20g and i is %d.', inf_str, rho, alpha, p, u, x_0, y_0, i);
                elseif (MMI_values(i,1) < -HPZ_Constants.MMI_threshold) || (MMI_values(i,1) > 1+HPZ_Constants.MMI_threshold)
                    warning('Criterion is not in the range [0,1]: Criterion is %.20g when Rho is %.20g, Alpha is %.20g, P is %.20g, u is %.20g, x_0 is %.20g, y_0 is %.20g and i is %d.', 1-MMI_values(i,1), rho, alpha, p, u, x_0, y_0, i);
                end
            end
            
        end   % end of loop over all observations of the subject 
        
    end   % end of CES
    
end   % end of other regarding



% MMI_values is a vector of minimal endowment (required to achieve 
% the same utility) for each observation, therefore for each observation
% the criterion is the real endowment (1) minus that minimal endowment

% in the function so far, we calculated the minimal endowment required to
% achieve the same utility, e.g. if the subject can't get this utiity with
% a smaller endowment, then MMI_values(i) == 1. but our criterion is one
% that is set to 0 when the subject chose in accordance to its preferences,
% therefore: criterions(i) = 1 - MMI_values(i)
criterions = ones(obs_num,1) - MMI_values;



end   % end of function