function [WARP_GARP_SARP_Mat, col_counter] = HPZ_WARP_GARP_SARP_residuals(GARP_flags, VIO_PAIRS, WARP, GARP, SARP, varargin)

% this matrix is a helper matrix for WARP/GARP/SARP out-of-sample residuals. 
% we performed the calculations earlier and temporarily stored them in
% this matrix, and now assign them to the main matrix.
% it is used/needed only when out-of-sample is chosen.
if ~isempty(varargin) && ~isempty(varargin{1})
    Mat_GARP = varargin{1};
else
    Mat_GARP = [];  % we don't need it for in-sample
    if GARP_flags(2) && GARP_flags(4)
        error('You must send Mat_GARP as varargin when attempting out-of-sample residuals'); 
    end
end


[obs_num,~] = size(WARP);


WARP_full = VIO_PAIRS(1);
GARP_full = VIO_PAIRS(2);
SARP_full = VIO_PAIRS(3);
WARP_VIO_PAIRS = VIO_PAIRS(1);
GARP_VIO_PAIRS = VIO_PAIRS(2);
SARP_VIO_PAIRS = VIO_PAIRS(3);


    
% explanation about num of columns: for each of WARP, GARP and SARP (flags 5,6,7), 
% if chosen, we need to print its full index, as well as: if in-sample (flag 3), 
% also the residual and residual(%), and if out-of-sample (flag 4), its
% out-of-sample index value and its difference in % from the full index.
num_of_columns = (1 + 2*GARP_flags(3) + 2*GARP_flags(4)) * sum(GARP_flags(5:7));
% initialization of residuals matrix
WARP_GARP_SARP_Mat = zeros(obs_num, num_of_columns);

% initialization of column counter
col_counter = 1;



%% WARP
if (GARP_flags(5) == 1 && WARP_full ~= 0)
    % full value
    WARP_GARP_SARP_Mat(:, col_counter) = WARP_full;
    % update the column
    col_counter = col_counter + 1;

    % in sample
    if GARP_flags(3)
        WARP_Pairs = triu(WARP|(WARP'));
        % this loop go through all the observations and finds the residual for each
        for i=1:obs_num
            for j=i:obs_num
                if WARP_Pairs(i,j) == 1
                    % if there is a GARP violation involving observations
                    % i and j, increase the in-sample residual for i and j
                    WARP_GARP_SARP_Mat(i, col_counter) = WARP_GARP_SARP_Mat(i, col_counter) + 1;  
                    WARP_GARP_SARP_Mat(j, col_counter) = WARP_GARP_SARP_Mat(j, col_counter) + 1;
                end
            end
        end

        % this loop finds the residual as (%) of the total WARP
        for i=1:obs_num
            WARP_GARP_SARP_Mat(i, col_counter+1) = WARP_GARP_SARP_Mat(i, col_counter) / WARP_VIO_PAIRS;
        end

        % update the column
        col_counter = col_counter + 2;
    end

    % out of sample
    if GARP_flags(4)
        for i=1:obs_num
            WARP_GARP_SARP_Mat(i, col_counter) = Mat_GARP(i, 1);
            WARP_GARP_SARP_Mat(i, col_counter+1) = Mat_GARP(i, 2);
        end

        % update the column
        col_counter = col_counter + 2;
    end
end

%% GARP
if (GARP_flags(6) == 1 && GARP_full ~= 0)
    % full value
    WARP_GARP_SARP_Mat(:, col_counter) = GARP_full;
    % update the column
    col_counter = col_counter + 1;

    % in sample
    if GARP_flags(3)
        GARP_Pairs = triu(GARP|(GARP'));
        % this loop go through all the observations and finds the residual for each
        for i=1:obs_num
            for j=i:obs_num
                if GARP_Pairs(i,j) == 1
                    % if there is a GARP violation involving observations
                    % i and j, increase the in-sample residual for i and j
                    WARP_GARP_SARP_Mat(i, col_counter) = WARP_GARP_SARP_Mat(i, col_counter) + 1;
                    WARP_GARP_SARP_Mat(j, col_counter) = WARP_GARP_SARP_Mat(j, col_counter) + 1;
                end
            end
        end

        % this loop finds the residual as (%) of the total GARP
        for i=1:obs_num
            WARP_GARP_SARP_Mat(i, col_counter+1) = WARP_GARP_SARP_Mat(i, col_counter) / GARP_VIO_PAIRS;
        end

        % update the column
        col_counter = col_counter + 2;
    end

    % out of sample
    if GARP_flags(4)
        for i=1:obs_num
            WARP_GARP_SARP_Mat(i, col_counter) = Mat_GARP(i, 3);
            WARP_GARP_SARP_Mat(i, col_counter+1) = Mat_GARP(i, 4);
        end

        % update the column
        col_counter = col_counter + 2;
    end
end

%% SARP
if (GARP_flags(7) == 1 && SARP_full ~= 0)
    % full value
    WARP_GARP_SARP_Mat(:, col_counter) = SARP_full;
    % update the column
    col_counter = col_counter + 1;

    % in sample 
    if GARP_flags(3)
        SARP_Pairs = triu(SARP|(SARP'));
        % this loop go through all the observations and finds the residual for each
        for i=1:obs_num
            for j=i:obs_num
                if SARP_Pairs(i,j) == 1
                    % if there is a GARP violation involving observations
                    % i and j, increase the in-sample residual for i and j
                    WARP_GARP_SARP_Mat(i, col_counter) = WARP_GARP_SARP_Mat(i, col_counter) + 1;
                    WARP_GARP_SARP_Mat(j, col_counter) = WARP_GARP_SARP_Mat(j, col_counter) + 1;
                end
            end
        end

        % this loop finds the residual as (%) of the total SARP
        for i=1:obs_num
            WARP_GARP_SARP_Mat(i, col_counter+1) = WARP_GARP_SARP_Mat(i, col_counter) / SARP_VIO_PAIRS;
        end

        % update the column
        col_counter = col_counter + 2;
    end

    % out of sample
    if GARP_flags(4)
        for i=1:obs_num
            WARP_GARP_SARP_Mat(i, col_counter) = Mat_GARP(i, 5);
            WARP_GARP_SARP_Mat(i, col_counter+1) = Mat_GARP(i, 6);
        end

        % update the column
        col_counter = col_counter + 2;
    end
end


end