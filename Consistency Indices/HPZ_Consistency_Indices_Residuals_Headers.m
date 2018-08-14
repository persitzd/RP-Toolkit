function col_headers = HPZ_Consistency_Indices_Residuals_Headers (GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags)

% this function returns a cell that contains all the column headers required 
% for the residuals files for consistency and inconsistency indices

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



%% step 1 : first we count how many headers/columns headers are needed
% column 1 is subject number, column 2 is observation number
num_of_columns = 2;
% counting the number of headers
if (GARP_flags(1) == 1)
    if (GARP_flags(2) == 1)
        % 1 column to print the original full value, for comparison
        num_of_columns = num_of_columns + sum(GARP_flags(5:7));
        if (GARP_flags(3) == 1)
            % 2 columns for in-sample : (1) residual value (2) residual % 
            num_of_columns = num_of_columns + 2 * sum(GARP_flags(5:7));
        end
        if (GARP_flags(4) == 1)
            % 2 columns for out-of-sample : (1) residaul value (2) residual difference % 
            num_of_columns = num_of_columns + 2 * sum(GARP_flags(5:7));
        end
    end
end
if (AFRIAT_flags(1) == 1)
    if (AFRIAT_flags(2) == 1)
        if (AFRIAT_flags(3) == 1)
            % in-sample - irrelevant (not implemented)
        end
        if (AFRIAT_flags(4) == 1)
            num_of_columns = num_of_columns + 3;
        end
    end
end
if (VARIAN_flags(1) == 1)
    if (VARIAN_flags(2) == 1)
        if (VARIAN_flags(3) == 1)
            num_of_columns = num_of_columns + 3;
        end
        if (VARIAN_flags(4) == 1)
            num_of_columns = num_of_columns + 14;
        end
    end
end
if (HOUTMAN_flags(1) == 1)
    if (HOUTMAN_flags(2) == 1)
        if (HOUTMAN_flags(3) == 1)
            % not implemented
        end
        if (HOUTMAN_flags(4) == 1)
            num_of_columns = num_of_columns + 6;
        end
    end
end
if (MPI_flags(1) == 1)
    if (MPI_flags(2) == 1)
        if (MPI_flags(3) == 1)
            % not implemented
        end
        if (MPI_flags(4) == 1)
            num_of_columns = num_of_columns + 6;
        end
    end
end



% create a row cell array for headers, now that we know what length it should be
col_headers = cell(1, num_of_columns);



%% step 2 : now we enter the headers to the cell
col_headers{1} = 'Subject';
col_headers{2} = 'Observation';
% current (next) column to put a header to
current_col = 3;
% setting all the headers
if (GARP_flags(1) == 1)
    if (GARP_flags(2) == 1)
        if (GARP_flags(5) == 1)
            col_headers{current_col} = 'WARP VIO-PAIRS full index';
            current_col = current_col + 1;
            if (GARP_flags(3) == 1)
                col_headers{current_col} = 'WARP VIO-PAIRS in-Sample';
                col_headers{current_col+1} = 'WARP VIO-PAIRS in-Sample (%)';
                current_col = current_col + 2;
            end
            if (GARP_flags(4) == 1)
                col_headers{current_col} = 'WARP VIO-PAIRS out-of-Sample index';
                col_headers{current_col+1} = 'WARP VIO-PAIRS difference (%) full from out-of-Sample';
                current_col = current_col + 2;
            end
        end
        if (GARP_flags(6) == 1)
            col_headers{current_col} = 'GARP VIO-PAIRS full index';
            current_col = current_col + 1;
            if (GARP_flags(3) == 1)
                col_headers{current_col} = 'GARP VIO-PAIRS in-Sample';
                col_headers{current_col+1} = 'GARP VIO-PAIRS in-Sample (%)';
                current_col = current_col + 2;
            end
            if (GARP_flags(4) == 1)
                col_headers{current_col} = 'GARP VIO-PAIRS out-of-Sample index';
                col_headers{current_col+1} = 'GARP VIO-PAIRS difference (%) full from out-of-Sample';
                current_col = current_col + 2;
            end
        end
        if (GARP_flags(7) == 1)
            col_headers{current_col} = 'SARP VIO-PAIRS full index';
            current_col = current_col + 1;
            if (GARP_flags(3) == 1)
                col_headers{current_col} = 'SARP VIO-PAIRS in-Sample';
                col_headers{current_col+1} = 'SARP VIO-PAIRS in-Sample (%)';
                current_col = current_col + 2;
            end
            if (GARP_flags(4) == 1)
                col_headers{current_col} = 'SARP VIO-PAIRS out-of-Sample index';
                col_headers{current_col+1} = 'SARP VIO-PAIRS difference (%) full from out-of-Sample';
                current_col = current_col + 2;
            end
        end
    end
end
if (AFRIAT_flags(1) == 1)
    if (AFRIAT_flags(2) == 1)
        if (AFRIAT_flags(3) == 1)
            % in-sample - irrelevant (not implemented)
        end
        if (AFRIAT_flags(4) == 1)
            col_headers{current_col} = 'AFRIAT full index';
            col_headers{current_col+1} = 'AFRIAT out-of-Sample index';
            col_headers{current_col+2} = 'AFRIAT difference full from out-of-Sample';
            current_col = current_col + 3;
        end
    end
end
if (VARIAN_flags(1) == 1)
    if (VARIAN_flags(2) == 1)
        if (VARIAN_flags(3) == 1)
            col_headers{current_col} = 'VARIAN Min in-Sample';
            col_headers{current_col+1} = 'VARIAN Mean in-Sample';
            col_headers{current_col+2} = 'VARIAN AVGSSQ in-Sample';
            current_col = current_col + 3;
        end
        if (VARIAN_flags(4) == 1)
            col_headers{current_col} = 'VARIAN Min full index';
            col_headers{current_col+1} = 'VARIAN Mean full index';
            col_headers{current_col+2} = 'VARIAN AVGSSQ full index';
            col_headers{current_col+3} = 'VARIAN full index is-exact';
            col_headers{current_col+4} = 'VARIAN Min out-of-Sample index';
            col_headers{current_col+5} = 'VARIAN Mean out-of-Sample index';
            col_headers{current_col+6} = 'VARIAN AVGSSQ out-of-Sample index';
            col_headers{current_col+7} = 'VARIAN out-of-Sample index is-exact';
            col_headers{current_col+8} = 'VARIAN Min difference full from out-of-Sample';
            col_headers{current_col+9} = 'VARIAN Min normalized difference full from out-of-Sample';
            col_headers{current_col+10} = 'VARIAN Mean difference full from out-of-Sample';
            col_headers{current_col+11} = 'VARIAN Mean normalized difference full from out-of-Sample';
            col_headers{current_col+12} = 'VARIAN AVGSSQ difference full from out-of-Sample';
            col_headers{current_col+13} = 'VARIAN AVGSSQ normalized difference full from out-of-Sample';
            current_col = current_col + 14;
        end
    end
end
if (HOUTMAN_flags(1) == 1)
    if (HOUTMAN_flags(2) == 1)
        if (HOUTMAN_flags(3) == 1)
            % not implemented
        end
        if (HOUTMAN_flags(4) == 1)
            col_headers{current_col} = 'HM full index';
            col_headers{current_col+1} = 'HM full index is-exact';
            col_headers{current_col+2} = 'HM residual index';
            col_headers{current_col+3} = 'HM residual index is-exact';
            col_headers{current_col+4} = 'HM difference full from residual';
            col_headers{current_col+5} = 'HM normalized difference full from residual';
            current_col = current_col + 6;
        end
    end
end
if (MPI_flags(1) == 1)
    if (MPI_flags(2) == 1)
        if (MPI_flags(3) == 1)
            % not implemented
        end
        if (MPI_flags(4) == 1)
            col_headers{current_col} = 'MPI Mean full index';
            col_headers{current_col+1} = 'MPI Median full index';
            col_headers{current_col+2} = 'MPI Mean out-of-Sample index';
            col_headers{current_col+3} = 'MPI Median out-of-Sample index';
            col_headers{current_col+4} = 'MPI Mean difference full from out-of-Sample';
            col_headers{current_col+5} = 'MPI Median difference full from out-of-Sample';
            current_col = current_col + 6; %#ok<NASGU>
        end
    end
end

end