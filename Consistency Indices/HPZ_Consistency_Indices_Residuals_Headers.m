function [col_headers , num_of_columns] = HPZ_Consistency_Indices_Residuals_Headers (GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags)

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
        num_of_columns = num_of_columns + GARP_flags(5) + 2*GARP_flags(6) + GARP_flags(7);
        if (GARP_flags(3) == 1)
            % 2 columns for in-sample : (1) residual value (2) residual % 
            num_of_columns = num_of_columns + 2*(GARP_flags(5) + 2*GARP_flags(6) + GARP_flags(7));
        end
        if (GARP_flags(4) == 1)
            % 2 columns for out-of-sample : (1) residaul value (2) residual difference % 
            num_of_columns = num_of_columns + 2*(GARP_flags(5) + 2*GARP_flags(6) + GARP_flags(7));
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
            num_of_columns = num_of_columns + sum(VARIAN_flags(5:7))*6; % = aggregates*(4+2)  = (num aggreagtes) * [(exact, lower, approximate, upper)*(Component) + (exact, approximate)*(Difference)]   
        end
        if (VARIAN_flags(4) == 1)
            num_of_columns = num_of_columns + sum(VARIAN_flags(5:7))*12; % aggregates*(4*2 + 2*2) = (num aggreagtes) * [(exact, lower, approximate, upper)*(full + out-of-sample) + (exact, approximate)*(residual + normalized-residual)] 
        end
    end
end
if (HOUTMAN_flags(1) == 1)
    if (HOUTMAN_flags(2) == 1)
        if (HOUTMAN_flags(3) == 1)
            % not implemented
        end
        if (HOUTMAN_flags(4) == 1)
            num_of_columns = num_of_columns + 4;
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
            col_headers{current_col} = 'WARP VIOLATIONS full index';
            current_col = current_col + 1;
            if (GARP_flags(3) == 1)
                col_headers{current_col} = 'WARP VIOLATIONS in-Sample';
                col_headers{current_col+1} = 'WARP VIOLATIONS in-Sample (%)';
                current_col = current_col + 2;
            end
            if (GARP_flags(4) == 1)
                col_headers{current_col} = 'WARP VIOLATIONS out-of-Sample index';
                col_headers{current_col+1} = 'WARP VIOLATIONS difference (%) full from out-of-Sample';
                current_col = current_col + 2;
            end
        end
        
        if (GARP_flags(6) == 1)
            col_headers{current_col} = 'GARP VIOLATIONS full index';
            current_col = current_col + 1;
            if (GARP_flags(3) == 1)
                col_headers{current_col} = 'GARP VIOLATIONS in-Sample';
                col_headers{current_col+1} = 'GARP VIOLATIONS in-Sample (%)';
                current_col = current_col + 2;
            end
            if (GARP_flags(4) == 1)
                col_headers{current_col} = 'GARP VIOLATIONS out-of-Sample index';
                col_headers{current_col+1} = 'GARP VIOLATIONS difference (%) full from out-of-Sample';
                current_col = current_col + 2;
            end
            col_headers{current_col} = 'GARP PAIRS full index';
            current_col = current_col + 1;
            if (GARP_flags(3) == 1)
                col_headers{current_col} = 'GARP PAIRS in-Sample';
                col_headers{current_col+1} = 'GARP PAIRS in-Sample (%)';
                current_col = current_col + 2;
            end
            if (GARP_flags(4) == 1)
                col_headers{current_col} = 'GARP PAIRS out-of-Sample index';
                col_headers{current_col+1} = 'GARP PAIRS difference (%) full from out-of-Sample';
                current_col = current_col + 2;
            end
        end
        
        if (GARP_flags(7) == 1)
            col_headers{current_col} = 'SARP VIOLATIONS full index';
            current_col = current_col + 1;
            if (GARP_flags(3) == 1)
                col_headers{current_col} = 'SARP VIOLATIONS in-Sample';
                col_headers{current_col+1} = 'SARP VIOLATIONS in-Sample (%)';
                current_col = current_col + 2;
            end
            if (GARP_flags(4) == 1)
                col_headers{current_col} = 'SARP VIOLATIONS out-of-Sample index';
                col_headers{current_col+1} = 'SARP VIOLATIONS difference (%) full from out-of-Sample';
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
            if (VARIAN_flags(5) == 1)
                col_headers{current_col}   = 'VARIAN Mean in-Sample Component';
                col_headers{current_col+1} = 'VARIAN Mean Lower Bound in-Sample Component';
                col_headers{current_col+2} = 'VARIAN Mean Approximate in-Sample Component';
                col_headers{current_col+3} = 'VARIAN Mean Upper Bound in-Sample Component';
                col_headers{current_col+4} = 'VARIAN Mean in-Sample Difference';
                col_headers{current_col+5} = 'VARIAN Mean Approximate in-Sample Difference';
                current_col = current_col + 6;
            end
            if (VARIAN_flags(6) == 1)
                col_headers{current_col}   = 'VARIAN AVGSSQ in-Sample Component';
                col_headers{current_col+1} = 'VARIAN AVGSSQ Lower Bound in-Sample Component';
                col_headers{current_col+2} = 'VARIAN AVGSSQ Approximate in-Sample Component';
                col_headers{current_col+3} = 'VARIAN AVGSSQ Upper Bound in-Sample Component';
                col_headers{current_col+4} = 'VARIAN AVGSSQ in-Sample Difference';
                col_headers{current_col+5} = 'VARIAN AVGSSQ Approximate in-Sample Difference';
                current_col = current_col + 6;
            end
            if (VARIAN_flags(7) == 1)
                col_headers{current_col}   = 'VARIAN Max in-Sample Component';
                col_headers{current_col+1} = 'VARIAN Max Lower Bound in-Sample Component';
                col_headers{current_col+2} = 'VARIAN Max Approximate in-Sample Component';
                col_headers{current_col+3} = 'VARIAN Max Upper Bound in-Sample Component';
                col_headers{current_col+4} = 'VARIAN Max in-Sample Difference';
                col_headers{current_col+5} = 'VARIAN Max Approximate in-Sample Difference';
                current_col = current_col + 6;
            end
        end
        if (VARIAN_flags(4) == 1)
            if (VARIAN_flags(5) == 1) % Mean Aggregator
                col_headers{current_col}    = 'VARIAN Mean full index';
                col_headers{current_col+1}  = 'VARIAN Mean Lower Bound full index';
                col_headers{current_col+2}  = 'VARIAN Mean Approximate full index';
                col_headers{current_col+3}  = 'VARIAN Mean Upper Bound full index';
                col_headers{current_col+4}  = 'VARIAN Mean out-of-Sample index';
                col_headers{current_col+5}  = 'VARIAN Mean Lower Bound out-of-Sample index';
                col_headers{current_col+6}  = 'VARIAN Mean Approximate out-of-Sample index';
                col_headers{current_col+7}  = 'VARIAN Mean Upper Bound out-of-Sample index';
                col_headers{current_col+8}  = 'VARIAN Mean difference full from out-of-Sample';
                col_headers{current_col+9}  = 'VARIAN Mean normalized difference full from out-of-Sample';
                col_headers{current_col+10} = 'VARIAN Mean Approximate difference full from out-of-Sample';
                col_headers{current_col+11} = 'VARIAN Mean Approximate normalized difference full from out-of-Sample';
                current_col = current_col + 12;
            end
            if (VARIAN_flags(6) == 1) % AVG(SSQ) Aggregator
                col_headers{current_col}    = 'VARIAN AVGSSQ full index';
                col_headers{current_col+1}  = 'VARIAN AVGSSQ Lower Bound full index';
                col_headers{current_col+2}  = 'VARIAN AVGSSQ Approximate full index';
                col_headers{current_col+3}  = 'VARIAN AVGSSQ Upper Bound full index';
                col_headers{current_col+4}  = 'VARIAN AVGSSQ out-of-Sample index';
                col_headers{current_col+5}  = 'VARIAN AVGSSQ Lower Bound out-of-Sample index';
                col_headers{current_col+6}  = 'VARIAN AVGSSQ Approximate out-of-Sample index';
                col_headers{current_col+7}  = 'VARIAN AVGSSQ Upper Bound out-of-Sample index';
                col_headers{current_col+8}  = 'VARIAN AVGSSQ difference full from out-of-Sample';
                col_headers{current_col+9}  = 'VARIAN AVGSSQ normalized difference full from out-of-Sample';
                col_headers{current_col+10} = 'VARIAN AVGSSQ Approximate difference full from out-of-Sample';
                col_headers{current_col+11} = 'VARIAN AVGSSQ Approximate normalized difference full from out-of-Sample';
                current_col = current_col + 12;
            end
            if (VARIAN_flags(7) == 1) % Maximum Aggregator 
                col_headers{current_col}    = 'VARIAN Max full index';
                col_headers{current_col+1}  = 'VARIAN Max Lower Bound full index';
                col_headers{current_col+2}  = 'VARIAN Max Approximate full index';
                col_headers{current_col+3}  = 'VARIAN Max Upper Bound full index';
                col_headers{current_col+4}  = 'VARIAN Max out-of-Sample index';
                col_headers{current_col+5}  = 'VARIAN Max Lower Bound out-of-Sample index';
                col_headers{current_col+6}  = 'VARIAN Max Approximate out-of-Sample index';
                col_headers{current_col+7}  = 'VARIAN Max Upper Bound out-of-Sample index';
                col_headers{current_col+8}  = 'VARIAN Max difference full from out-of-Sample';
                col_headers{current_col+9}  = 'VARIAN Max normalized difference full from out-of-Sample';
                col_headers{current_col+10} = 'VARIAN Max Approximate difference full from out-of-Sample';
                col_headers{current_col+11} = 'VARIAN Max Approximate normalized difference full from out-of-Sample';
                current_col = current_col + 12;
            end
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
            col_headers{current_col+1} = 'HM residual index';
            col_headers{current_col+2} = 'HM difference full from residual';
            col_headers{current_col+3} = 'HM normalized difference full from residual';
            current_col = current_col + 4;
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