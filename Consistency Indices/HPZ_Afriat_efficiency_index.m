function Index = HPZ_Afriat_efficiency_index (expenditure, index_threshold)

% this function calculates the Afriat inconsistency index,
% by using a binary search approach ("Lion in the Desert")



num_of_iterations = 30;

AFRIAT_UPPER = 1; % The maximum value for Afriat index.

AFRIAT_LOWER = 0; % The minimum value for Afriat index.

AFRIAT = 1/2;

[rows,cols] = size(expenditure); % rows=cols=number of observations 

for i=1:num_of_iterations
    
    % diag(diag(expenditure)) is a matrix with the same diagonal 
    % as expenditure, and zero at the cell in the i'th row and 
    % the j'th row column (i~=j).
    % af_exp has at the cell in the i'th row and the 
    % j'th column, where i doesn't equal j, expenditure(i,j), 
    % otherwise (i=j),AFRIAT*expenditure(i,i).        

    af_exp = expenditure - (diag(diag(expenditure))*(1-AFRIAT));
   
    %The matrix REF has at the cell in the i'th row and the j'th
    %column, in case i doesn't equal j, the difference between 
    %the value of the bundle that was chosen in observation i
    %multiplied by AFRIAT and the bundle that was chosen in observation j given the 
    %prices of observation i, and 0 on the diagonal.

    REF = diag(af_exp)*ones(rows,1)' - af_exp;

    %Denote v the following m dimension vector: 
    % v=(AFRIAT,AFRIAT,…,AFRIAT).

    %The matrix DRP has at the cell in the i'th row and the j'th
    %column, 1 if and only if the bundle that was chosen in 
    %observation i is directly v-revealed preferred to the bundle that was chosen 
    %in observation j (the corresponding value in REF is greater than -1).
    
    % DRP is the round up of the matrix K that satisfies: 
    % K*(max(max(abs(REF+index_threshold)))+1)= REF+index_threshold.
    % So, DRP(i,j)=1 iff 0<K(i,j)<=1.

    % Note that:
    % a. max(max(abs(REF+index_threshold))) is a positive number, 
    %    because   REF(i,i)=0 and  index_threshold>0. Therefore, 
    %    max(max(abs(REF+index_threshold)))+1 >1.  
    % b. For every i,j:
    %    max(max(abs(REF+index_threshold)))+1 >   (REF+index_threshold)(i,j) 
    %    conclusion:(REF+index_threshold)(i,j)>0 iff 0<K(i,j)<1.  

    DRP = ceil((REF+index_threshold)/(max(max(abs(REF+index_threshold)))+1));

    %The matrix SDRP has at the cell in the i'th row and the j'th
    %column, 1 if and only if the bundle that was chosen in 
    %observation i is strictly directly v-revealed preferred to 
    %the bundle that was chosen 
    %in observation j (the corresponding value in REF is greater than -1).

    SDRP = ceil((REF-index_threshold)/(max(max(abs(REF-index_threshold)))+1));

    % statement needed for the graph theory external package

    set_matlab_bgl_default(struct('full2sparse',1));

    %The matrix NS_RP has at the cell in the i'th row and the j'th
    %column, Inf if and only if the bundle that was chosen in 
    %observation i is not v-revealed preferred to the bundle that was chosen 
    %in observation j. Otherwise it includes a positive integer.

    NS_RP = all_shortest_paths(DRP);

    %The matrix RP has at the cell in the i'th row and the j'th
    %column, 1 if and only if the bundle that was chosen in 
    %observation i is v-revealed preferred to the bundle that was chosen 
    %in observation j. Otherwise, it equals 0. 
    
    RP = zeros(rows);

    for j=1:rows
        for k=1:cols               % going through all NS_RP’s cells.
            if ~isinf(NS_RP(j,k))
                
                %if the length of the path from j to k is finite, meaning there   
                %is a path from j to k (or in other words bundle j is v-revealed 
                %preferred to bundle k), then RP(j,k)=1, otherwise RP(j,k) stays 0.
                
                RP(j,k)=1;
            end
        end
    end

    %To test for GARP we will do the following: for every pair of
    %choices x and y if xR(AFRIAT)y then not yP0(AFRIAT)x. We will take RP and the transpose of
    %SDRP and multiply element by element. Every cell with the value 1 corresponds to 
    %xRy and yP0x, otherwise the value is 0. The final matrix is the zero matrix if and only if 
    %GARP_AFRIAT is satisfied. 

    GARP = RP.*(SDRP');
    
    %sums the values of all matrix GARP's cells. Counting the number of violations (notice 
    %that a pair of bundles x^i,x^j might be counted as two violations; one for (i,j), and 
    %the other for (j,i)).

    GARP_ERRORS = sum(sum(GARP));

    %If the GARP matrix is only zeros then the data satisfies GARP
    
    if GARP_ERRORS==0
        
        %if the data satisfies GARP_v, then it narrows the range to 
        %[AFRIAT, AFRIAT_UPPER], so in the next iteration, in case it's not the 
        %last one, the algorithm will try to look for a smaller parallel shifting 
        %(which means a bigger AFRIAT where the data satisfy GARP_v).    

        AFRIAT_LOWER = AFRIAT;
        
    else
        
        %if the data doesn't satisfy GARP_v, then it narrows the range to 
        %[AFRIAT_LOWER, AFRIAT] because if the data doesn’t satisfy GARP_v, 
        %it won't satisfy GARP_u for every u>v

        AFRIAT_UPPER = AFRIAT;
        
    end 

    AFRIAT = (1/2)*(AFRIAT_LOWER + AFRIAT_UPPER);        
end

%AFRIAT_LOWER is the lowest rate found where there are no  
%violations of GARP under a parallel shifting in that rate. 

AFRIAT = AFRIAT_LOWER;      

Index = 1 - AFRIAT;

end

