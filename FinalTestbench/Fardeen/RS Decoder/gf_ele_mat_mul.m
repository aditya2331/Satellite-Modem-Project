% This function can ONLY BE USED FOR 2D ELEMENTWISE MATRIX MULTIPLICATION

function [C] = gf_ele_mat_mul(A,B,q)
    persistent mul_table_cache q_cache;
    
    % Load the multiplication table only if it hasn't been loaded or if 'q' has changed
    if isempty(mul_table_cache) || q_cache ~= q
        data = load("Tables.mat");
        mul_table_cache = data.(sprintf('MT%d', q));
        q_cache = q;  % Cache the current 'q' value
    end
    
    mul_table = mul_table_cache;

    A_rows = size(A,1);
    A_cols = size(A,2);

    B_rows = size(B,1);
    B_cols = size(B,2);
    
    C = zeros(A_rows, B_cols);  

    % Check if multiplication is valid
    if(A_rows ~= B_rows)
        disp('Invalid Matrix Multiplication -  Row Dimension Mismatch');
        return;
    end

    if(A_cols ~= B_cols)
        disp('Invalid Matrix Multiplication -  Column Dimension Mismatch');
        return;
    end

    for i=1:A_rows
        for j=1:B_cols
            C(i,j) = mul_table(A(i,j) + 1, B(i,j) + 1);
        end
    end
end

