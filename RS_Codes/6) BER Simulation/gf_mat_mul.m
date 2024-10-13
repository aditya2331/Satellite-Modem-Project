% This function can ONLY BE USED FOR 2D MATRIX MULTIPLICATION

function [C] = gf_mat_mul(A, B, q)
    A_rows = size(A,1);
    A_cols = size(A,2);

    B_rows = size(B,1);
    B_cols = size(B,2);
    
    C = zeros(A_rows, B_cols);  

    % Check if multiplication is valid
    if(A_cols ~= B_rows)
        disp('Invalid Matrix Multiplication -  Dimension Mismatch');
        return;
    end
    
    data = load("Tables.mat");
    add_table = data.(sprintf('AT%d',q));
    mul_table = data.(sprintf('MT%d',q));

    for i=1:A_rows
        for j=1:B_cols
            row = A(i,:);  % Extract the ith row of A 
            col = B(:,j);  % Extract the jth column of B 

            final_mul = 0;
            for ind=1:A_cols
                curr_mul = mul_table(row(ind)+1,col(ind)+1);
                final_mul = add_table(curr_mul+1, final_mul+1); % For Addition, we need to add 1 to the indices
            end
            C(i,j) = final_mul;
        end
    end
end

