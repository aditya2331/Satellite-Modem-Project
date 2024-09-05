function [C] = galois_mat_mul(A,B,primitive_poly)
    % Because we are working with 3 bit codewords, ie, 0 to 7
    A = mat_dec_to_bin(A,3); 
    B = mat_dec_to_bin(B,3);

    A_rows = size(A,1);
    A_cols = size(A,2);
    A_height = size(A,3); % deg(vector in A) + 1

    B_rows = size(B,1);
    B_cols = size(B,2);
    B_height = size(B,3); % deg(vector in B) + 1
    
    % length of vector in C = deg(A) + deg (B) + 1 = len(A) + len(B) - 1
    C = zeros(A_rows, B_cols, 3);  % 3 because we are working with 3 bit codewords

    % Check if multiplication is valid
    if(A_cols ~= B_rows)
        disp('Invalid Matrix Multiplication -  Dimension Mismatch');
        return;
    end
    if(A_height ~= B_height)
        disp('Invalid Matrix Multiplication - Matrix Length is different')
    end

    for i=1:A_rows
        for j=1:B_cols
            row = squeeze(A(i,:,:));  % Extract the ith row of A (1x3x3 -> 3x3)
            col = squeeze(B(:,j,:));  % Extract the jth column of B (3x1x3 -> 3x3)

            final_mul = zeros(1,3); % Because we are working with 3 bit codewords
            for ind=1:A_cols
                mul_a = row(ind,:); % Row vectors are stored in binary format rowwise, Row 1 contains binary rep of element 1, row 2 contains binary rep of element 2 and so on
                mul_b = col(ind,:); % Column vectors are stored in binary format rowise
                
                curr_mul = poly_mul(mul_a,mul_b,primitive_poly);
                final_mul = poly_sum(curr_mul, final_mul);
            end
            final_mul = reshape(final_mul, [1,1,3]); % To match the dimension of C(i,j,:)
            C(i,j,:) = final_mul;
        end
    end
    
    C = mat_bin_to_dec(C);
end

