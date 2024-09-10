function [C] = galois_mat_add(A, B, q)
    % Check if addition is invalid
    if((size(A,1) ~= size(B,1)) || (size(A,2) ~= size(B,2)))
        disp("Dimensions don't match");
        return;
    end
    
    % Messages can take values 0 to 2^q-1
    A = mat_dec_to_bin(A,q); 
    B = mat_dec_to_bin(B,q);

    for i=1:size(A,1)
        for j=1:size(A,2)
            row = squeeze(A(i,j,:));  % Extract the ith row of A (1xqxq -> 1xq)
            col = squeeze(B(i,j,:));  % Extract the jth column of B (qx1xq -> 1xq)
            final_add = poly_add(row,col);
            C(i,j,:) = final_add;
        end
    end
    C = mat_bin_to_dec(C);
end

