% This function can ONLY BE USED FOR 2D MATRIX ADDITION

function [C] = gf_mat_add(A, B, q)
    persistent add_table_cache q_cache;
    
    % Load the addition table only if it hasn't been loaded or if 'q' has changed
    if isempty(add_table_cache) || q_cache ~= q
        data = load("Tables.mat");
        add_table_cache = data.(sprintf('AT%d', q));
        q_cache = q;  % Cache the current 'q' value
    end

    add_table = add_table_cache;

    % Check if addition is invalid
    if((size(A,1) ~= size(B,1)) || (size(A,2) ~= size(B,2)))
        disp("Dimensions don't match");
        return;
    end

    C = zeros(size(A,1),size(A,2));

    for i=1:size(A,1)
        for j=1:size(A,2)
            C(i,j) = add_table(A(i,j)+1,B(i,j)+1); % Value in Table is at Index + 1
        end
    end
end

