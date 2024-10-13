function [mul_res] = poly_mul(poly_a, poly_b, poly_prim, q)
    % Base Case - One of the vectors is zero
    if(all(poly_a == 0) || all(poly_b == 0))
        mul_res = zeros(1,q);  % Return 000..0 (q times) because we are working with q bit codewords
        return; 
    end

    % It is assumed that input polynomials are in Big Endian format
    % ie, [1 1 0 1] = x^3 + x^2 + 1

    % For Multiplication, it is easier to perform computations in little
    % endian format, ie, [1 1 0] = 1 + x
    poly_a_little = flip(poly_a);
    poly_b_little = flip(poly_b);

    % POLYNOMIALS MUST NOT HAVE REDUNDANT HIGHER ORDER COEFFECIENTS -
    % Otherwise code doesn't work
    % Removing trailing zeros
    for i=length(poly_a_little):-1:1
        if poly_a_little(end) == 0
            poly_a_little(end) = [];
        else
            break;
        end
    end
    
    for i=length(poly_b_little):-1:1
        if poly_b_little(end) == 0
            poly_b_little(end) = [];
        else
            break;
        end
    end

    % Degree calculation is done on the basis that last element in the
    % polynomial is a 1 after redundant terms are added

    deg_a = length(poly_a_little) - 1;
    deg_b = length(poly_b_little) - 1;
    
    deg_res = 2*q;                 % Length of Multiplication is 2 * q
    res = zeros(1, deg_res + 1);   % Array's length is 1 more than it's degree
    curr_mul = zeros(1,deg_res+1); % Stores current multiplication result
    
    % Set the polynomial with less number of non zero elements as the multiplier so that there is less
    % number of XOR operations

    % No: of XOR operations = number of non zero elements of smaller multiplier
    if nnz(poly_a_little) <= nnz(poly_b_little)
        multiplier = poly_a_little;
        multiplicand = poly_b_little;
    else
        multiplier = poly_b_little;
        multiplicand = poly_a_little;
    end

    for i = 1:length(multiplier)
        curr_mul(:) = 0;
        if multiplier(i) == 1
            for j=1:length(multiplicand)
                if multiplicand(j) == 1
                    curr_mul(i+j-1) = 1;
                end
            end
        else
            curr_mul(:) = 0;
        end
        res = poly_add(res,curr_mul); % res is in Little Endian format
    end

    % Division with primitive polynomial

    % If degree of multiplication result is more than degree of primitive
    % polynomial, no need to divide with it - CHECK IF THIS IS CORRECT

    % Division is done with polynomials in Big Endian format
    if(length(res) < length(poly_prim))
        remainder = flip(res);
    else
        res_big = flip(res); % In Big Endian format, ie, [1 1 0] = x^2 + x 
        num_loops = length(res_big)-length(poly_prim)+1;

        for i = 1:num_loops
        % XOR only if leading factor in res_bin = 1 (Just like normal divison)
            if res_big(i) == 1
                for j = i:length(poly_prim)+i-1
                    res_big(j) = poly_add(res_big(j), poly_prim(j-i+1));
                end
            end
        end
        remainder = res_big;
    end

    % Remainder is already in Big Endian format
    if length(remainder) >= q % Only keep last q elements as final answer should be less than 2^q
        mul_res = remainder(end-q:end);
    else 
        mul_res = [zeros(1,q-length(remainder)), remainder]; % Prepending zeros to make the final answer into q bits
    end
end
