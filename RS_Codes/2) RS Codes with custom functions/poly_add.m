% Function takes 2 polynomials and performs bitwise XOR
% Assumes that input polynomials are in Big Endian Format

function [res] = poly_add(poly_a, poly_b)
    deg_a = length(poly_a) - 1;
    deg_b = length(poly_b) - 1;

    if deg_a >= deg_b
        deg_min = deg_b;
        res = poly_a;
    else
        deg_min = deg_a;
        res = poly_b;
    end
    
    % Converting to Little Endian Format
    poly_a = flip(poly_a);
    poly_b = flip(poly_b);
    res = flip(res);
    for i=0:deg_min
        res(i+1) = bitxor(poly_a(i+1), poly_b(i+1));
    end
    res = flip(res); % Converting back to Big Endian Format
end