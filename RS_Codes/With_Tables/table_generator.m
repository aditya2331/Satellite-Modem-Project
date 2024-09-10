function [AT,MT] = table_generator(q, prim)
    AT = zeros(2^q,2^q);
    for i=1:2^q
        for j=1:2^q
            bin_i = arrayfun(@(x) str2double(x), dec2bin(i-1,q));
            bin_j = arrayfun(@(x) str2double(x), dec2bin(j-1,q));
            bin_add_res = flip(poly_add(bin_i,bin_j)); % bin_add_res is in Little Endian
            for k=1:length(bin_add_res)
                AT(i,j) = AT(i,j)+bin_add_res(k)*(2^(k-1));
            end
        end
    end
    
    MT = zeros((2^q),(2^q));
    for i=1:2^q
        for j=1:2^q
            bin_i = arrayfun(@(x) str2double(x), dec2bin(i-1,q));
            bin_j = arrayfun(@(x) str2double(x), dec2bin(j-1,q));
            bin_mul_res = flip(poly_mul(bin_i,bin_j,prim,q)); % bin_mul_res is in Little Endian
            for k=1:length(bin_mul_res)
                MT(i,j) = MT(i,j)+bin_mul_res(k)*(2^(k-1));
            end
        end
    end
end

