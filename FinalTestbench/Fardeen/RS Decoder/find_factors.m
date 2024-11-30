function factors = find_factors(n)
    factors = [];  % Initialize an empty array to store factors 
    for i = 1:sqrt(n)
        if mod(n, i) == 0
            factors = [factors, i];  % Add the factor
            if i ~= n/i
                factors = [factors, n/i];  % Add the corresponding pair factor
            end
        end
    end
    factors = sort(factors); 
end

