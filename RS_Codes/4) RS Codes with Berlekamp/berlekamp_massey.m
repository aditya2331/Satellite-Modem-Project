function [Lambda, L] = berlekamp_massey(syndromes, t, p, pp)
    n = length(syndromes);
    Lambda = gf([1 zeros(1, n)], p, pp);  % Initialize Lambda(x) = 1
    B = gf([1 zeros(1, n)], p, pp);  % B(x) = 1
    L = 0;  % Degree of the error locator polynomial
    m = 1;  % Step counter
    b = gf(1, p, pp);  % Initialize the discrepancy scalar
    for r = 1:n
        delta = syndromes(r);
        for i = 1:L
            delta = delta + Lambda(i+1) * syndromes(r-i);
        end
        if delta == 0
            m = m + 1;
        else
            T = Lambda;
            Lambda = Lambda - delta / b * [zeros(1, m) B(1:end-m)];
            if 2 * L <= r - 1
                B = T;
                L = r - L;
                b = delta;
                m = 1;
            else
                m = m + 1;
            end
        end
    end
    Lambda = Lambda(1:L+1);
end