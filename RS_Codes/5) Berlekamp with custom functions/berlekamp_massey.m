% Function that computes Error Locator Polynomial
function [Lambda, L] = berlekamp_massey(syndromes, ~, q, power_to_ele_dict,ele_to_power_dict)
    n = length(syndromes);
    % Lambda = gf([1 zeros(1, n)], p, pp);  % Initialize Lambda(x) = 1
    Lambda = [1 zeros(1,n)];
    % B = gf([1 zeros(1, n)], p, pp);  % B(x) = 1
    B = [1 zeros(1,n)];
    L = 0;  % Degree of the error locator polynomial
    m = 1;  % Step counter
    % b = gf(1, p, pp);  % Initialize the discrepancy scalar
    b = 1;
    data = load("Tables.mat");
    add_table = data.(sprintf('AT%d',q));
    mul_table = data.(sprintf('MT%d',q));
    for r = 1:n
        delta = syndromes(r);
        for i = 1:L
            % delta = delta + Lambda(i+1) * syndromes(r-i);
            delta = add_table(delta+1, mul_table(Lambda(i+1)+1, syndromes(r-i)+1)+1);
        end
        
        if delta == 0
            m = m + 1;
        else
            T = Lambda;
            % Lambda = Lambda - delta / b * [zeros(1, m) B(1:end-m)];
            temp_vec = [zeros(1,m) B(1:end-m)];
            temp_ele = mul_table(delta+1, gf_inverse(b,q,power_to_ele_dict,ele_to_power_dict)+1);
            for i=1:length(temp_vec)
                temp_vec(i) = mul_table(temp_ele + 1, temp_vec(i)+1);
            end
            % Lambda = Lambda - temp_vec; - Adddition and Subtraction are the same in GF(2)
            Lambda = gf_mat_add(Lambda,temp_vec,q);
            
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