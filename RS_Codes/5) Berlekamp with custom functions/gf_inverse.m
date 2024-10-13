% Calculates the multiplicative inverse of a Galois Field Element

function inv = gf_inverse(gf_ele,m,power_to_ele_dict,ele_to_power_dict)
    % gf_ele1 = gf_ele.x;
    n = ele_to_power_dict(gf_ele);
    n_inv = 2^m - 1 - n;
    inv = power_to_ele_dict(n_inv);
    % inv = gf(inv,m,pp);
end