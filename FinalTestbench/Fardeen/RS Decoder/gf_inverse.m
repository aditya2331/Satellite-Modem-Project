% Calculates the multiplicative inverse of a Galois Field Element

function inv = gf_inverse(gf_ele,m)
    persistent power_to_ele_dict_cache ele_to_power_dict_cache;
    % Load dictionaries if not already loaded
    if isempty(power_to_ele_dict_cache) || isempty(ele_to_power_dict_cache)
        dicts = load("Dicts.mat");  % Load dictionaries once
        power_to_ele_dict_cache = dicts.power_to_ele_dict;
        ele_to_power_dict_cache = dicts.ele_to_power_dict;
    end

    power_to_ele_dict = power_to_ele_dict_cache;
    ele_to_power_dict = ele_to_power_dict_cache;

    % gf_ele1 = gf_ele.x;
    n = ele_to_power_dict(gf_ele);
    n_inv = 2^m - 1 - n;
    inv = power_to_ele_dict(n_inv);
    % inv = gf(inv,m,pp);
end