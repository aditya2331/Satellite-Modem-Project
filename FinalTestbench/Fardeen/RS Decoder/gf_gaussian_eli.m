% Gaussian Elimination in Galois Field
function M = gf_gaussian_eli(G,p)
    persistent add_table_cache mul_table_cache p_cache;
    
    % Load tables only if they haven't been loaded or if 'p' has changed
    if isempty(add_table_cache) || isempty(mul_table_cache) || p_cache ~= p
        data = load("Tables.mat");
        add_table_cache = data.(sprintf('AT%d', p));
        mul_table_cache = data.(sprintf('MT%d', p));
        p_cache = p;  % Cache the current 'p' value
    end

    add_table = add_table_cache;
    mul_table = mul_table_cache;
    
    persistent power_to_ele_dict_cache ele_to_power_dict_cache;

    % Load dictionaries if not already loaded
    if isempty(power_to_ele_dict_cache) || isempty(ele_to_power_dict_cache)
        dicts = load("Dicts.mat");  % Load dictionaries once
        power_to_ele_dict_cache = dicts.power_to_ele_dict;
        ele_to_power_dict_cache = dicts.ele_to_power_dict;
    end

    power_to_ele_dict = power_to_ele_dict_cache;
    ele_to_power_dict = ele_to_power_dict_cache;

    M = G;
    [m,n] = size(M);
    row = 1;
    pivot = 1;

    while(row<=m && pivot<=m)
        % if M(row,pivot) == gf(0,p,pp)
        if M(row, pivot) == 0
            m1 = transpose(M(row+1:end,pivot));
            nz_indx = find(m1~=0);
            if isempty(nz_indx)
                pivot = pivot+1;
                continue;
            else
                % swapping row to bring non-zero pivot to the current row location
                nz_indx = nz_indx(1,1);
                x = M(row,:);
                M(row,:) = M(nz_indx,:);
                M(nz_indx,:) = x;
            end
        end
        pv = M(row,pivot);
        pv_inv = gf_inverse(pv,p);
        % M(row,:) = pv_inv*M(row,:);
        for i = 1:n
            M(row,i) =  mul_table(pv_inv+1,M(row,i)+1);
        end

        for j = 1:m
            if j ~= row
                if M(j,pivot) ~= 0
                    alp = M(j,pivot);
                    % M(j,:) = M(j,:) + alp*M(row,:);
                    for k = 1:n
                        M(j,k) = add_table(M(j,k)+1, mul_table(alp+1, M(row,k)+1)+1);
                    end
                end
            end
        end
        row = row + 1;
        pivot = pivot + 1;
    end
end