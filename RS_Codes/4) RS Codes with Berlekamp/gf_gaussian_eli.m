function M = gf_gaussian_eli(G,p,pp,power_to_ele_dict,ele_to_power_dict)
    M = G;
    [m,n] = size(M);
    row = 1;
    pivot = 1;
    while(row<=m && pivot<=m)
        if M(row,pivot) == gf(0,p,pp)
            m1 = transpose(M(row+1:end,pivot));
            nz_indx = find(m1~=0);
            if isempty(nz_indx)
                pivot = pivot +1;
                continue;
            else
                % swapping row to bring non-zero pivot to the current row
                % location
                nz_indx = nz_indx(1,1);
                x = M(row,:);
                M(row,:) = M(nz_indx,:);
                M(nz_indx,:) = x;
            end
            
        end
        pv = M(row,pivot);
        pv_inv = gf_inverse(pv,p,pp,power_to_ele_dict,ele_to_power_dict);
        M(row,:) = pv_inv*M(row,:);
        for j = 1:m
            if j ~= row
                if M(j,pivot)~=gf(0,p,pp)
                    alp = M(j,pivot);
                    M(j,:) = M(j,:) + alp*M(row,:);
                end
            end
        end
        row = row + 1;
        pivot = pivot + 1;
    end
end