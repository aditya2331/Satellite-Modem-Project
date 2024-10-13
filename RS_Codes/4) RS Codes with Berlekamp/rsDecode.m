function [rxcodeWord,flag] = rsDecode(rx_word,alpha,n,t,p,pp,power_to_ele_dict,ele_to_power_dict)
    % Flag bit
    % 0 - Decoding Successful (<=t errors case)
    % 1 - Decoding Failed     (> t  errors case)
    rx_word1 = flip(rx_word); 
    syn_a = syndrome_finder(rx_word1,alpha,n,t,p,pp);
    [Lambda, L] = berlekamp_massey(syn_a, t, p, pp);

    num_list = 0:1:length(Lambda)-1;
    error_indx = zeros(1,0);
    for i = 0:n-1
        al = alpha^i;
        al = gf_inverse(al,p,pp,power_to_ele_dict,ele_to_power_dict);
        al_list = al.^num_list;
        res_x = al_list.*Lambda;
        res= 0;
        for j=1:length(res_x)
            res = res + res_x(1,j);
        end
        if res==0
            error_indx(1,end+1) = i+1;
        end
    end
    error_indx_1 = error_indx - 1;
    leng = length(error_indx_1);
    if leng~=0
        G = gf(zeros(leng,leng),p,pp);
        for i = 1:leng
            G(i,:) = (alpha^i).^error_indx_1;
        end
        G_prime = [G transpose(syn_a(1,1:leng))];
        M = gf_gaussian_eli(G_prime,p,pp,power_to_ele_dict,ele_to_power_dict);
        Y = transpose(M.x(:,end));
        mapped_code_word = rx_word1;
        for i=1:leng
            err_i = error_indx(1,i);
            err = Y(1,i);
            mapped_code_word(1,err_i) = mapped_code_word(1,err_i) + err;
        end
    else
        mapped_code_word = rx_word1;
    end
    % Checking whether correction is done
    flag = 1; % 1 means Error Correction is failed
    syn_a_corr = syndrome_finder(mapped_code_word,alpha,n,t,p,pp);
    dum_val = sum(syn_a_corr.x);
    if dum_val == 0
        flag = 0; % Correction is successful
    end
    rxcodeWord = flip(mapped_code_word);
end