% Berlekamp Algorithm
function [rxcodeWord,flag] = rsDecode(rx_word,alpha,n,t,p)
    % Flag bit
    % 0 - Decoding Successful (<=t errors case)
    % 1 - Decoding Failed     (> t  errors case)
    
    persistent add_table_cache p_cache;
    
    % Load the addition table only if it hasn't been loaded or if 'p' has changed
    if isempty(add_table_cache) || p_cache ~= p
        data = load("Tables.mat");
        add_table_cache = data.(sprintf('AT%d', p));
        p_cache = p;  % Cache the current 'p' value
    end

    add_table = add_table_cache;

    persistent power_to_ele_dict_cache ele_to_power_dict_cache;
    % Load dictionaries if not already loaded
    if isempty(power_to_ele_dict_cache) || isempty(ele_to_power_dict_cache)
        dicts = load("Dicts.mat");  % Load dictionaries once
        power_to_ele_dict_cache = dicts.power_to_ele_dict;
        ele_to_power_dict_cache = dicts.ele_to_power_dict;
    end

    power_to_ele_dict = power_to_ele_dict_cache;
    ele_to_power_dict = ele_to_power_dict_cache;

    rx_word = flip(rx_word);
    syn_a = syndrome_finder(rx_word, alpha, n, t, p);
    [Lambda, ~] = berlekamp_massey(syn_a, t, p);

    num_list = 0:1:length(Lambda)-1;
    error_indx = zeros(1,0);

    for i = 0:n-1
        % al = alpha^i;
        al = gf_ele_exponent(alpha,i,p);
        al = gf_inverse(al,p);
        % al_list = al.^num_list;
        al_list = zeros(1,length(num_list));
        for j=1:length(num_list)
            al_list(j) = gf_ele_exponent(al,num_list(j),p);
        end
        % res_x = al_list.*Lambda;
        res_x = gf_ele_mat_mul(al_list, Lambda, p);
        res= 0;

        for j=1:length(res_x)
            % res = res + res_x(1,j);
            res = add_table(res+1, res_x(1,j) + 1);
        end

        if res==0
            error_indx(1,end+1) = i+1;
        end
    end

    error_indx_1 = error_indx - 1;
    leng = length(error_indx_1);

    if leng~=0
        % G = gf(zeros(leng,leng),p,pp);
        G = zeros(leng, leng);

        for i = 1:leng
            % G(i,:) = (alpha^i).^error_indx_1;
            temp = gf_ele_exponent(alpha,i,p);
            for j = 1:length(error_indx_1)
                G(i,j) = gf_ele_exponent(temp,error_indx_1(j),p);
            end
        end

        G_prime = [G transpose(syn_a(1,1:leng))];
        M = gf_gaussian_eli(G_prime,p);
        Y = transpose(M(:,end));
        mapped_code_word = rx_word;
        for i=1:leng
            err_i = error_indx(1,i);
            err = Y(1,i);
            % mapped_code_word(1,err_i) = mapped_code_word(1,err_i) + err;
            mapped_code_word(1,err_i) = add_table(mapped_code_word(1,err_i)+1, err+1);
        end
    else
        mapped_code_word = rx_word;
    end

    % Checking whether correction is done
    flag = 1; % 1 means Error Correction is failed
    syn_a_corr = syndrome_finder(mapped_code_word,alpha,n,t,p);
    dum_val = sum(syn_a_corr);
    if dum_val == 0
        flag = 0; % Correction is successful
    end
    rxcodeWord = flip(mapped_code_word);
end