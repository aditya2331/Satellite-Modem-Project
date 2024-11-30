% Syndrome Finder
function syn_a  = syndrome_finder(rx_word, alpha, n, t, p)
    % rx_word : Received code word
    % alpha   : primitive element of choice
    % n       : length of the codeword
    % t       : Error Correcting capability of the RS Coding Scheme
    
    persistent add_table_cache p_cache;
    
    % Load the addition table only if it hasn't been loaded or if 'p' has changed
    if isempty(add_table_cache) || p_cache ~= p
        data = load("Tables.mat");
        add_table_cache = data.(sprintf('AT%d', p));
        p_cache = p;  % Cache the current 'p' value
    end

    add_table = add_table_cache;

    % syn_a = gf(zeros(1,length(2*t)),p,pp);
    syn_a = zeros(1, length(2*t));
    for i = 1:2*t
        % al  = (alpha^i).^(0:1:n-1);
        temp = gf_ele_exponent(alpha,i,p);
        al = zeros(1,n);
        for j = 1:n
            al(j) = gf_ele_exponent(temp,j-1,p);
        end

        % res = al.*rx_word;
        res = gf_ele_mat_mul(al, rx_word, p);

        % sum_a = gf(0,p,pp);
        sum_a = 0;

        for j = 1:length(res)
            % sum_a = sum_a + res(1,j);
            sum_a = add_table(sum_a+1, res(1,j)+1);
        end

        syn_a(1,i) = sum_a;
    end
end