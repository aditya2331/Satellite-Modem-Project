function syn_a  = syndrome_finder(rx_word,alpha,n,t,p,pp)
    % rx_word : Received code word
    % alpha   : primitive element of choice
    % n       : length of the codeword
    % t       : Error Correcting capability of the RS Coding Scheme
    syn_a = gf(zeros(1,length(2*t)),p,pp);
    for i = 1:2*t
        al  = (alpha^i).^(0:1:n-1);
        res = al.*rx_word;
        sum_a = gf(0,p,pp);
        for j = 1:length(res)
            sum_a = sum_a + res(1,j);
        end
        syn_a(1,i) = sum_a;
    end
    
end