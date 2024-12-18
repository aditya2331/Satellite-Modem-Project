function corrected_message_bits = RSdecoder(CC_decoded_bits, RS_n, RS_k, GF_LUT)
    % CC_decoded_bits     :   Output bits of Soft Decision Viterbi Decoder
    % RS_n                :     Number of bits per codeword block
    % RS_k                :     Number of bits per message block
    % GF_LUT              :   Galois Field Lookup Table (Not being used in Current Implementation)
    
    % RS Code Parameters
    RS_q = log2(RS_n+1);        % GF(2^q)-field
    RS_t = (RS_n - RS_k)/2;     % 

    % Computing alpha
    % Finding all possible primitive elements
    ff = 0:1:2^RS_q -1;  % All elements in GF(2^q)
    ff_star = ff(1,2:end);       % All elements in GF(2^q) except 0
    factor_n = find_factors(RS_n);  
    ord_element = zeros(1,length(ff_star));
    
    for i = 1:length(ff_star)
        for j = 1:2^RS_q-1
            % if ff_star(1,i)^j == 1
            if gf_ele_exponent(ff_star(1,i),j,RS_q) == 1
                ord_element(1,i) = j;
                break;
            end
        end
    end
    
    % Putting All primitive element into a list
    % prim_ele = gf(zeros(1,0),q,primitive_poly);
    prim_ele = zeros(1,0); % Empty Matrix
    for i = 1:length(ord_element)
        if ord_element(1,i)==RS_n
            prim_ele(1,end+1) = ff_star(1,i);
        end
    end
    
    % Power Table;
    % alpha - Chosen primitive Element
    alpha = prim_ele(1,1);

    % Reshaping bit sequence into the original messsage matrix dimension of k rows and appropriate columns
    binary_matrix = reshape(CC_decoded_bits, RS_q, []).';            % Matrix of q rows and appropriate columns, each column represents 1 q bit chunk
    rec_codeword_stream = bin2dec(char(binary_matrix + '0'));        % q-bit chunk is converted to a decimal number in GF(2^q)
    rec_codeword_matrix = reshape(rec_codeword_stream, RS_n, []);    % Matrix of n columns and appropriate columns, each row represents 1 n block codeword
    rec_codeword_matrix = transpose(rec_codeword_matrix);

    corrected_codewords = zeros(size(rec_codeword_matrix,1), size(rec_codeword_matrix,2));  % To store final answer
    
    for i=1:size(rec_codeword_matrix,1)
        corrected_codewords(i,:) = rsDecode(rec_codeword_matrix(i,:),alpha,RS_n,RS_t,RS_q);
    end

    corrected_codewords = corrected_codewords';                                                % Corrected codewords     
    corrected_codeword_bits = reshape(dec2bin(corrected_codewords(:),RS_q).' - '0', [], 1);    % Corrected codewords in bits

    % Decode the message from the systematic codeword
    corrected_message = corrected_codewords(1:RS_k,:);                                         % Corrected messages
    corrected_message_bits = reshape(dec2bin(corrected_message(:), RS_q).' - '0', [], 1);      % Converting messages in bits
end

