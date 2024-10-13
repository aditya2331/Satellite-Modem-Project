% Takes Demodulated bits and performs Syndrome Decoding and returns
% corrected received message as bits

function [corrected_rec_codeword_bits, corrected_rec_message_bits] = rs_decoder(rec_bits,n,k,q,t,alpha,power_to_ele_dict,ele_to_power_dict,add_table)
    % Reshaping bit sequence into the original messsage matrix dimension of
    % k rows and appropriate columns
    bin_matrix = reshape(rec_bits, q, []).';
    rec_message_stream = bin2dec(char(bin_matrix + '0'));
    rec_message = reshape(rec_message_stream, n, []);
    % num_cols = size(rec_message,2);
    
    rec_message_transpose = rec_message;
    rec_message = rec_message';
    
    corr_rec_message = zeros(size(rec_message,1), size(rec_message,2));

    for i=1:size(rec_message,1)
        corr_rec_message(i,:) = rsDecode(rec_message(i,:),alpha,n,t,q,power_to_ele_dict,ele_to_power_dict,add_table);
    end

    corrected_rec_codeword = corr_rec_message';
    corrected_rec_codeword_bits = reshape(dec2bin(corrected_rec_codeword(:),q).' - '0', [], 1);
    
    code_dec = corrected_rec_codeword;

    % Decode the codewords
    corrected_rec_message = corrected_rec_codeword(1:k,:); 
    % corrected_rec_message = corrected_rec_codeword(n-k+1:n,:); 
    corrected_rec_message_bits = reshape(dec2bin(corrected_rec_message(:), q).' - '0', [], 1); % Converting final_rec_message to a binary sequence
   
    mes_out = corrected_rec_message;
end

