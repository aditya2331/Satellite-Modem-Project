% Takes Input Message Stream in GF(2^q) and returns codewords as a bit sequence

function [message_bits, c_bits] = rs_encoder(message, k, q, G)
    message_bits = reshape((dec2bin(message,q)-'0')',[],1);           % Input message in GF(2^m) as bits
    message = reshape(message, k, []);            % Reshaping message into a matrix with k rows
    c = gf_mat_mul(G, message, q);                % Generating codewords
    c_bits = reshape((dec2bin(c,q) - '0')', [], 1); % Converting codewords to a binary sequence
end