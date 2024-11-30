function c_bits = RSencoder(scrambled_bits, RS_n, RS_k, GF_LUT)
    % RS_n         :     Number of bits per codeword block
    % RS_k         :     Number of bits per message block
    % GF_LUT       :     Lookup Table for Galois Field Addition and Multiplication

    % Currently, I am not using GF_LUT as in this implementation, GF_LUT is
    % being used from a table inside the function itself based on RS_q

    % The code assumes that the number of input bits is a multiple of
    % (RS_q*RS_k)

    % Because q-bit chunks are converted to decimal numbers 
    % And k-bit blocks of GF(2^q) are then converted to n-bit blocks by
    % multiplying it with the Generator Matrix

    % RS Code parameters
    RS_q = log2(RS_n+1);        % GF(2^q)-field

    % Loading the Generator Matrix
    data_g = load("galois_matrices.mat");                   % Contains Generator Matrix for different (n,k) 
    RS_G = data_g.(sprintf('G%d',RS_q));

    binary_matrix = reshape(scrambled_bits, RS_q, []).';    % Reshaping the input bit sequence to a matrix of q bit columns (each column has q bits)
    dec_num_stream = bin2dec(char(binary_matrix + '0'));
    dec_num = reshape(dec_num_stream, RS_k, []);            % Each column represents 1 k bit block in GF(2^q)

    c = gf_mat_mul(RS_G, dec_num, RS_q);                    % Generating codewords

    c_bits = reshape((dec2bin(c,RS_q) - '0')', [], 1);      % Converting codewords to a column vector of bits
end

