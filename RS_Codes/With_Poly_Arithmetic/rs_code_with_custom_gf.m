clear;
%% Generator Matrix, Codewords, Parity Check Matrix
% RS Code parameters
k = 3; % Bits/message
n = 7; % Bits/codeword
q = 3; % n = 2^q-1
t = 2; % Error correction capability
primitive_poly = [1 0 1 1]; % In Big Endian Format

% Random message
length = 10*k;
bitstream = randi([0 n], 1, length);

% Sampling k bits and creating column vectors
num_cols = length/k;
message = reshape(bitstream, k, []);
inp_bits = reshape(dec2bin(message,k)-'0', [], 1);

% Define the Galois field matrix in GF(2^3)
data = load('galois_matrices.mat');
G = data.G3;

% Generating codewords for every k message bits
c = galois_mat_mul(G, message, primitive_poly, q); 
c_without_noise = c;
c_bits = reshape(dec2bin(c, k) - '0', [], 1); % Converting codewords to a binary sequence

% Constructing Parity Check Matrix
p_mat = G(k+1:n, :); % (n-k)*k matrix
p_mat_inverse = p_mat; % Additive inverse of an element is the same as itself for GF(2)
i2_mat = eye(n-k);
H = [p_mat_inverse i2_mat];

%% Adding Noise
% Adding noise to the transmitted message
noise = zeros(n,num_cols);
for i=1:num_cols
    noise_c = zeros(n,1);
    noise_pos = randperm(n,t);      % Choose 2 values from 1 to 7
    noise_val = randperm(n+1,t)-1;  % Choose 2 values from 0 to 7
    noise_c(noise_pos) = noise_val;
    noise(:,i) = noise_c;
end
c = galois_mat_add(c,noise,q);

%% Error Vectors
% Set of all possible error vectors
err_vectors = [];

% m represents number of possible errors at a time
for num_of_errors=1:t
    positions = nchoosek(1:n,num_of_errors);              % Each row represents a set of indices
    values = combnk_with_replacement(1:n, num_of_errors); % Each row represents a set of values that can be taken
    
    for l_ind_3 = 1:size(positions,1)
        for l_ind_4 = 1:size(values,1)
            vec = zeros(1,n);
            vec(positions(l_ind_3,:)) = values(l_ind_4,:);
            err_vectors = [err_vectors; vec];
        end
    end
end

%% Lookup Table
% Mapping between e and H.e'
lookup_table = containers.Map('KeyType', 'char', 'ValueType', 'any');

for l_ind_5 = 1:size(err_vectors,1)
    e = err_vectors(l_ind_5,:);          % Error Vector
    e = e';
    % e = gf(e,3);                   
    % syndrome = H*e';                     % Parity Matrix * Error Vector
    syndrome = galois_mat_mul(H, e, primitive_poly, q);
    
    % Convert numeric array to a string key using sprintf
    key = sprintf('%d ', syndrome);
    key = strtrim(key); % Remove any extra whitespace
    
    if isKey(lookup_table, key)    % Check if key already exists in lookup table
        current_best = lookup_table(key);  % If the key exists, update if the new e has fewer non-zero elements
        
        if nnz(e) < nnz(current_best)      % nnz - No: of non zero elements
            lookup_table(key) = e;
        end
    else
        lookup_table(key) = e;     % Otherwise, add the new key and error vector
    end
end

%% Syndrome Decoding
% Calculate syndrome for each codeword
s = galois_mat_mul(H,c,primitive_poly,q);
er = zeros(n,num_cols);

for i=1:num_cols
    s_c = s(:, i);
    
    % Convert syndrome to a string
    s_c_key = sprintf('%d ', s_c);
    s_c_key = strtrim(s_c_key);

    if(all(s_c_key ~= 0)) % If syndrome is not equal to 0
        if isKey(lookup_table, s_c_key)
            er_c = lookup_table(s_c_key);
            er_c = er_c';
            er(:,i) = er_c;
        end
    end
end

corrected_c = galois_mat_add(c,er,q); % Addition and Subtraction are the same for GF(2)

% Decode the codewords
rec_message = corrected_c(1:k,:);