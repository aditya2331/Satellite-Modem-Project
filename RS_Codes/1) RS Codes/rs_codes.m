clear;
close all;
%% Generator Matrix, Codewords, Parity Check Matrix
% RS Code parameters
q = 3; % GF(2^q)-field
t = 2; % Error correction capability
n = 2^q -1; % Bits/codeword 
k = n - 2*t; % Bits/message

% Random message
len = 1*k;
bitstream = randi([0 n], 1, len);

% Sampling k bits and creating column vectors
num_cols = len/k;
message = reshape(bitstream, k, []);
message = gf(message, k);
inp_bits = reshape(dec2bin(message.x,k)-'0', [], 1);

% Define the Galois field matrix in GF(2^3)
data = load("galois_matrices.mat");
G = data.(sprintf('G%d',q));
G = gf(G, k);

% Generating codewords for every k message bits
c = G*message;
c_without_noise = c;
c_bits = reshape(dec2bin(c.x, k) - '0', [], 1); % Converting codewords to a binary sequence

% Constructing Parity Check Matrix
p_mat = G(k+1:n, :); % (n-k)*k matrix
p_mat_inverse = -p_mat;
i2_mat = eye(n-k);
H = [p_mat_inverse i2_mat];

%% Adding Noise
% Adding noise to the transmitted message
noise = zeros(n,num_cols);
for i=1:num_cols
    noise_c = zeros(n,1);
    noise_pos = randperm(n,t);      % Choose t values from 1 to n-1
    noise_val = randperm(n+1,t)-1;  % Choose t values from 0 to n-1
    noise_c(noise_pos) = noise_val;
    noise(:,i) = noise_c;
end
noise = gf(noise,k);
c = c+noise;

%% Error Vectors
% Set of all possible error vectors
err_vectors = [];

% m represents number of possible errors at a time
for num_of_errors=1:t
    positions = nchoosek(1:n,num_of_errors);              % Each row represents a set of indices
    values = comb_nk(1:n, num_of_errors); % Each row represents a set of values that can be taken
    
    for i = 1:size(positions,1)
        for j = 1:size(values,1)
            vec = zeros(1,n);
            vec(positions(i,:)) = values(j,:);
            err_vectors = [err_vectors; vec];
        end
    end
end

%% Lookup Table
% Mapping between e and H.e'
lookup_table = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i = 1:size(err_vectors,1)
    e = err_vectors(i,:);          % Error Vector
    e = gf(e,k);                   
    syndrome = H*e';                     % Parity Matrix * Error Vector
    
    % Convert numeric array to a string key using sprintf
    key = sprintf('%d ', syndrome.x);
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
s = H*c;
er = gf(zeros(n,num_cols),k);

for i=1:num_cols
    s_c = s(:, i);
    
    % Convert syndrome to a string
    s_c_key = sprintf('%d ', s_c.x);
    s_c_key = strtrim(s_c_key);

    if(all(s_c_key ~= 0)) % If syndrome is not equal to 0
        if isKey(lookup_table, s_c_key)
            er_c = lookup_table(s_c_key);
            er_c = gf(er_c,k);
            er_c = er_c';
            er(:,i) = er_c;
        end
    end
end

corrected_c = c-er;

% Decode the codewords
rec_message = corrected_c(1:k,:);