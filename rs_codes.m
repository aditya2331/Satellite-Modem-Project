% RS Code parameters
k = 3; % Bits/message
n = 7; % Bits/codeword
t = 2; % Error correction capability

% Random message
length = 30;
bitstream = randi([0 n], 1, length);

% Sampling k bits and creating column vectors
num_cols = length/k;
message = zeros(k, num_cols);
for i = 1:num_cols
    start_idx = (i-1)*k + 1;
    end_idx = i*k;
    bit_chunk = bitstream(start_idx:end_idx).';  % Transpose to form a column vector
    message(:,i) = bit_chunk;
end

message = gf(message, 3);

% Galois Field Elements
% g = gf(0:n,3); 
% Generator Matrix   
% g_prime = g(g ~= 0);
% G = gf(zeros(n,k),3);
% for i = 1:n
%    for j = 1:k
%        G(i, j) = g_prime(i) ^ (j-1);
%    end
% end

% Define the Galois field matrix in GF(2^3)
G = [1 0 0 6 1 6 7;
     0 1 0 4 1 5 5;
     0 0 1 3 1 2 3];

G = gf(G, 3);
G = G.';

% Generating codewords for every k message bits
c = G*message;
c_without_noise = c;

% Adding noise to the transmitted message
noise = zeros(n,num_cols);
for i=1:num_cols
    noise_c = zeros(n,1);
    noise_pos = randperm(n,t);      % Choose 2 values from 1 to 7
    noise_val = randperm(n+1,t)-1;  % Choose 2 values from 0 to 7
    noise_c(noise_pos) = noise_val;
    noise(:,i) = noise_c;
end
noise = gf(noise,3);
c = c+noise;

% Constructing Parity Check Matrix
p_mat = G(k+1:n, :); % (n-k)*k matrix
p_mat_inverse = -p_mat;
i2_mat = eye(n-k);
H = [p_mat_inverse i2_mat];
% K = H*G;

% Create a set of all possible error vectors
err_vectors = [];

% Custom function
function result = combnk_with_replacement(v, k)
    result = [];
    
    % Total number of possible combinations
    num_combinations = numel(v)^k;
    
    % Generate all combinations
    for i = 1:num_combinations
        % Convert index to combination
        index = i - 1;
        combination = zeros(1, k);
        for j = k:-1:1
            combination(j) = v(mod(index, numel(v)) + 1);
            index = floor(index / numel(v));
        end
        result = [result; combination];
    end
end

% m represents number of possible errors at a time
for m=1:t
    positions = nchoosek(1:n,m);     %  Each row represents a set of indices
    values = combnk_with_replacement(1:n, m); % Each row represents a set of values that can be taken
    
    for i = 1:size(positions,1)
        for j = 1:size(values,1)
            vec = zeros(1,n);
            vec(positions(i,:)) = values(j,:);
            err_vectors = [err_vectors; vec];
        end
    end
end

% Create a mapping between e and H.e'
lookup_table = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i = 1:size(err_vectors,1)
    e = err_vectors(i,:);          % Error Vector
    e = gf(e,3);                   
    syndrome = H*e';               % Parity Matrix * Error Vector
    
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

% Lookup table content
table_keys = lookup_table.keys;
table_values = lookup_table.values;

% Calculate syndrome for each codeword
s = H*c;
er = gf(zeros(n,num_cols),3);

for i=1:num_cols
    s_c = s(:, i);
    
    % Convert syndrome to a string
    s_c_key = sprintf('%d ', s_c.x);
    s_c_key = strtrim(s_c_key);

    if(all(s_c_key ~= 0)) % If syndrome is not equal to 0
        if isKey(lookup_table, s_c_key)
            er_c = lookup_table(s_c_key);
            er_c = gf(er_c,3);
            er_c = er_c';
            er(:,i) = er_c;
        end
    end
end

corrected_c = c-er;

% Decode the codewords
rec_message = corrected_c(1:k,:);