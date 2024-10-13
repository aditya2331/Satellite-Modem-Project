clear;
close all;
%% Parameters
% RS Code parameters
q = 8; % GF(2^q)-field
t = 16; % Error correction capability
n = 2^q -1; % Bits/codeword 
k = n - 2*t; % Bits/message

% Loading Addition and Multiplication Tables
data = load("Tables.mat");
add_table = data.(sprintf('AT%d',q));
mul_table = data.(sprintf('MT%d',q));

%% Finding all possible primitive elements
ff = 0:1:2^q -1;  % All elements in GF(2^q)
ff_star = ff(1,2:end);       % All elements in GF(2^q) except 0
factor_n = find_factors(n);  
ord_element = zeros(1,length(ff_star));

for i = 1:length(ff_star)
    for j = 1:2^q-1
        % if ff_star(1,i)^j == 1
        if gf_ele_exponent(ff_star(1,i),j,q) == 1
            ord_element(1,i) = j;
            break;
        end
    end
end

% Putting All primitive element into a list
% prim_ele = gf(zeros(1,0),q,primitive_poly);
prim_ele = zeros(1,0); % Empty Matrix
for i = 1:length(ord_element)
    if ord_element(1,i)==n
        prim_ele(1,end+1) = ff_star(1,i);
    end
end

%% Power Table;
% alpha - Chosen primitive Element
alpha = prim_ele(1,1);
power_ele = zeros(1,2^q);

for i = 1:2^q
    % ele = alpha^(i-1);
    ele = gf_ele_exponent(alpha,i-1,q);  % Galois Element Exponentiation
    power_ele(1,i) = ele;
end

power = 0:1:n;
% power_ele = power_ele.x;
power_to_ele_dict = containers.Map(power,power_ele);
ele_to_power_dict = containers.Map(power_ele,power);

%% Generating Message, Codeword
% Random message
len = 1*k;
bitstream = randi([0 n], 1, len);

% Sampling k bits and creating column vectors
num_cols = len/k;
message = reshape(bitstream, k, []);
inp_bits = reshape(dec2bin(message,k)-'0', [], 1);

% Define the Generator matrix in GF(2^q)
data_g = load('galois_matrices.mat');
G = data_g.(sprintf('G%d',q));

% Generating codewords for every k message bits
c = gf_mat_mul(G, message, q); 
c_without_noise = c';
c_bits = reshape(dec2bin(c, k) - '0', [], 1); % Converting codewords to a binary sequence
message = message';

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
c = gf_mat_add(c,noise,q);
noise = noise';
c = c';

%% Decoding
[corrected_c,flag] = rsDecode(c,alpha,n,t,q,power_to_ele_dict,ele_to_power_dict,add_table);
corrected_message = corrected_c(1:q);