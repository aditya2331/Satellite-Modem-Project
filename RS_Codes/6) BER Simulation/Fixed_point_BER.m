clear;
close all;
%% Parameters
% RS Code parameters
rs_q = 8; % GF(2^q)-field
rs_t = 16; % Error correction capability
rs_n = 2^rs_q - 1; % Bits/codeword
rs_k = rs_n - 2*rs_t; % Bits/message
data_g = load('galois_matrices.mat');
rs_G = data_g.(sprintf('G%d',rs_q));

% Modulation Parameters
alp = 0.4;          % roll-off factor
burst = 1;          % numsber of bursts
sym_rate = 4.8e3;
k = 8;               % oversampling_factor
M = 4;               % M-PSK modulation

% Fixed point Parameters
wordLength = 12;         
fractionLength = 9;
srrc = srrcGen(alp,sym_rate,10,k,'normalized',wordLength,fractionLength);  % Generating SRRC pulse

% EbN0 Range
EbN0_db = 1;
EbN0_db = EbN0_db - 10*log10(rs_n/rs_k);
EbN0_n = 10.^(EbN0_db/10);

EbN0_db_scaled = EbN0_db + 10*log10(rs_n/rs_k);
EbN0_n_scaled = 10.^(EbN0_db_scaled/10);
BER_t = qfunc(sqrt(2*EbN0_n_scaled));
BER = zeros(burst,length(BER_t));

% Loading Addition and Multiplication Tables
data = load("Tables.mat");
add_table = data.(sprintf('AT%d',rs_q));
mul_table = data.(sprintf('MT%d',rs_q));

%% Finding all possible primitive elements
ff = 0:1:2^rs_q -1;  % All elements in GF(2^q)
ff_star = ff(1,2:end);       % All elements in GF(2^q) except 0
factor_n = find_factors(rs_n);  
ord_element = zeros(1,length(ff_star));

for i = 1:length(ff_star)
    for j = 1:2^rs_q-1
        % if ff_star(1,i)^j == 1
        if gf_ele_exponent(ff_star(1,i),j,rs_q) == 1
            ord_element(1,i) = j;
            break;
        end
    end
end

% Putting All primitive element into a list
% prim_ele = gf(zeros(1,0),q,primitive_poly);
prim_ele = zeros(1,0); % Empty Matrix
for i = 1:length(ord_element)
    if ord_element(1,i)==rs_n
        prim_ele(1,end+1) = ff_star(1,i);
    end
end

%% Power Table;
% alpha - Chosen primitive Element
alpha = prim_ele(1,1);
power_ele = zeros(1,2^rs_q);

for i = 1:2^rs_q
    % ele = alpha^(i-1);
    ele = gf_ele_exponent(alpha,i-1,rs_q);  % Galois Element Exponentiation
    power_ele(1,i) = ele;
end

power = 0:1:rs_n;
power_to_ele_dict = containers.Map(power,power_ele);
ele_to_power_dict = containers.Map(power_ele,power);

%%
tic
for j = 1:burst         
    len = 2*rs_k; % Length of i/p message
    message = randi([0 rs_n], 1, len); % Generating a random input sequence in GF(2^q)

    [message_bits, enc_inp_signal_bits] = rs_encoder(message,rs_k,rs_q,rs_G); % RS Encoding on the input message

    txSig = qpskMod1(enc_inp_signal_bits);              % QPSK Modulation, in fixed point
    filt_inp = upsample(txSig,k);
    int_sig1 = conv(filt_inp,srrc,'same');     % Pulse shaping at transmitter
    Ps = mean(double(int_sig1).*conj(double(int_sig1)))/2;  % Signal power
    for i = 1:1:1
        % simulating AWGN channel, 12-bit fixed point
        int_sig2 = awgn_channel(int_sig1,EbN0_n(i),Ps,k,wordLength,fractionLength);
        % Output of the matched filter
        filt_out = fi(conv(int_sig2,srrc,'same'), 1, wordLength, fractionLength, 'SumMode', 'SpecifyPrecision', 'SumWordLength', wordLength, 'SumFractionLength', fractionLength, 'ProductMode', 'SpecifyPrecision', 'ProductWordLength', wordLength, 'ProductFractionLength', fractionLength);
        rxSig = downsample(filt_out,k);
        demod_signal = qpskDemod1(rxSig);
        
        [out_codeword_bits, out_message_bits] = rs_decoder(demod_signal,rs_n,rs_k,rs_q,rs_t,alpha,power_to_ele_dict,ele_to_power_dict,add_table); % RS Decoding on the demodulated bits
        
        BER(j,i) = sum(out_message_bits ~= message_bits)/(len*rs_q);        
        % BER_rs(j,i) = sum(out_message_bits ~= message_bits)/(len*rs_q); % BER of i/p and o/p message with RS Coding
    end
end
toc
BER_f = mean(BER);
% BER_rs_f = mean(BER_rs);
%%
EbN0_db_graph = 0:1:10;
BER_t_graph = qfunc(sqrt(2*(10.^(EbN0_db_graph/10))));

figure;  % Comparing the simulation performance with theorotical behaviour
semilogy(EbN0_db_graph, BER_t_graph, '-ro', 'LineWidth', 1); % Theoretical BER
hold on;
semilogy(EbN0_db_scaled, BER, '-go', 'LineWidth', 1);      
% semilogy(EbN0_db_scaled, BER_rs, 'g', 'LineWidth', 2);   % Simulated BER with RS Encoding
hold off
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend('Theoretical BER', 'Simulated BER with RS Encoding') %,'Simulated BER with floating point');
title('BER Performance of QPSK over AWGN Channel');
grid on;