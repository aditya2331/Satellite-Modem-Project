% The following code analyses the BER performance for an AWGN channel, with
% filter coefficients in fixed point. Monte Carlo method is used for 
% comparing the performance with theorotical output. 
%%
alp = 0.4;                         % roll-off factor
burst = 10;                        % number of bursts
EbN0_db = 0:1:8;                   % Range of SNR
EbN0_n = 10.^(EbN0_db/10);    
BER_t = qfunc(sqrt(2*EbN0_n));     % Calculating theorotical values of BER 
BER = zeros(burst,length(BER_t));
sym_rate = 4.8e3;                  % symbol rate
k = 8;                             % oversampling_factor
M = 4;                             % M-PSK modulation
N = 2048;                          % No. of symbols
%%
%Fixed point Parameters
wordLength = 12;         
fractionLength = 9;
%Generating SRRC pulse
srrc = srrcGen(alp,sym_rate,10,k,'normalized',wordLength,fractionLength);
%%
tic
parfor j = 1:1:burst
    inp_signal = randi([0,1],2*N,1);           % generating a random input sequence
    txSig = qpskMod1(inp_signal);              % qpsk modulation, in fixed point
    filt_inp = upsample(txSig,k);              % upsampling the modulated signal, to account for oversampling
    int_sig1 = conv(filt_inp,srrc,'same');     % pulse shaping at transmitter
    for i = 1:1:9
        % simulating AWGN channel, 12-bit fixed point
        int_sig2 = awgn_channel(int_sig1,EbN0_n(i),k,wordLength,fractionLength);
        % Output of the matched filter
        filt_out = fi(conv(int_sig2,srrc,'same'), 1, wordLength, fractionLength, 'SumMode', 'SpecifyPrecision', 'SumWordLength', wordLength, 'SumFractionLength', fractionLength, 'ProductMode', 'SpecifyPrecision', 'ProductWordLength', wordLength, 'ProductFractionLength', fractionLength);
        rxSig = downsample(filt_out,k);
        % Demodulating the received signal
        out_signal = qpskDemod1(rxSig);
        % BER calculation
        BER(j,i) = sum(out_signal ~= inp_signal)/(2*N);
    end
end
toc
BER_f = mean(BER);   % Average BER over all the bursts
%%
figure;  % Comparing the simulation performance with theorotical behaviour
semilogy(EbN0_db, BER_t, 'r', 'LineWidth', 2); % Theoretical BER
hold on;
semilogy(EbN0_db, BER_f, 'b', 'LineWidth', 2);   % Simulated BER
hold off
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend('Theoretical BER','Simulated BER with normalised SRRC and 1-2-9') %,'Simulated BER with floating point');
title('BER Performance of QPSK over AWGN Channel');
grid on;