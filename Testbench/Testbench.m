% Simulation Parameters

% Fixed Point parameters---------------------------------------------------
% Here wordLength denotes the total number of bits for each coefficient and
% fractionLength denotes how many of these bits are for fractional part.

wordLength = 12;
fractionLength = 9;

% Matched filter parameters (Square Root Raised Cosine Pulse)--------------

roll_off = 0.4;  
truncation_length = 10;
oversampling_factor = 8;
pulse_type = 'normalized'; 

% Signal parameters--------------------------------------------------------

symbol_rate = 4.8e3;
N = 2048;             % Number of symbols
M = 4;                % Number of constellation points for M-PSK modulation

% Channel parameters-------------------------------------------------------

Eb_N0 = 2; 
% Transmitter End

% Message signal-----------------------------------------------------------
message_signal = randi([0,1],2*N,1);           % Generating a random input sequence
modulated_signal = qpskMod1(message_signal);   % Message signal after QPSK modulation
%scatterplot(modulated_signal);

% Pulse Shaping Filter-----------------------------------------------------
srrc = srrcGen(roll_off,symbol_rate,truncation_length,oversampling_factor,pulse_type,wordLength,fractionLength);

% Pulse shaping to generate transimtted signal ----------------------------
transmitted_signal = pulseShaper(modulated_signal,oversampling_factor,srrc,wordLength,fractionLength);
% Channel simulation

% Introducing AWGN --------------------------------------------------------
received_sig = awgn_channel(transmitted_signal,Eb_N0,oversampling_factor,wordLength,fractionLength);
% Receiver End

% Matched Filter-----------------------------------------------------------
received_message = matchedFilter(received_sig,oversampling_factor,srrc,wordLength,fractionLength);
%scatterplot(double(received_message));

% Demodulation and BER calculation-----------------------------------------
output_message = qpskDemod1(received_message);
bit_error_rate = sum(output_message ~= message_signal)/(2*N)