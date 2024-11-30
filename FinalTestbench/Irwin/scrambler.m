
%% (De)Scrambler functional modules(Ordinary version)(Not Persistent version)

function scrambled_bits = scrambler(bits)
    % Input  : bits (column vector)
    % Output : scrambled_bits (column vector) 
    bits = transpose(bits);
    sr      = [1 0 0 1 1 0 1 1 0 0 1 1 0 0 1 1 0 1 1 0]; % shift register
    count   = ones(1,5); % counter
    scrambled_bits = zeros(1,length(bits));
    N_clk = length(bits);
    for i = 1:N_clk
        % performs scrambling bit-by-bit
        [scrambled_bits(1,i),sr,count] = scrambler_serial(bits(1,i),sr,count);
    end
    scrambled_bits = transpose(scrambled_bits); % output scrambled bits
end



%% Supporting Functions

function [out,sr_out,count_out] = scrambler_serial(inp_bits,sr,count)

    % Monopulser ON condition or Reset Condition
    if(bitxor(sr(1,1),sr(1,9))==1)
        count = ones(1,5);
    end


    % Combinational Logic
    nor_1 = ~(count(1,1)|~count(1,2)|~count(1,3)|~count(1,4)|~count(1,5));
    xnor_1 = ~(bitxor(sr(1,3),sr(1,20)));
    xnor_2 = ~(bitxor(xnor_1,nor_1));
    inp = ~bitxor(inp_bits, xnor_2);
    % Output
    out = inp;

    % Clock Update
    sr_out    = shift_register(sr,inp);
    count_out = counter_rev(count,5);
end

% Shift Register Update Function
function [sr_out] = shift_register(sr_in,inp)
    % Right Shift the contents of the Register
    sr_out = [inp sr_in(1,1:end-1)];
end

function [counter_out] = counter_rev(counter_in, bit_len)
    max_val = 2^bit_len;

    % Convert binary vector to decimal manually
    c_val_in = sum(counter_in .* (2.^(bit_len-1:-1:0)));

    % Calculate the next counter value (decrement with modulo)
    c_val_out = mod(c_val_in - 1, max_val);

    % Convert decimal back to binary vector
    counter_out = bitget(c_val_out, bit_len:-1:1);
end