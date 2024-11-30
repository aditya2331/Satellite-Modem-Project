function CC_encoded_bits = CCencoder(input_bits, g1, g2)
    % % Define generator polynomials (octal representation)
    % g1 = [1 1 1 1 0 0 1];  % 171 in octal
    % g2 = [1 0 1 1 0 1 1];  % 133 in octal
    % Check if the input is a column vector
    if iscolumn(input_bits)
        % Transpose to make it a row vector
        input_bits = input_bits.';
    end
    if iscolumn(g1)
        % Transpose to make it a row vector
        g1 = g1.';
    end
    if iscolumn(g2)
        % Transpose to make it a row vector
        g2 = g2.';
    end
    if length(g1) ~= length(g2)
        error('g1 and g2 must have same length');
    end
    
    constraint_length = length(g1);

    % Initialize the encoder
    register = zeros(1, constraint_length - 1);  % 6 memory elements for constraint length 7

    % Preallocate output
    CC_encoded_bits = zeros(1, 2*length(input_bits));

    % Encoding process
    for i = 1:length(input_bits)
        
        % Generate two output bits
        output1 = mod(sum([input_bits(i) register] .* g1), 2);
        output2 = mod(sum([input_bits(i) register] .* g2), 2);

        % Shift register
        register = [input_bits(i) register(1:5)];
        % disp(register)
        
        % Store encoded bits
        CC_encoded_bits(2*i-1) = output1;
        CC_encoded_bits(2*i) = output2;
    end
    
    % Ensure the output is a column vector
    CC_encoded_bits = CC_encoded_bits.';

end