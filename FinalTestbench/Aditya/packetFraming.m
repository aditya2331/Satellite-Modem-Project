function framed_data = packetFraming(modulated_symbols, modulated_preamble, modulated_preamble_size, packet_size)
    
    if iscolumn(modulated_symbols)
        % Transpose to make it a row vector
        modulated_symbols = modulated_symbols.';
    end
    if iscolumn(modulated_preamble)
        % Transpose to make it a row vector
        modulated_preamble = modulated_preamble.';
    end

    data_length = length(modulated_symbols);
    if mod(data_length, packet_size - modulated_preamble_size) ~= 0
        error('number of data symbols not divisible by packet_size - preamble_size');
    elseif length(modulated_preamble) ~= modulated_preamble_size
       error('wrong modulated_preamble_size');
    end
    num_packets = data_length/(packet_size - modulated_preamble_size);
    framed_data = zeros(1, packet_size*num_packets);
    for packet = 1:num_packets
        framed_data((packet-1)*(packet_size) + 1 : packet*(packet_size)) = [modulated_preamble modulated_symbols((packet-1)*(packet_size-modulated_preamble_size) + 1 : packet*(packet_size-modulated_preamble_size))]; % remember to use modulated preamble
    end

    % Ensure the output is a column vector
    framed_data = framed_data.';



