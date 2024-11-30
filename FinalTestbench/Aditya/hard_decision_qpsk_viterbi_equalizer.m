function [decoded_symbols] = hard_decision_qpsk_viterbi_equalizer(received_signal_downsampled, packet_size, channel_taps, known_preamble_bits, traceback_length)
    % QPAM Viterbi Equalizer with limited trellis storage
    % Inputs:
    %   received_signal: Complex received signal samples
    %   channel_taps: Known channel impulse response of form [h(0) h(1) ... h(12)]
    %   known_initial_state: Initial state sequence
    %   traceback_length: Length of trellis history to maintain

    if iscolumn(received_signal_downsampled)
        % Transpose to make it a row vector
        received_signal_downsampled = received_signal_downsampled.';
    end
    if iscolumn(channel_taps)
        % Transpose to make it a row vector
        channel_taps = channel_taps.';
    end

    
    % Set default traceback length if not specified
    if nargin < 4
        traceback_length = 15; % Default traceback length
    end
    
    % QPAM constellation
    constellation = [1+1i, 1-1i, -1+1i, -1-1i] / sqrt(2);
    
    % Parameters
    L = length(channel_taps);
    N = length(received_signal_downsampled);
    num_states = 4^(L-1);  % Number of states based on channel memory
    
    % % Initialize trellis structures with limited length
    % path_metrics = Inf(num_states, traceback_length + 1);
    % branch_metrics = Inf(num_states, 4, traceback_length);
    % survivors = zeros(num_states, traceback_length);
    decoded_symbols = zeros(1, N);
    
    % Calculate initial state from known preamble
    preamble_state_bits = known_preamble_bits(end:-1:end - L + 2);

    % Initialize starting state
    inistate = 0;
    for i = length(preamble_state_bits):-1:1
        inistate = inistate + preamble_state_bits(i)*4^(i-1);
    end
    
    % State transition matrix for QPAM symbol transitions
    state_transitions = zeros(num_states/4, 4);
    for state = 0:length(state_transitions)-1
        state_symbols = zeros(1, L-2);
        temp_state = state;
        for i = L-2:-1:1
            state_symbols(i) = mod(temp_state, 4);
            temp_state = floor(temp_state/4);
        end

        for new_symbol = 0:3
            next_symbols = [new_symbol, state_symbols];
            next_state = 0;
            for i = 1:L-1
                next_state = next_state + next_symbols(i) * 4^(L-1-i);
            end
            state_transitions(state+1, new_symbol+1) = next_state;
        end
    end

    coded_data_packet_size = packet_size - length(known_preamble_bits);

    for j = 1:length(received_signal_downsampled)/coded_data_packet_size
        
        % Initialize trellis structures with limited length
        path_metrics = Inf(num_states, traceback_length + 1);
        % branch_metrics = Inf(num_states, 4, traceback_length);
        survivors = zeros(num_states, traceback_length);
        path_metrics(inistate + 1, 1) = 0;


        decoded_packet = zeros(1, coded_data_packet_size);
        received_packet = received_signal_downsampled((j-1)*coded_data_packet_size + 1: j*coded_data_packet_size);

        % Viterbi algorithm with sliding window
        for n = 1:coded_data_packet_size
            % Current position in the trellis window
            trellis_idx = mod(n-1, traceback_length) + 1;
            next_trellis_idx = mod(n, traceback_length) + 1;
            
            % Reset next column of path metrics
            path_metrics(:, next_trellis_idx) = Inf;
            
            for current_state = 0:num_states-1
                % Extract QPAM symbols from state
                state_symbols = zeros(1, L-1);
                temp_state = current_state;
                for i = L-1:-1:1
                    state_symbols(i) = mod(temp_state, 4);
                    temp_state = floor(temp_state/4);
                end
    
                state_qpsk = constellation(state_symbols + 1);
    
                if n < L     
                    state_qpsk(n:L-1) = 2*preamble_state_bits(n:L-1) - 1;
                end
                
                % Try each possible input symbol
                for sym_idx = 1:4
                    next_symbol = constellation(sym_idx);
                    symbols = [next_symbol, state_qpsk];
                    
                    % Calculate expected output
                    expected = sum(symbols .* channel_taps);
                    
                    % Branch metric
                    branch_metric = abs(received_packet(n) - expected)^2;
                    % branch_metrics(current_state+1, sym_idx, trellis_idx) = branch_metric;
                    
                    % Next state
                    first3bits = floor(current_state/4);
                    next_state = state_transitions(first3bits+1, sym_idx);
                    
                    % Update path metric
                    new_metric = path_metrics(current_state+1, trellis_idx) + branch_metric;
                    
                    if new_metric < path_metrics(next_state+1, next_trellis_idx)
                        path_metrics(next_state+1, next_trellis_idx) = new_metric;
                        survivors(next_state+1, trellis_idx) = current_state;
                    end
                end
            end
            
            % Perform traceback and decode if we have enough history
            if n >= traceback_length
                % Find best ending state
                [~, current_state] = min(path_metrics(:, next_trellis_idx));
                current_state = current_state - 1;
                
                % Traceback through the window
                for k = 0:traceback_length-1
                    idx = mod(n-k-1, traceback_length) + 1;
                    decoded_idx = (current_state - mod(current_state, 4^(L-2)))/4^(L-2);
                    if k == traceback_length-1
                        % Only save the oldest decoded symbol
                        decoded_packet(n-traceback_length+1) = constellation(decoded_idx + 1);
                    end
                    current_state = survivors(current_state+1, idx);
                end
            end
        end
        
        % Process final symbols
        final_symbols = zeros(1, traceback_length-1);
        [~, current_state] = min(path_metrics(:, next_trellis_idx));
        current_state = current_state - 1;
        
        for k = 0:traceback_length-2
            idx = mod(coded_data_packet_size-k-1, traceback_length) + 1;
            decoded_idx = (current_state - mod(current_state, 4^(L-2)))/4^(L-2);
            final_symbols(traceback_length-1-k) = constellation(decoded_idx + 1);
            current_state = survivors(current_state+1, idx);
        end
        
        % Add final symbols to output
        decoded_packet(coded_data_packet_size-traceback_length+2:coded_data_packet_size) = final_symbols;
        decoded_symbols((j-1)*coded_data_packet_size + 1: j*coded_data_packet_size) = decoded_packet;
    end
    % Ensure the output is a column vector
    decoded_symbols = decoded_symbols.';
end