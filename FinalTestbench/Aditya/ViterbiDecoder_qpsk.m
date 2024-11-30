function [decoded_bits, trellis, survivors] = ViterbiDecoder_qpsk(soft_equalized_received_symbols, traceback_len, g1, g2, symbol_mapping, num_quantization_bits)

    if iscolumn(soft_equalized_received_symbols)
        % Transpose to make it a row vector
        soft_equalized_received_symbols = soft_equalized_received_symbols.';
    end
    
    if length(g1) ~= length(g2)
        error('g1 and g2 must have same length');
    end
    
    constraint_length = length(g1);

    % Number of states in the trellis (2^6 for constraint length 7)
    num_states = 2^(constraint_length-1);

    % Initialize 
    trellis = inf(num_states, traceback_len + 1);
    trellis(1, 1) = 0;  % Start in state 0
    survivors = zeros(num_states, traceback_len + 1);
    
    if num_quantization_bits ~= 0
        quantization = true;
        % disp('quantization = true')
    else
        quantization = false;
    end

    decoded_bits = zeros(1, length(soft_equalized_received_symbols));

    output_dictionary = zeros([64 2]);
    prevstate_dictionary = zeros(64,2);
    for s = 0:num_states-1
        input = s >= 32;
        % Calculate expected output for transition from prev_states(1)
        register_0 = [bitget(s,5:-1:1), 0];
        expected_output_0 = [mod(sum([input register_0] .* g1), 2), mod(sum([input register_0] .* g2), 2)];
        register_1 = [bitget(s,5:-1:1), 1];
        expected_output_1 = [mod(sum([input register_1] .* g1), 2), mod(sum([input register_1] .* g2), 2)];
        prev_state_0 = 2*mod(s, 32);
        prev_state_1 = 2*mod(s, 32) + 1;
        output_dictionary(s+1, :) = [2*expected_output_0(1)+expected_output_0(2) 2*expected_output_1(1)+expected_output_1(2)];
        prevstate_dictionary(s+1, :) = [prev_state_0 prev_state_1];
    end

    % Decoding process
    if trellis(1, 1) == 0
        
        for t = 2:length(soft_equalized_received_symbols)+1
            for s = 0:num_states-1
    
                % Calculate expected output for transition from prev_states(1)
                expected_outputs = output_dictionary(s+1, :);

                prev_states = prevstate_dictionary(s+1, :);
       
                metric_0 = (abs(soft_equalized_received_symbols(t-1) - symbol_mapping(expected_outputs(1) + 1)))^2;
                metric_1 = (abs(soft_equalized_received_symbols(t-1) - symbol_mapping(expected_outputs(2) + 1)))^2;

                if quantization == true
                    metric_0 = generalQuantize(metric_0, 0, 30, num_quantization_bits);
                    metric_1 = generalQuantize(metric_1,0, 30, num_quantization_bits);
                end

                if (t > traceback_len + 1)
                    % Update trellis and survivors
                    path_metric_0 = metric_0 + trellis(prev_states(1)+1, end);
                    path_metric_1 = metric_1 + trellis(prev_states(2)+1, end);
        
                    if path_metric_0 <= path_metric_1
                        trellis(s+1, 1) = path_metric_0;
                        survivors(s+1,1) = prev_states(1);
                    else
                        trellis(s+1, 1) = path_metric_1;
                        survivors(s+1,1) = prev_states(2);
                    end
                else
                    % Update trellis and survivors
                    path_metric_0 = metric_0 + trellis(prev_states(1)+1, t-1);
                    path_metric_1 = metric_1 + trellis(prev_states(2)+1, t-1);
        
                    if path_metric_0 <= path_metric_1
                        trellis(s+1, t) = path_metric_0;
                        survivors(s+1, t) = prev_states(1);
                    else
                        trellis(s+1, t) = path_metric_1;
                        survivors(s+1, t) = prev_states(2);
                    end
                end
            end
    
            if (t > traceback_len + 1)
                [~, best_state] = min(trellis(:, end));
                temp_decoded_trellis = traceback_decoder(best_state, survivors);
                % disp(t - traceback_len-1)
                if length(soft_equalized_received_symbols) >  traceback_len
                    decoded_bits(t - traceback_len-1) =  temp_decoded_trellis(1,1);
                end
                % disp(trellis)
                trellis(:, end) = trellis(:, 1);
                % disp(trellis)
                survivors = circshift(survivors, -1, 2);
            end
 
        end
    end
    % Traceback and circ shift 
    if length(soft_equalized_received_symbols) >  traceback_len
        [~, best_state] = min(trellis(:, end));
        decoded_bits(end - traceback_len + 1 : end) =  traceback_decoder(best_state, survivors);
    else
        [~, best_state] = min(trellis(:, length(soft_equalized_received_symbols)+1));
        decoded_bits =  traceback_decoder(best_state, survivors(:, 1:length(soft_equalized_received_symbols)+1));
    end

    % Ensure the output is a column vector
    decoded_bits = decoded_bits.';
end

function traceback_decoded_bits = traceback_decoder(best_state, survivors)

    traceback_decoded_bits = zeros(1, length(survivors(1, :))-1);
   
    % for t = length(trellis(1,:)):-1:2
    for t = length(survivors(1,:)):-1:2
        % disp(best_state)
        prev_state = survivors(best_state, t); % (0 to 63)
        traceback_decoded_bits(t-1) = best_state-1 >= 32; % ??
        best_state = prev_state + 1;
    end

end