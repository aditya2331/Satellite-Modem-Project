function traceback_decoded_bits = traceback_decoder(best_state, survivors)

    % traceback_decoded_bits = zeros(1, length(trellis(1, :))-1);
    traceback_decoded_bits = zeros(1, length(survivors(1, :))-1);

    % [~, current_state] = min(trellis(:, end)); % TODO only using the end, so we need not store entire trellis
    % [~, current_state] = min(trellis); % TODO only using the end, so we need not store entire trellis
    % disp(current_state)
   
    % for t = length(trellis(1,:)):-1:2
    for t = length(survivors(1,:)):-1:2
        % disp(best_state)
        prev_state = survivors(best_state, t); % (0 to 63)
        traceback_decoded_bits(t-1) = best_state-1 >= 32; % ??
        best_state = prev_state + 1;
    end

end