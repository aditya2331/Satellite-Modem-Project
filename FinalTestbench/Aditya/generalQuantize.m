function q_val = generalQuantize(x, min_val, max_val, num_bits)
    % Quantizes a floating point value x into num_bits levels
    % between min_val and max_val

    % Ensure the value x is within the range
    x = min(max(x, min_val), max_val);

    % Compute the quantization step size
    step_size = (max_val - min_val) / (2^num_bits - 1);

    % Quantize the input value
    q_val = round((x - min_val) / step_size) * step_size + min_val;
end