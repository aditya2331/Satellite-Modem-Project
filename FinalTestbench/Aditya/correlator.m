function rx_correlation = correlator(rx_match_filtered_packet, Corr_FIR_coefficients)
    filter_length = length(Corr_FIR_coefficients);
    data_length = length(rx_match_filtered_packet);
    rx_correlation = zeros(1, data_length - filter_length + 1); % assuming packetwise correlation
    for m = 1:(data_length - filter_length + 1)
        sumProduct = 0;
        for n = 1:filter_length
            sumProduct = sumProduct + real(rx_match_filtered_packet(m + n - 1)) * Corr_FIR_coefficients(n);
        end
        rx_correlation(m) = rx_correlation(m) + abs(sumProduct);
    end
