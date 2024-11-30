function rx_downsampled_packet = downSample(rx_match_filtered, current_pkt_pointer, sampling_index, samplesPerSymbol, data_symbols_per_packet, coderate)
  
    if sampling_index <= 0
        error('negative sampling_index')
    elseif (current_pkt_pointer + sampling_index - 1 + (coderate*data_symbols_per_packet - 1) * samplesPerSymbol) > length(rx_match_filtered) % end of signal reached
        rx_downsampled_packet = 0; 
        disp('end of signal reached')
    else
        rx_downsampled_packet = rx_match_filtered(current_pkt_pointer + sampling_index - 1:samplesPerSymbol:current_pkt_pointer + sampling_index - 1 + (coderate*data_symbols_per_packet - 1) * samplesPerSymbol);
    end