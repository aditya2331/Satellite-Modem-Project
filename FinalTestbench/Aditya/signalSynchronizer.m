function synchronized_downsampled_data = signalSynchronizer(rx_match_filtered, packetsForStartDetection, packetsForAveraging, syncFrequency, num_data_packets, upsampling_factor, matched_preamble, coded_packet_size, data_packet_size, threshold, coderate, samplingOffset)
    
    % Check if the input is a column vector
    if iscolumn(rx_match_filtered)
        % Transpose to make it a row vector
        rx_match_filtered = rx_match_filtered.';
    end
    % Check if the input is a column vector
    if iscolumn(matched_preamble)
        % Transpose to make it a row vector
        matched_preamble = matched_preamble.';
    end
    
    
    start_found = false;
    current_pkt_pointer = 1;
    matched_prb_len = length(matched_preamble);
    data_pkts_seen = 0;
    preamble_correlation = zeros(1, coded_packet_size*upsampling_factor - matched_prb_len + 1);
    synchronized_downsampled_data = zeros(1, num_data_packets*data_packet_size*coderate);

    while current_pkt_pointer + upsampling_factor*coded_packet_size - 1 <= length(rx_match_filtered)
        
        if start_found == false % search for start of burst
            avging_pkt_pointer = current_pkt_pointer;
            rx_match_filtered_packet = rx_match_filtered(1, avging_pkt_pointer : avging_pkt_pointer + upsampling_factor*coded_packet_size - 1);

            for x = 1: packetsForStartDetection
                preamble_correlation = preamble_correlation + correlator(rx_match_filtered_packet, matched_preamble);
                avging_pkt_pointer = avging_pkt_pointer + upsampling_factor*coded_packet_size;
                if avging_pkt_pointer + upsampling_factor*coded_packet_size - 1 <= length(rx_match_filtered)
                    rx_match_filtered_packet = rx_match_filtered(1, avging_pkt_pointer : avging_pkt_pointer + upsampling_factor*coded_packet_size - 1);
                else
                    break
                end
            end
            preamble_correlation = preamble_correlation/x;
            start_sampling_index = packetSynchronizer(preamble_correlation, matched_prb_len, threshold, samplingOffset);
            if start_sampling_index  ~= 0 % start of burst found
                start_found = true;
                disp('start of burst found');
                plot(preamble_correlation)
            else
                current_pkt_pointer = current_pkt_pointer + upsampling_factor*coded_packet_size;
            end 
        else % start synchronizing packets
            avging_pkt_pointer = current_pkt_pointer;
            rx_match_filtered_packet = rx_match_filtered(1, avging_pkt_pointer : avging_pkt_pointer + upsampling_factor*coded_packet_size - 1);
            for x = 1: packetsForAveraging
                % Loop to calculate the correlation at each point
                preamble_correlation = preamble_correlation + correlator(rx_match_filtered_packet, matched_preamble);
                avging_pkt_pointer = avging_pkt_pointer + upsampling_factor*coded_packet_size;
                if avging_pkt_pointer + upsampling_factor*coded_packet_size - 1 <= length(rx_match_filtered)
                    rx_match_filtered_packet = rx_match_filtered(1, avging_pkt_pointer : avging_pkt_pointer + upsampling_factor*coded_packet_size - 1);
                else
                    break
                end
            end
            preamble_correlation = preamble_correlation/x;
            sampling_index = packetSynchronizer(preamble_correlation, matched_prb_len, threshold, samplingOffset);
            if sampling_index  ~= 0 % start of burst found
                for y = data_pkts_seen:data_pkts_seen+syncFrequency-1
                    synchronized_downsampled_data((y*coderate*data_packet_size) + 1 : (y+1)*(coderate*data_packet_size)) = downSample(rx_match_filtered, current_pkt_pointer, sampling_index, upsampling_factor, data_packet_size, coderate);
                    current_pkt_pointer = current_pkt_pointer + upsampling_factor*coded_packet_size;
                end
                data_pkts_seen = data_pkts_seen + syncFrequency;
            else
                error('discontinuous burst, could not average, stopping')
            end
        end
    end

    % Ensure the output is a column vector
    synchronized_downsampled_data = synchronized_downsampled_data.';