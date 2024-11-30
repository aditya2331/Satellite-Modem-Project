function h_est = channelEstimatorModified(rx_sampled, preamble_len, ch_len)
    % rx_sampled(column vec): inital 26 bits must be the convolution of channel with
    % the pilot sequence
    % ch_len : length of the ISI channel 
    % premable_len : length of preamble sequence, 26 in our case
    % preamble sequence

    preamble_seqs = [0 0 1 0 0 1 0 1 1 1 0 0 0 0 1 0 0 0 1 0 0 1 0 1 1 1; 0 0 1 0 1 1 0 1 1 1 0 1 1 1 1 0 0 0 1 0 1 1 0 1 1 1; 0 1 0 0 0 0 1 1 1 0 1 1 1 0 1 0 0 1 0 0 0 0 1 1 1 0; 0 1 0 0 0 1 1 1 1 0 1 1 0 1 0 0 0 1 0 0 0 1 1 1 1 0; 0 0 0 1 1 0 1 0 1 1 1 0 0 1 0 0 0 0 0 1 1 0 1 0 1 1; 0 1 0 0 1 1 1 0 1 0 1 1 0 0 0 0 0 1 0 0 1 1 1 0 1 0; 1 0 1 0 0 1 1 1 1 1 0 1 1 0 0 0 1 0 1 0 0 1 1 1 1 1; 1 1 1 0 1 1 1 1 0 0 0 1 0 0 1 0 1 1 1 0 1 1 1 1 0 0];
    % Preamble sequence chosen for channel estimation
    g_seq = circshift(preamble_seqs(2,:),19);
    for i=1:length(g_seq)
        if(g_seq(1,i)==0) 
            g_seq(1,i) = -1;
        end
    end
    % pilot matrix computation
    tp_mat = zeros([preamble_len-ch_len+1,ch_len]);
    for i = 1:length(tp_mat(:,1))
        tp_mat(i,:) = flip(g_seq(1,i:ch_len -1 +i));
    end 
    pilot_matrix = pinv(tp_mat);
    rx_sampled = transpose(rx_sampled);
    h_est = transpose(pilot_matrix*transpose(rx_sampled(1,ch_len:preamble_len))); % 
    h_est = transpose(h_est);
end
