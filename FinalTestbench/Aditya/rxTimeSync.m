function [sampled_signal, sampled_pilot] = rxTimeSync(rx_data, pilot, packet_size, threshold_factor, mismatch_tolerance, data_sampling_range, pilot_sampling_range, offset)


peaks = zeros(3, 1);                    % Unsigned integer + fractional
indices = [0; 3*mismatch_tolerance; 0]; % Unsigned integer 
thresholds = zeros(3, 1);               % Unsigned integer + fractional
start_at = 0;               % Unsigned integer 
sampling_index = 1;         % Unsigned integer 
is_synced = 0;              % Bit
packets_2extract = 0;       % Unsigned integer 

data_window = zeros(3*packet_size + offset, 1);


sampled_pilot = [];
sampled_signal = [];


% Correlator - output of correlator fills the following window with 1
%              packet size worth of samples at a time
cross_corr = conv(rx_data, pilot);
abs_cross_corr = abs(cross_corr);

plot(abs_cross_corr);
rx_data = [rx_data; zeros(offset,1)]; 


while (start_at+1)*packet_size <= numel(abs_cross_corr)

    % Buffer to contain:
    % correlation window of length = packet_size and 
    %     rx_data window of length = 3*packet_size + offset
    % Slide in (packet_size) number of samples in each iteration
    range = (1+start_at*packet_size):((start_at+1)*packet_size);
    corr_window = abs_cross_corr(range); 
    data_window = [data_window((packet_size+1):(3*packet_size + offset)); rx_data(range)];



    % Find peaks in each of the 3 blocks
    peaks(1:2) = peaks(2:3);
    indices(1:2) = indices(2:3);
    thresholds(1:2) = thresholds(2:3);
    [peaks(3), indices(3), thresholds(3)] = findPeak(corr_window, threshold_factor); %((2*packet_size+1):3*packet_size)


    
    % Check if preambles are found in all 3 blocks
    [is_synced, packets_2extract] = checkSync(is_synced, packets_2extract, peaks, indices, thresholds, mismatch_tolerance);



    % If preamble is found, then update the sampling point
    if is_synced == 1
        sampling_index = round(mean(indices));
    end



    % Extract data from packets
    if packets_2extract > 0
        sampled_pilot = [sampled_pilot;  data_window(sampling_index + offset + pilot_sampling_range)];
        sampled_signal = [sampled_signal; data_window(sampling_index + offset + data_sampling_range)];
        packets_2extract = packets_2extract - 1;

        fprintf('One packet extracted at %d and %d more left in buffer\n', sampling_index + (start_at-2)*packet_size, packets_2extract);

    end



    % Slide the window by one packet
    start_at = start_at + 1;
end


end




%----------------------------

function [peak, index, threshold] = findPeak(signal, threshold_factor)
    [peak, index] = max(signal);
    threshold = threshold_factor*mean(signal);
end

%----------------------------


function [is_synced, packets_2extract] = checkSync(is_synced, packets_2extract, peaks, indices, thresholds, mismatch_tolerance)
    if is_synced == 1 
        if abs(indices(3) - indices(2)) > mismatch_tolerance %|| peaks(3) < thresholds(3)
            is_synced = 0; % Sync lost
        else
            packets_2extract = packets_2extract + 1;
        end
    elseif peaks(1) >= thresholds(1) && peaks(2) >= thresholds(2) && peaks(3) >= thresholds(3) && ...
            abs(indices(2) - indices(1)) <= mismatch_tolerance && abs(indices(3) - indices(2)) <= mismatch_tolerance
        is_synced = 1;  % Sync'ed
        packets_2extract = 3;
    else
        is_synced = 0;  % No sync
    end
end

%----------------------------

