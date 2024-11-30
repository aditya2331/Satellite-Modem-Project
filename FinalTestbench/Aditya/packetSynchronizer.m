function start_sampling_index = packetSynchronizer(preamble_correlation, matched_preamble_length, threshold, samplingOffset)
    
    if ~exist('samplingOffset','var')
        samplingOffset = 0;
    end

    [maxcorr, maxindex] = max(preamble_correlation);
    if maxcorr > threshold
        start_sampling_index = maxindex + matched_preamble_length + samplingOffset;
    else
        start_sampling_index = 0;
        disp('start not found, looking at next packet');
    end
    
    