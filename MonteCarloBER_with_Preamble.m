function qpskModemSimulation()
    % Parameters
    fixedPrecision = false;
    symbolRate = 4800; % 4.8 Ksymbols/sec
    samplesPerSymbol = 8; % Oversampling factor
    numSymbols = 1024;
    alpha = 0.4; % Roll-off factor for SRRC
    Eb_N0_dB = 0:2:12; % Eb/N0 range in dB
    numBursts = 1; % Number of bursts for averaging
    prb = [0 0 1 0 0 1 0 1 1 1 0 0 0 0 1 0 0 0 1 0 0 1 0 1 1 1];
    samplingOffsets = [0];

    % QPSK symbol mapping (Gray coded)
    symbolMapping = [1+1i, 1-1i, -1+1i, -1-1i] / sqrt(2);
    wordLength = 12;         % Total number of bits
    fractionLength = 7;      % Number of fractional bits
    
    if fixedPrecision == true
        symbolMapping = fi(symbolMapping, 1, wordLength, fractionLength, 'SumMode', 'SpecifyPrecision', 'SumWordLength', wordLength, 'SumFractionLength', fractionLength, 'ProductMode', 'SpecifyPrecision', 'ProductWordLength', wordLength, 'ProductFractionLength', fractionLength);
    end

    % Generate SRRC pulse
    [srrcPulse, t] = generateSRRC(symbolRate, samplesPerSymbol, alpha, fixedPrecision);

    % (a) Constellation plot for Eb/N0 = 6 dB
    % plotConstellation(symbolMapping, symbolRate, samplesPerSymbol, srrcPulse, 10, fixedPrecision);
    
    for i = 1:length(samplingOffsets)
        % (b) & (c) BER/SER simulation
        [ber, ser, EbN0] = simulateBERSER_preamble(prb,samplingOffsets(i), symbolMapping, symbolRate, samplesPerSymbol, srrcPulse, Eb_N0_dB, numBursts, numSymbols, fixedPrecision);
    
        % (d) Plot BER vs Eb/N0 and compare with analytical
        plotBERComparison(EbN0, ber, samplingOffsets(i));
    
        % Plot SER vs Es/N0
        % plotSER(EbN0 + 3, ser); % Es/N0 = Eb/N0 + 3dB for QPSK
    end
end

function [srrcPulse, t] = generateSRRC(symbolRate, samplesPerSymbol, alpha, fixedPrecision)
    Ts = 1 / symbolRate;
    t = -5*Ts:Ts/samplesPerSymbol:5*Ts;
    
    % Create fixed-point objects with specific precision
    % Parameters
    wordLength = 12;         % Total number of bits
    fractionLength = 7;      % Number of fractional bits
    srrcPulse = zeros(size(t));
    if fixedPrecision == true
        srrcPulse = fi(zeros(size(t)), 1, wordLength, fractionLength, 'SumMode', 'SpecifyPrecision', 'SumWordLength', wordLength, 'SumFractionLength', fractionLength, 'ProductMode', 'SpecifyPrecision', 'ProductWordLength', wordLength, 'ProductFractionLength', fractionLength);
    end
    
    for i = 1:length(t)
        if t(i) == 0
            srrcPulse(i) = 1 - alpha + 4*alpha/pi;
        elseif abs(t(i)) == Ts/(4*alpha)
            srrcPulse(i) = (alpha/sqrt(2)) * ((1+2/pi)*sin(pi/(4*alpha)) + (1-2/pi)*cos(pi/(4*alpha)));
        else
            num = sin(pi*t(i)/Ts*(1-alpha)) + 4*alpha*t(i)/Ts.*cos(pi*t(i)/Ts*(1+alpha));
            den = pi*t(i)/Ts .* (1 - (4*alpha*t(i)/Ts).^2);
            srrcPulse(i) = num / den;
        end
    end
    srrcPulse = srrcPulse / sqrt(sum(srrcPulse.^2));
end

function plotConstellation(symbolMapping, symbolRate, samplesPerSymbol, srrcPulse, Eb_N0_dB, fixedPrecision)
    numSymbols = 2048;
    wordLength = 12;         % Total number of bits
    fractionLength = 7;      % Number of fractional bits
    symbols = symbolMapping(randi(4, 1, numSymbols));
    
    % Upsample and apply pulse shaping
    upsampled = upsample(symbols, samplesPerSymbol);
    if fixedPrecision == true
        upsampled = fi(upsampled, 1, wordLength, fractionLength, 'SumMode', 'SpecifyPrecision', 'SumWordLength', wordLength, 'SumFractionLength', fractionLength, 'ProductMode', 'SpecifyPrecision', 'ProductWordLength', wordLength, 'ProductFractionLength', fractionLength);
    end
    transmittedSignal = conv(upsampled, srrcPulse, 'same');
    transmittedSignal2 = conv(transmittedSignal, srrcPulse, 'same');
    
    % Generate and shape noise
    Eb_N0 = 10^(Eb_N0_dB/10);
    noisePower = 1 / (2 * Eb_N0);
    if fixedPrecision == true
        noisePower = fi(noisePower, 1, wordLength, fractionLength, 'SumMode', 'SpecifyPrecision', 'SumWordLength', wordLength, 'SumFractionLength', fractionLength, 'ProductMode', 'SpecifyPrecision', 'ProductWordLength', wordLength, 'ProductFractionLength', fractionLength);
    end

    noise = (randn(size(transmittedSignal)) + 1i*randn(size(transmittedSignal)));
    if fixedPrecision == true
        noise = fi(noise, 1, wordLength, fractionLength, 'SumMode', 'SpecifyPrecision', 'SumWordLength', wordLength, 'SumFractionLength', fractionLength, 'ProductMode', 'SpecifyPrecision', 'ProductWordLength', wordLength, 'ProductFractionLength', fractionLength);
    end

    shapedNoise = conv(noise, srrcPulse, 'same');
    % Pn = sum(abs(shapedNoise).^2)/length(shapedNoise);
    shapedNoise = shapedNoise*sqrt(noisePower);
    
    % Add shaped noise to signal
    receivedSignal = transmittedSignal2 + shapedNoise;
    
    % receivedSignal = transmittedSignal;
    
    % Matched filter and downsample
    % matchedFilterOutput = conv(receivedSignal, srrcPulse, 'same');
    downsampled = receivedSignal(1:samplesPerSymbol:end);
    
    figure;
    scatter(real(downsampled), imag(downsampled), '.');
    hold on;
    scatter(real(symbolMapping), imag(symbolMapping), 'ro', 'filled');
    title(['QPSK Constellation (Eb/N0 = ' num2str(Eb_N0_dB) ' dB)']);
    xlabel('In-phase');
    ylabel('Quadrature');
    grid on;
    axis equal;
end

function [ber, ser, EbN0] = simulateBERSER_preamble(prb,samplingOffset, symbolMapping, symbolRate, samplesPerSymbol, srrcPulse, Eb_N0_dB, numBursts, numSymbols, fixedPrecision)
    ber = zeros(size(Eb_N0_dB));
    ser = zeros(size(Eb_N0_dB));
    EbN0 = 10.^(Eb_N0_dB/10);

    bpsk_mod_prb = 2*prb-1;
    % Length of the preamble
    upsampled_prb = upsample(bpsk_mod_prb, samplesPerSymbol);
    matched_prb = conv(conv(upsampled_prb, srrcPulse, 'same'), srrcPulse, 'same');
    matched_prb_len = length(matched_prb);

    wordLength = 12;        
    fractionLength = 7;

    syncFrequency = 20;
    packetsForDetection = 5;
    numPackets = 100;
    bitsPerSymbol = 2;

    for i = 1:length(Eb_N0_dB)
        totalBitErrors = 0;
        totalSymbolErrors = 0;

        % Preallocate an array to hold the QPSK-modulated symbols and bits
        packets = zeros(1, (length(bpsk_mod_prb)+numSymbols)*numPackets);
        bits = zeros(1, bitsPerSymbol*numSymbols*numPackets);
        % Map symbols to QPSK constellation points using Gray coding
        for packet = 1:numPackets
            pktsymbols = randi([1, 4], 1, numSymbols);
            packets((packet-1)*(length(bpsk_mod_prb)+numSymbols) + 1 : packet*(length(bpsk_mod_prb)+numSymbols)) = [bpsk_mod_prb symbolMapping(pktsymbols)];
            pktbits = de2bi(pktsymbols - 1, 2, 'left-msb')';
            pktbits = pktbits(:)';
            bits((packet-1)*numSymbols*bitsPerSymbol + 1 : packet*numSymbols*bitsPerSymbol) = pktbits;
        end
        % packets(1:1050) = randn(1,1050);
        % Upsample and apply pulse shaping
        txSig = upsample(packets, samplesPerSymbol);
        if fixedPrecision == true
            upsampled = fi(upsampled, 1, wordLength, fractionLength, 'SumMode', 'SpecifyPrecision', 'SumWordLength', wordLength, 'SumFractionLength', fractionLength, 'ProductMode', 'SpecifyPrecision', 'ProductWordLength', wordLength, 'ProductFractionLength', fractionLength);
        end

        txSig = conv(txSig, srrcPulse, 'same');
        
        % Generate and shape noise
        meanval = mean(txSig.*conj(txSig))/2;
        noisePower =  (samplesPerSymbol*meanval) / (2 * EbN0(i));
        if fixedPrecision == true
            noisePower = fi(noisePower, 1, wordLength, fractionLength, 'SumMode', 'SpecifyPrecision', 'SumWordLength', wordLength, 'SumFractionLength', fractionLength, 'ProductMode', 'SpecifyPrecision', 'ProductWordLength', wordLength, 'ProductFractionLength', fractionLength);
        end
        noise = (randn(size(txSig)) + 1i*randn(size(txSig)));
        if fixedPrecision == true
            noise = fi(noise, 1, wordLength, fractionLength, 'SumMode', 'SpecifyPrecision', 'SumWordLength', wordLength, 'SumFractionLength', fractionLength, 'ProductMode', 'SpecifyPrecision', 'ProductWordLength', wordLength, 'ProductFractionLength', fractionLength);
        end
        noise = noise*sqrt(noisePower);

        rxSig = conv(txSig + noise, srrcPulse, 'same');
        % fprintf('length of rxSig: %i\n', length(rxSig))

        pkt_len = (length(bpsk_mod_prb)+numSymbols);
        current_pkt_pointer = 1;
        current_bits_pointer = 1;
        data_pkts_seen = 0;
        start_found = false;
        while current_pkt_pointer < length(rxSig)
            
            if start_found == false 
                currentPkt = rxSig(1, current_pkt_pointer : current_pkt_pointer + samplesPerSymbol*pkt_len - 1);
                % Initialize the correlation array
                correlation = zeros(length(currentPkt) - matched_prb_len + 1, 1);
                % Loop to calculate the correlation at each point
                for m = 1:(length(currentPkt) - matched_prb_len + 1)
                    sumProduct = 0;
                    for n = 1:matched_prb_len
                        sumProduct = sumProduct + real(currentPkt(m + n - 1)) * matched_prb(n);
                    end
                    correlation(m) = abs(sumProduct);
                end

                [maxcorr, maxindex] = max(correlation);
                if maxcorr > 130 % start of burst found
                    % sampleIdx = maxindex + 26 * samplesPerSymbol + samplingOffset;
                    start_found = true;
                    fprintf('start of burst found at index %i\n', maxindex);
                else
                    current_pkt_pointer = current_pkt_pointer + samplesPerSymbol*pkt_len;
                end

            else
                % Initialize the correlation array 
                correlation = zeros(length(currentPkt) - matched_prb_len + 1, 1);
                avging_pkt_pointer = current_pkt_pointer;
                for x = 1: packetsForDetection
                    currentPkt = rxSig(1, avging_pkt_pointer : avging_pkt_pointer + samplesPerSymbol*pkt_len - 1);
                    % Loop to calculate the correlation at each point
                    for m = 1:(length(currentPkt) - matched_prb_len + 1)
                        sumProduct = 0;
                        for n = 1:matched_prb_len
                            sumProduct = sumProduct + real(currentPkt(m + n - 1)) * matched_prb(n);
                        end
                        correlation(m) = correlation(m) + (abs(sumProduct)/packetsForDetection);
                    end
                    avging_pkt_pointer = avging_pkt_pointer + samplesPerSymbol*pkt_len;
                end

                [maxcorr, maxindex] = max(correlation);
                % fprintf('maxcorr: %f\n', maxcorr)
                if maxcorr > 130 % start of burst found
                    sampleIdx = maxindex + 26 * samplesPerSymbol + samplingOffset;

                    for y = 1:syncFrequency
                        downsampled = rxSig(current_pkt_pointer + sampleIdx - 1:samplesPerSymbol:current_pkt_pointer + sampleIdx - 1 + (numSymbols - 1) * samplesPerSymbol);
                        % Detect bits based on the sign of the real and imaginary parts
                        detectedBitsReal = real(downsampled) < 0; % 1 if real part is negative, 0 otherwise
                        detectedBitsImag = imag(downsampled) < 0; % 1 if imaginary part is negative, 0 otherwise
                        % Combine the real and imaginary parts to get the detected bits
                        detectedBits = [detectedBitsReal; detectedBitsImag];
                        % Reshape the detected bits to match the original format
                        detectedBits = detectedBits(:)';
                        trueBits = bits(current_bits_pointer:current_bits_pointer + numSymbols*bitsPerSymbol - 1);
                        % Count errors
                        totalBitErrors = totalBitErrors + sum(trueBits ~= detectedBits);
                        current_pkt_pointer = current_pkt_pointer + samplesPerSymbol*pkt_len;
                        % fprintf('current_pkt_pointer: %i\n', current_pkt_pointer)
                        current_bits_pointer = current_bits_pointer + numSymbols*bitsPerSymbol;
                        % fprintf('current_bits_pointer: %i\n', current_bits_pointer)
                        data_pkts_seen = data_pkts_seen + 1;
                    end
                else
                    error('Correlation avging failed');
                end
            end
        end
    % totalBitErrors
    % totalSymbolErrors
    ber(i) = totalBitErrors / (2 * numSymbols * numPackets);
    % ser(i) = totalSymbolErrors / (numSymbols * numBursts);
    end
end

function plotBERComparison(EbN0, ber, offset)
    figure;
    semilogy(10*log10(EbN0), ber, 'bo-', 'LineWidth', 2);
    hold on;
    EbN0_theory = 10.^((-1:0.1:10)/10);
    ber_theory = qfunc(sqrt(2*EbN0_theory));
    semilogy(10*log10(EbN0_theory), ber_theory, 'r-', 'LineWidth', 2);
    grid on;
    xlabel('Eb/N0 (dB)');
    ylabel('Bit Error Rate (BER)');
    title('BER vs Eb/N0 for QPSK, Offset = ', offset);
    legend('Simulation', 'Theoretical');
end

function plotSER(EsN0_dB, ser)
    figure;
    semilogy(EsN0_dB, ser, 'bo-', 'LineWidth', 2);
    grid on;
    xlabel('Es/N0 (dB)');
    ylabel('Symbol Error Rate (SER)');
    title('SER vs Es/N0 for QPSK');
end