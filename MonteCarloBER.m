function qpskModemSimulation()
    % Parameters
    symbolRate = 4800; % 4.8 Ksymbols/sec
    samplesPerSymbol = 8; % Oversampling factor
    numSymbols = 2048;
    alpha = 0.4; % Roll-off factor for SRRC
    Eb_N0_dB = 0:2:12; % Eb/N0 range in dB
    numBursts = 500; % Number of bursts for averaging

    % QPSK symbol mapping (Gray coded)
    symbolMapping = [1+1i, -1+1i, 1-1i, -1-1i] / sqrt(2);

    % Generate SRRC pulse
    [srrcPulse, t] = generateSRRC(symbolRate, samplesPerSymbol, alpha);

    % (a) Constellation plot for Eb/N0 = 6 dB
    plotConstellation(symbolMapping, symbolRate, samplesPerSymbol, srrcPulse, 10);

    % (b) & (c) BER/SER simulation
    [ber, ser, EbN0] = simulateBERSER(symbolMapping, symbolRate, samplesPerSymbol, srrcPulse, Eb_N0_dB, numBursts, numSymbols);

    % (d) Plot BER vs Eb/N0 and compare with analytical
    plotBERComparison(EbN0, ber);

    % Plot SER vs Es/N0
    plotSER(EbN0 + 3, ser); % Es/N0 = Eb/N0 + 3dB for QPSK
end

function [srrcPulse, t] = generateSRRC(symbolRate, samplesPerSymbol, alpha)
    Ts = 1 / symbolRate;
    t = -5*Ts:Ts/samplesPerSymbol:5*Ts;
    
    % Create fixed-point objects with specific precision
    % Parameters
    wordLength = 12;         % Total number of bits
    fractionLength = 6;      % Number of fractional bits
    srrcPulse = fi(zeros(size(t)), 1, wordLength, fractionLength, 'SumMode', 'SpecifyPrecision', 'SumWordLength', wordLength, 'SumFractionLength', fractionLength, 'ProductMode', 'SpecifyPrecision', 'ProductWordLength', wordLength, 'ProductFractionLength', fractionLength);

    % srrcPulse = zeros(size(t));
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

function plotConstellation(symbolMapping, symbolRate, samplesPerSymbol, srrcPulse, Eb_N0_dB)
    numSymbols = 2048;
    symbols = symbolMapping(randi(4, 1, numSymbols));
    
    % Upsample and apply pulse shaping
    upsampled = upsample(symbols, samplesPerSymbol);
    transmittedSignal = conv(upsampled, srrcPulse, 'same');
    
    % Add noise
    Eb_N0 = 10^(Eb_N0_dB/10);
    noise = sqrt(1/(2*Eb_N0)) * (randn(size(transmittedSignal)) + 1i*randn(size(transmittedSignal)));
    % shapedNoise = conv(noise, srrcPulse, 'same');
            
    % Add shaped noise to signal
    receivedSignal = transmittedSignal + noise;
    % receivedSignal = transmittedSignal;
    
    % Matched filter and downsample
    matchedFilterOutput = conv(receivedSignal, srrcPulse, 'same');
    downsampled = matchedFilterOutput(1:samplesPerSymbol:end);
    
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

function [ber, ser, EbN0] = simulateBERSER(symbolMapping, symbolRate, samplesPerSymbol, srrcPulse, Eb_N0_dB, numBursts, numSymbols)
    ber = zeros(size(Eb_N0_dB));
    ser = zeros(size(Eb_N0_dB));
    EbN0 = 10.^(Eb_N0_dB/10);
    
    for i = 1:length(Eb_N0_dB)
        totalBitErrors = 0;
        totalSymbolErrors = 0;
        
        for burst = 1:numBursts
            % Generate random symbols
            symbolIndices = randi(4, 1, numSymbols);
            symbols = symbolMapping(symbolIndices);
            bits = de2bi(symbolIndices - 1, 2, 'left-msb')';
            bits = bits(:)';
            
            % Upsample and apply pulse shaping
            upsampled = upsample(symbols, samplesPerSymbol);
            transmittedSignal = conv(upsampled, srrcPulse, 'same');
            
            % Generate and shape noise
            noisePower = 1 / (2 * EbN0(i));
            noise = sqrt(noisePower) * (randn(size(transmittedSignal)) + 1i*randn(size(transmittedSignal)));
            % shapedNoise = conv(noise, srrcPulse, 'same');
            
            % Add shaped noise to signal
            receivedSignal = transmittedSignal + noise;
            % receivedSignal = transmittedSignal;
            
            % Matched filter
            matchedFilterOutput = conv(receivedSignal, srrcPulse, 'same');
            
            % Downsample
            downsampled = matchedFilterOutput(1:samplesPerSymbol:end);
            
            % Detect symbols
            [~, detectedSymbolIndices] = min(abs(downsampled.' - symbolMapping), [], 2);
            detectedSymbols = symbolMapping(detectedSymbolIndices);
            detectedBits = de2bi(detectedSymbolIndices - 1, 2, 'left-msb')';
            detectedBits = detectedBits(:)';
            
            % Count errors
            totalBitErrors = totalBitErrors + sum(bits ~= detectedBits);
            totalSymbolErrors = totalSymbolErrors + sum(symbols ~= detectedSymbols);
        end
        % totalBitErrors
        % totalSymbolErrors
        ber(i) = totalBitErrors / (2 * numSymbols * numBursts);
        ser(i) = totalSymbolErrors / (numSymbols * numBursts);
    end
end

function plotBERComparison(EbN0, ber)
    figure;
    % semilogy(10*log10(EbN0) - 3, ber, 'bo-', 'LineWidth', 2);
    hold on;
    EbN0_theory = 10.^((-1:0.1:10)/10);
    ber_theory = qfunc(sqrt(2*EbN0_theory));
    semilogy(10*log10(EbN0_theory), ber_theory, 'r-', 'LineWidth', 2);
    grid on;
    xlabel('Eb/N0 (dB)');
    ylabel('Bit Error Rate (BER)');
    title('BER vs Eb/N0 for QPSK');
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