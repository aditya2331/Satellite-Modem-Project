function [noisy_signal, noise_sd] = addNoise(signal, EbN0dB, upsampling_factor)

    EbN0 = 10^(EbN0dB/10);
    Eb = mean(double(signal).*conj(double(signal)))/2;
    noise_var = upsampling_factor*Eb/EbN0;
    noise_sd = sqrt(noise_var/2);
    noise = noise_sd*(randn(numel(signal),1) + 1i*randn(numel(signal),1));
    
    noisy_signal = signal + noise;  