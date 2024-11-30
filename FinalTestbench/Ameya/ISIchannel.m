function [faded_signal, cir] = ISIchannel(pulseShaped_signal, upsampling_factor, no_paths)

    bandwidth = 4.8e3;
    delay_spread = 1e-4;
    L = no_paths - 1;

    cir = exp(-(0:L)'/(bandwidth*delay_spread)).*(randn(L+1, 1) + 1i*randn(L+1, 1));
    cir = cir/norm(cir); %[1 0 0 0 0];%

    upsampled_cir = zeros((L+1)*upsampling_factor, 1);
    upsampled_cir(1:upsampling_factor:end) = cir;

    faded_signal = conv(upsampled_cir, pulseShaped_signal);
