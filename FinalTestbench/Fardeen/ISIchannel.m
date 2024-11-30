function fadedSignal = ISIchannel(pulseShaped_signal, ISI_FIR_coeffecients)
    % pulseShaped_signal    :  Analog Signal that has been upsampled and convolved with SRRC at the Tx
    % ISI_FIR_coeffecients  :  FIR of the channel modelled by an exponentially decaying Normal Distribution 
    % FIR_coeffecients are assumed to be [h0, h1, ... ] where h0 is for no delay, h1 is for 1 time step delay and so on
    
    fadedSignal = zeros(1, size(pulseShaped_signal));                     % Defining output signal
    shift_reg = zeros(1, size(ISI_FIR_coeffecients));                     % Defining shift register to perform computation

    for i = 1:length(pulseShaped_signal)
        shift_reg = [pulseShaped_signal(i) shift_reg(1:end-1)];           % New element is added to the start of the shift register
        fadedSignal(i) = shift_reg.*ISI_FIR_coeffecients;
    end

    % fadedSignal = filter(ISI_FIR_coeffecients,1,pulseShaped_signal);
end

