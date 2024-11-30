function rx_decimated = polyphaseDecimator(rx_signal, decimation_factor, Dec_FIR_coefficients)
%This function does decimation by a factor of 25 on the received signal

output_len = ceil((length(rx_signal) + length(Dec_FIR_coefficients) - 1)/decimation_factor);
rx_decimated = zeros(output_len, 1);

for i = 1 : decimation_factor
    if i == 1
        inp_poly = rx_signal(1:decimation_factor:end);
        poly_filt = Dec_FIR_coefficients(i:decimation_factor:end);
        temp = conv(poly_filt,inp_poly);
        rx_decimated(1:length(temp)) = rx_decimated(1:length(temp)) + temp;
    else
        inp_poly = rx_signal(27-i:decimation_factor:end);
        poly_filt = Dec_FIR_coefficients(i:decimation_factor:end);
        temp = conv(poly_filt,inp_poly);
        rx_decimated(2:length(temp)+1) = rx_decimated(2:length(temp)+1) + temp;
    end

end

end