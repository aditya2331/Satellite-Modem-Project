function rx_matchedfilter_output = matchedFilter(fo_corrected_signal, MF_FIR_coefficients)
%This function performs matched filtering on the input signal

rx_matchedfilter_output = conv(fo_corrected_signal, MF_FIR_coefficients);

end