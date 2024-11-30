function agc_signal = automaticGainControl(signal)

agc_signal = signal/max(abs(signal));
