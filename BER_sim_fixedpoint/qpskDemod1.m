function out_sig = qpskDemod1(in_sig)

out_sig = zeros(size(in_sig*2));

for i = 1:1:length(in_sig)
    k = in_sig(i);
    if real(k) > 0 && imag(k) > 0
        out_sig(2*i-1) = 0;
        out_sig(2*i) = 0;
    elseif real(k) <= 0 && imag(k) > 0 
        out_sig(2*i-1) = 0;
        out_sig(2*i) = 1;
    elseif real(k) <= 0 && imag(k) <= 0
        out_sig(2*i-1) = 1;
        out_sig(2*i) = 1;
    else
        out_sig(2*i-1) = 1;
        out_sig(2*i) = 0;
    end
end

end