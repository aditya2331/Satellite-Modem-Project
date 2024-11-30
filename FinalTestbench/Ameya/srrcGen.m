function out = srrcGen(alpha,truncation_length,oversampling_factor,type,wordLength,fractionLength)

Ts = 1;                                   % Symbol period
T = Ts/oversampling_factor;               % Sampling period
L = truncation_length/2;
t = -L*Ts : T : L*Ts;
out = zeros(size(t));
for i = 1:1:length(t)
    if t(i) == 0
        out(i) = 1 - alpha + (4 * alpha / pi);
    elseif abs(t(i)) == Ts/(4*alpha)
        out(i) = (alpha/sqrt(2))*((1 + 2/pi)*sin(pi/(4*alpha)) + (1 - 2/pi)*cos(pi/(4*alpha)));
    else
        out(i) = (sin(pi*(1-alpha)*t(i)/Ts) + (4*alpha*t(i)/Ts)*cos(pi*(1+alpha)*t(i)/Ts)) / ((pi*t(i)/Ts)*(1 - (4*alpha*t(i)/Ts)^2));
    end
end

if strcmp(type, 'normalized')
    out = out / sqrt(sum(out.^2)); % Normalize the output
elseif strcmp(type, 'original')
    % Keep output as it is
    out = out;
else
    error('Invalid type. Use ''normalized'' or ''original''.');
end

out = fi(out, 1, wordLength, fractionLength, 'SumMode', 'SpecifyPrecision', 'SumWordLength', wordLength, 'SumFractionLength', fractionLength, 'ProductMode', 'SpecifyPrecision', 'ProductWordLength', wordLength, 'ProductFractionLength', fractionLength);

end