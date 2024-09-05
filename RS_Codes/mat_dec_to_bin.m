function out = mat_dec_to_bin(inp, m)
    out = zeros(size(inp, 1), size(inp, 2), m);

    for id = 1:numel(inp)
        binStr = dec2bin(inp(id), m);                  % Convert to binary string
        binVec = arrayfun(@(x) str2double(x), binStr); % Convert string to numeric array
        out(id + (0:m-1)*numel(inp)) = binVec;         % Place the binary digits in the 3rd dimension of the output matrix
    end
end