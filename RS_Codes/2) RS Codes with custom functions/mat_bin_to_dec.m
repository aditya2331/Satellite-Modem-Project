function out = mat_bin_to_dec(inp)
    out = zeros(size(inp,1), size(inp,2));
    
    % Binary number is assumed to be stored in Big Endian Format, ie, [1 1 0] = x^2 + x
    for i=1:size(inp,1)
        for j=1:size(inp,2)
            bin_num = inp(i,j,:);
            bin_num = flip(bin_num(:)); % Flattening it to a row vector  and storing in Little Endian Format
            dec_num = 0;
            for k=1:length(bin_num)
                dec_num = dec_num + bin_num(k)*(2^(k-1));
            end
            out(i,j) = dec_num;
        end
    end
end