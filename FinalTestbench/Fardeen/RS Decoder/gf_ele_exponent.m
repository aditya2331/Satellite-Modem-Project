% Assumes that inp, l, q are decimal numbers 
function out = gf_ele_exponent(inp, l, q)
    persistent mul_table_cache q_cache;
    
    % Load the multiplication table only if it hasn't been loaded or if the 'q' value has changed
    if isempty(mul_table_cache) || q_cache ~= q
        data = load("Tables.mat");
        mul_table_cache = data.(sprintf('MT%d', q));
        q_cache = q;  % Cache the current 'q' value
    end
    
    mul_table = mul_table_cache;
    
    out = 1;
    for i=1:l
        temp = out;
        out = mul_table(inp+1,temp+1);
    end
end

