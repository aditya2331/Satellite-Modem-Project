% Assumes that inp, l, q are decimal numbers 
function out = gf_ele_exponent(inp, l, q)
    data = load("Tables.mat");
    mul_table = data.(sprintf('MT%d',q));
    
    out = 1;
    for i=1:l
        temp = out;
        out = mul_table(inp+1,temp+1);
    end
end

