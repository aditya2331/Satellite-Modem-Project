q = [3,4,5,6,7,8];
prim = {[1 0 1 1], [1 0 0 1 1], [1 0 0 1 0 1], [1 0 0 0 0 1 1], [1 0 0 0 1 0 0 1], [1 0 0 0 1 1 1 0 1]}; 
% Use Cell Arrays to store arrays of different sizes, types or dimensions

[AT3, MT3] = table_generator(q(1), prim{1});
[AT4, MT4] = table_generator(q(2), prim{2});
[AT5, MT5] = table_generator(q(3), prim{3});
[AT6, MT6] = table_generator(q(4), prim{4});
[AT7, MT7] = table_generator(q(5), prim{5});
[AT8, MT8] = table_generator(q(6), prim{6});

save("Tables.mat",'AT3','MT3','AT4','MT4','AT5','MT5','AT6','MT6','AT7','MT7','AT8','MT8');
