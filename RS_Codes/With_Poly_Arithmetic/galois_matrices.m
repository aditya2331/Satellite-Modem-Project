5% GF(3) Matrix
msg3 = gf(eye(3),3);
G3 = rsenc(msg3,7,3)';
G3 = G3.x;

% GF(4) Matrix
msg4 = gf(eye(7),4);
G4 = rsenc(msg4,15,7)';
G4 = G4.x;

% GF(5) Matrix
msg5 = gf(eye(19),5);
G5 = rsenc(msg5,31,19)';
G5 = G5.x;

% GF(6) Matrix
msg6 = gf(eye(47),6);
G6 = rsenc(msg6,63,47)';
G6 = G6.x;

% GF(7) Matrix
msg7 = gf(eye(103),7);
G7 = rsenc(msg7,127,103)';
G7 = G7.x;

% GF(8) Matrix
msg8 = gf(eye(223),8);
G8 = rsenc(msg8,255,223)';
G8 = G8.x;

save('galois_matrices.mat','G3','G4','G5','G6','G7','G8');