A = [1 0 0 6 1 6 7;
     0 1 0 4 1 5 5;
     0 0 1 3 1 2 3];

A = A.';

B = [6 2 3; 7 4 5; 1 7 1];

prim = [1 0 1 1];
q = 3;

mul = galois_table_mat_mul(A,B,q);
mul_2 = galois_mat_mul(A,B,prim,q);
disp(mul);
disp(mul_2);