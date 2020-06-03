%% Problem 5

M_Bazier= [-1 3 -3 1; 3 -6 3 0; -3 3 0 0 ; 1 0 0 0]; % 3rd degree Bazier Matrices
B= [ -5 0; -5 2.7613; -2.7613 5; 0 5];

u= 0:0.025:1; % value of u varies from 0-1 with an interval of 0.025 
u_vector= [u.^3; u.^2; u; ones(1, length(u))]';

P= u_vector*M_Bazier*B;

plot(P(:,1), P(:,2))

