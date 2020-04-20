%% Problem 4

%% Part A

B= [ 1 1 1 ;1 3 3; 3 2 1; 2 4 2]; % Control Points have been represented in B Matrix
M_H=[-1 3 -3 1; 3 -6 3 0; -3 3 0 0; 1 0 0 0]; % Hermite Matrix for the 3rd Degree curve i.e with 4 control points

u= 0:0.025:1; % value of u varies from 0-1 with an interval of 0.025 
u_vector= [u.^3; u.^2; u; ones(1, length(u))]';

P= u_vector*M_H*B;
figure()
plot3(P(:, 1), P(:, 2), P(:, 3))

%% Part B

B1= [ 1 1 1 ;1 3 3; 3 2 1; 3 2 1; 2 4 2];% New Control Points have been represented in B1 Matrix
M_H_1= [1 -4 6 -4 1; -4 12 -12 4 0; 6 -12 6 0 0; -4 4 0 0 0; 1 0 0 0 0]; % Hermite Matrix for the 4th Degree curve i.e with 5 control points

u1= 0:0.005:1; % value of u1 varies from 0-1 with an interval of 0.005
u_vector_1= [ u1.^4; u1.^3; u1.^2; u1; ones(1, length(u1))]';

P1= u_vector_1*M_H_1*B1

figure()
plot3(P1(:, 1), P1(:, 2), P1(:, 3))

% tan vector & unit tan vector
%% u=0
u = 0;
uVec_1 = [ 4*(u^3) 3*(u^2) 2*u 1 0]
PsupU_1 = uVec_1*M_H_1*B1;
PsupUNorm_1 = PsupU_1./norm(PsupU_1);
fprintf('Problem 3a - Tangent Vector = [%.2f, %.2f, %.2f], Unit Tangent Vector = [%.2f, %.2f, %.2f]\n', PsupU_1(1), PsupU_1(2), PsupU_1(3), PsupUNorm_1(1), PsupUNorm_1(2), PsupUNorm_1(3))

%% u=0.5
u = 0.5;
uVec_2 = [ 4*(u^3) 3*(u^2) 2*u 1 0]
PsupU_2 = uVec_2*M_H_1*B1;
PsupUNorm_2 = PsupU_2./norm(PsupU_2);
fprintf('Problem 3a - Tangent Vector = [%.2f, %.2f, %.2f], Unit Tangent Vector = [%.2f, %.2f, %.2f]\n', PsupU_2(1), PsupU_2(2), PsupU_2(3), PsupUNorm_2(1), PsupUNorm_2(2), PsupUNorm_2(3))

%% u=0.75
u = 0.75;
uVec_3 = [ 4*(u^3) 3*(u^2) 2*u 1 0]
PsupU_3 = uVec_3*M_H_1*B1;
PsupUNorm_3 = PsupU_3./norm(PsupU_3);
fprintf('Problem 3a - Tangent Vector = [%.2f, %.2f, %.2f], Unit Tangent Vector = [%.2f, %.2f, %.2f]\n', PsupU_3(1), PsupU_3(2), PsupU_3(3), PsupUNorm_3(1), PsupUNorm_3(2), PsupUNorm_3(3))

%% u=1
u = 1;
uVec_4 = [ 4*(u^3) 3*(u^2) 2*u 1 0]
PsupU_4 = uVec_4*M_H_1*B1;
PsupUNorm_4 = PsupU_4./norm(PsupU_4);
fprintf('Problem 3a - Tangent Vector = [%.2f, %.2f, %.2f], Unit Tangent Vector = [%.2f, %.2f, %.2f]\n', PsupU_4(1), PsupU_4(2), PsupU_4(3), PsupUNorm_4(1), PsupUNorm_4(2), PsupUNorm_4(3))


%% Part C

B2= [ 1 1 1 ;0 5 -2; 3 2 1; 3 2 1; 2 4 2]; % New Control Points have been represented in B2 Matrix
M_H_2=[1 -4 6 -4 1; -4 12 -12 4 0; 6 -12 6 0 0; -4 4 0 0 0; 1 0 0 0 0]; % Hermite Matrix for the 4th Degree curve i.e with 5 control points

u2= 0:0.0025:1; % value of u2 varies from 0-1 with an interval of 0.0025
u_vector_2= [ u2.^4; u2.^3; u2.^2; u2; ones(1, length(u2))]';

P2= u_vector_2*M_H_2*B2

figure()
plot3(P2(:, 1), P2(:, 2), P2(:, 3))

figure()
plot3(P1(:, 1), P1(:, 2), P1(:, 3))
hold on
plot3(P2(:, 1), P2(:, 2), P2(:, 3))


