clear all, close all, clc
%% Q1
clear all

Bsub0 = [1, 2];
Bsub1 = [2, 7];
Bsub2 = [5, 12];
Bsub3 = [3, 2];

% plotting
u = 0:0.0125:1;
Psub0 = Bsub0'*(1-u.^3-3*u+3*u.^2);
Psub1 = Bsub1'*(3*u+3*u.^3-6*u.^2);
Psub2 = Bsub2'*(3*u.^2-3*u.^3);
Psub3 = Bsub3'*(u.^3);
P = (Psub0 + Psub1 + Psub2 + Psub3)'; % composed P with individual terms

figure
plot(P(:, 1), P(:, 2))
title('4th Order Open Non-Periodic B-Spline Curve')
grid on

% solve at specfic u values
u = [.3, .6];
Psub0 = Bsub0'*(1-u.^3-3*u+3*u.^2);
Psub1 = Bsub1'*(3*u+3*u.^3-6*u.^2);
Psub2 = Bsub2'*(3*u.^2-3*u.^3);
Psub3 = Bsub3'*(u.^3);
P = (Psub0 + Psub1 + Psub2 + Psub3)'; % composed P with individual terms

disp('Problem 1')
fprintf('Point @ (u = %.2f) = [%.3f, %.3f]\n', u(1), P(1, 1), P(1, 2))
fprintf('Point @ (u = %.2f) = [%.3f, %.3f]\n', u(2), P(2, 1), P(2, 2))
%% Q2
clear all

P0 = [-1, 2];
P1 = [1.75, 4];
P2 = [2, 1];
P3 = [2.25, 4];
P4 = [5, 2];
P5 = [2, -1];
P6 = [2, -1];

M4 = 1/6*[
    -1 3 -3 1;
    3 -6 3 0;
    -3 0 3 0;
    1 4 1 0]; % 4x4

u=0:0.0125:1;
uVec = [u.^3; u.^2; u; ones(1, length(u))]';
P_1 = uVec*M4*[P0; P1; P2; P3];
P_2 = uVec*M4*[P1; P2; P3; P4];
P_3 = uVec*M4*[P2; P3; P4; P5];
P_4 = uVec*M4*[P3; P4; P5; P6];
P_5 = uVec*M4*[P4; P5; P6; P0];
P_6 = uVec*M4*[P5; P6; P0; P1];
P_7 = uVec*M4*[P6; P0; P1; P2];

figure(), hold on, grid on
plot(P_1(:, 1), P_1(:, 2))
plot(P_2(:, 1), P_2(:, 2))
plot(P_3(:, 1), P_3(:, 2))
plot(P_4(:, 1), P_4(:, 2))
plot(P_5(:, 1), P_5(:, 2))
plot(P_6(:, 1), P_6(:, 2))
plot(P_7(:, 1), P_7(:, 2))
title('4th order Closed B-Spline Curve')
%% Q3
clear all

Bsub0 = [0, 0];
Bsub1 = [0, 2];
Bsub2 = [2, 2];
Bsub3 = [4, 2];
Bsub4 = [8, 5];
Bsub5 = [3, 7];

u = 0:0.0125:3;
P = zeros(length(u), 2);

for i = 1:length(u)
    if u(i) >=0 && u(i) <= 1
        P(i, :) = Bsub0*(1-u(i))^2 + Bsub1*(2*u(i)*(1-u(i))) + Bsub2*u(i)^2;
    elseif u(i) >= 1 && u(i) <= 2
        P(i, :) = Bsub2*(2-u(i))^2 + Bsub3*( (2-u(i))*(u(i)-1)+.5*(3-u(i))*(u(i)-1) ) + Bsub4*(.5*(u(i)-1)^2);
    elseif u(i) >= 2 && u(i) <= 3
        P(i, :) = Bsub3*(.5*(3-u(i))^2) + Bsub4*(.5*(u(i)-1)*(3-u(i))+(3-u(i))*(u(i)-2)) + Bsub5*(u(i)-2)^2;
    else
        P(i, :) = 0;
    end
end

figure
plot(P(:, 1), P(:, 2))
title('Composite Curve Plot')
grid on


u = [1.5, 2];
P = zeros(length(u), 2);
for i = 1:length(u)
    if u(i) >=0 && u(i) <= 1
        P(i, :) = Bsub0*(1-u(i))^2 + Bsub1*(2*u(i)*(1-u(i))) + Bsub2*u(i)^2;
    elseif u(i) >= 1 && u(i) <= 2
        P(i, :) = Bsub2*(2-u(i))^2 + Bsub3*( (2-u(i))*(u(i)-1)+.5*(3-u(i))*(u(i)-1) ) + Bsub4*(.5*(u(i)-1)^2);
    elseif u(i) >= 2 && u(i) <= 3
        P(i, :) = Bsub3*(.5*(3-u(i))^2) + Bsub4*(.5*(u(i)-1)*(3-u(i))+(3-u(i))*(u(i)-2)) + Bsub5*(u(i)-2)^2;
    else
        P(i, :) = 0;
    end
end

disp('Problem 3')
fprintf('Point @ (u = %.2f) = [%.3f, %.3f]\n', u(1), P(1, 1), P(1, 2))
fprintf('Point @ (u = %.2f) = [%.3f, %.3f]\n', u(2), P(2, 1), P(2, 2))