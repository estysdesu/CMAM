%% Computational Methods in Additive Manufacturing - Assignment 2
% Tyler Estes (estests@mail.uc.edu)
% Utkarsh Wani (wanius@mail.uc.edu)

%% Problem 3 
clc, clear all, close all

P1 = [4, 2, 6];
M1 = [3, 1, -1];
P2 = [2, 8, 4];
M2 = [-1, 1, -1];

MsubH = [2, -2, 1, 1;
    -3, 3, -2, -1;
    0, 0, 1, 0;
    1, 0, 0, 0];
MsubHsupU = [0, 0, 0, 0;
    6, -6, 3, 3;
    -6, 6, -4, -2;
    0, 0, 1, 0];
B = [P1; P2; M1; M2];

% plotting
uStart = 0;
uStop = 1;
step = .05;
n = uStop/step;
P = zeros(n, 3);
u = linspace(uStart, uStop, n);
for n = 1:length(u)
    uVec = [u(n)^3, u(n)^2, u(n), 1];
    P(n, :) = uVec*MsubH*B;
end
plot3(P(:, 1), P(:, 2), P(:, 3))
title('Problem 3a - Hermite Cubic Curve')

% tan vector & unit tan vector
u = .6;
uVec = [u^3, u^2, u, 1];
PsupU = uVec*MsubHsupU*B;
PsupUNorm = PsupU./norm(PsupU);
fprintf('Problem 3a - Tangent Vector = [%.2f, %.2f, %.2f], Unit Tangent Vector = [%.2f, %.2f, %.2f]\n', PsupU(1), PsupU(2), PsupU(3), PsupUNorm(1), PsupUNorm(2), PsupUNorm(3))

%% Problem 4

clc, clear all, close all
B= [ 1 1 1 ;1 3 3; 3 2 1; 2 4 2];
M_H=[-1 3 -3 1; 3 -6 3 0; -3 3 0 0; 1 0 0 0];
% u= [0 0.5 0.75 1]
% U= [u^3 u^2 u 1];

% P(u)= U*M_H*B;

u= 0:0.025:1;
uvec= [u.^3; u.^2; u; ones(1, length(u))]';

 P= uvec*M_H*B;
figure()
plot3(P(:, 1), P(:, 2), P(:, 3))

B1= [ 1 1 1 ;1 3 3; 3 2 1; 3 2 1; 2 4 2];
M_H=[-1 3 -3 1; 3 -6 3 0; -3 3 0 0; 1 0 0 0];
Msub_H= [1 -4 6 -4 1; -4 12 -12 4 0; 6 -12 6 0 0; -4 4 0 0 0; 1 0 0 0 0];
uvec1= [ u.^4; u.^3; u.^2; u; ones(1, length(u))]';

P1= uvec1*Msub_H*B1
figure()
plot3(P1(:, 1), P1(:, 2), P1(:, 3))

B2= [ 1 1 1 ;0 5 -2; 3 2 1; 3 2 1; 2 4 2];


%% Problem 5

clc, clear all, close all
v = [4, 6, -5; %
    14, 6, -5;
    4, 16, -5;
    4, 6, 5; %
    14, 16, 5;
    4, 16, 5;
    4, 6, 15; %
    14, 6, 15;
    14, 16, 15];

% face 7, 6, 3, 1, 4 -> face 7, 6, 3 + face 7, 3, 1 + face 7, 1, 4
% face 3, 6, 9, 5 -> face 3, 6, 9 + 3, 9, 5
f = [1, 3, 2;
    1, 2, 4;
    2, 3, 5;
    7, 6, 3;
    7, 3, 1;
    7, 1, 4;
    2, 5, 4;
    3, 6, 9;
    3, 9, 5;
    4, 5, 8;
    4, 8, 7;
    5, 9, 8;
    6, 7, 9;
    7, 8, 9];

figure()
h = patch('vertices', v, 'faces', f, 'facecolor', 'c')
grid; box;
view(3)


% total_vol = 0;
% for n = 1:length(f)
%     vol = volumeFromNormalAndFace(v( f(n, :), :))
%     total_vol = total_vol + vol;
% end
% total_vol = total_vol*1/6

v(4, :) - v(7, :)
v1 = v(1, :)
v4 = v(4, :)
v7 = v(7, :)
volumeFromNormalAndFace([v1; v4; v7])

