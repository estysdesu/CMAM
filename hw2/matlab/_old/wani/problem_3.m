% % Computational Methods in Additive Manufacturing - Assignment 2
% Tyler Estes (estests@mail.uc.edu)
% Utkarsh Wani (wanius@mail.uc.edu)

%% Problem 3 

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


P3= [2 8 4];
M3= [ -1 1 -1];
P4= [-2 5 4];
M4= [1 2 -1];

B1= [P3; P4; M3; M4];

% plotting
u = 0:0.0125:1
for n = 1:length(u)
    U_vector = [u(n)^3, u(n)^2, u(n), 1];
    P1(n, :) = U_vector*MsubH*B1;
end

figure()
plot3(P1(:, 1), P1(:, 2), P1(:, 3))
title('Problem 3a - Hermite Cubic Curve')

figure()
plot3(P(:, 1), P(:, 2), P(:, 3))
hold on
plot3(P1(:, 1), P1(:, 2), P1(:, 3))
title('Problem 3a - Hermite Cubic Curve combined plot')

