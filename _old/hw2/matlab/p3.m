% #####
% Homework 2 Problem 3
% Tyler Estes
% Utkarsh Wani
% #####

%% p3a
% solve for corresponding cartesian points from u curve points
P1 = [4, 2, 6];
M1 = [3, 1, -1];
P2 = [2, 8, 4];
M2 = [-1, 1, -1];
B = [P1; P2; M1; M2];
MsubH = [2, -2, 1, 1;
    -3, 3, -2, -1;
    0, 0, 1, 0;
    1, 0, 0, 0];
uVals = 0:0.025:1;
U = [uVals.^3; uVals.^2; uVals; ones(1, length(uVals))]';
P = U*MsubH*B;

% plotting
figure()
lineA = plot3(P(:, 1), P(:, 2), P(:, 3));
grid on
title('Problem 3a - Hermite Cubic Curve')

% tangent vector at specific u value
u = .6;
UsupU = [3*u^2, 2*u, 1, 0]; % derivative of U vector
PsupU = UsupU*MsubH*B;
PsupUNorm = PsupU./norm(PsupU);
fprintf('Problem 3a - Tangent Vector (@ u = %.2f) = [%.2f, %.2f, %.2f], Unit Tangent Vector = [%.2f, %.2f, %.2f]\n', u, PsupU(1), PsupU(2), PsupU(3), PsupUNorm(1), PsupUNorm(2), PsupUNorm(3))

%% p3b
% solve for corresponding cartesian points from u curve points
P3 = P2;
M3 = M2;
P4 = [-2, 5, 4];
M4 = [1, 2, -1];
B = [P3; P4; M3; M4];
P = U*MsubH*B;

% plotting
figure()
plot3(P(:, 1), P(:, 2), P(:, 3))
grid on
title('Problem 3b - Hermite Cubic Curve w/ C1 Continuity from 3a')

figure()
plot3(get(lineA, 'xdata'), get(lineA, 'ydata'), get(lineA, 'zdata')), hold on % reusing plot data from part B
plot3(P(:, 1), P(:, 2), P(:, 3)), hold off
grid on
title('Problem 3b - Hermite Cubic Curve w/ C1 Continuity from 3a')