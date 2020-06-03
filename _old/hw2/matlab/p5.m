% #####
% Homework 2 Problem 5
% Tyler Estes
% Utkarsh Wani
% #####

%% p5
C = [0, 0, 0, 0]; % center
P0 = [-5, 0, 0, 0]; % left endpoint
P3 = [0, 5, 0, 0]; % right endpoint

k = 2.7613; % solved for manually
P1 = [-5, k, 0, 0];
P2 = [-k, 5, 0, 0];

B = [P0; P1; P2; P3];
MsubB = [-1, 3, -3, 1;
    3, -6, 3, 0;
    -3, 3, 0, 0;
    1, 0, 0, 0];
uVals = 0:.025:1;
U = [uVals.^3; uVals.^2; uVals; ones(1, length(uVals))]';
P = U*MsubB*B;

% plotting
plot(P(:, 1), P(:, 2))
title('Problem 5 - Quarter Circle Approximation w/ Cubic Bezier Curve')

u = [0.25, 0.75];
radErr = sqrt(5.89*u.^6 - 17.67*u.^5 + 19.13*u.^4 - 8.8140*u.^3 + 1.46*u.^2 + 25) - 5;
fprintf('Problem 5 - Radial Error (@ u = %.2f) = %.5f\n', u(1), radErr(1))
fprintf('Problem 5 - Radial Error (@ u = %.2f) = %.5f\n', u(2), radErr(2))

