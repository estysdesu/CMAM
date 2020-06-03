% #####
% Homework 2 Problem 4
% Tyler Estes
% Utkarsh Wani
% #####

%% p4a
% solve for corresponding cartesian points from u curve points
pointA = [1, 1, 1];
pointB = [1, 3, 3];
pointC = [3, 2, 1];
pointE = [2, 4, 2];
B = [pointA; pointB; pointC; pointE];
MsubB = [-1, 3, -3, 1;
    3, -6, 3, 0;
    -3, 3, 0, 0;
    1, 0, 0, 0];
uVals = 0:0.025:1;
U = [uVals.^3; uVals.^2; uVals; ones(1, length(uVals))]';
P = U*MsubB*B;

% plotting
figure()
plot3(P(:, 1), P(:, 2), P(:, 3))
grid on
title('Problem 4a - Cubic Bezier Curve')

%% p4b
% solve for corresponding cartesian points from u curve points
pointD = pointC; % pointD is identical to pointC
B = [B(1:3, :); pointD; B(4:end, :)]; % insert pointD
MsubB= [1, -4, 6, -4, 1;
    -4, 12, -12, 4, 0;
    6, -12, 6, 0, 0;
    -4, 4, 0, 0, 0;
    1, 0, 0, 0, 0];
uVals = 0:0.025:1;
U = [uVals.^4; uVals.^3; uVals.^2; uVals; ones(1, length(uVals))]';
P = U*MsubB*B;

% plotting
figure()
lineB = plot3(P(:, 1), P(:, 2), P(:, 3));
grid on
title('Problem 4b - Fourth Order Bezier Curve')

% tangent vectors @ specific u values
for u = [0, 0.5, 0.75, 1]
    UsupU = [4*u^3, 3*u^2, 2*u, 1, 0]; % derivative of U vector
    V = UsupU*MsubB*B;
    fprintf('Problem 4b - Tangent Vector (@ u = %.2f) = [%.2f, %.2f, %.2f]\n', u, V(1), V(2), V(3))
end

%% p4c
% solve for corresponding cartesian points from u curve points
pointB = [0, 5, -2]; % update pointB
B(2, :) = pointB; % change B vector to new pointB
MsubB= [1, -4, 6, -4, 1;
    -4, 12, -12, 4, 0;
    6, -12, 6, 0, 0;
    -4, 4, 0, 0, 0;
    1, 0, 0, 0, 0];
uVals = 0:0.025:1;
U = [uVals.^4; uVals.^3; uVals.^2; uVals; ones(1, length(uVals))]';
P = U*MsubB*B;

% plotting
figure()
subplot(1, 2, 1), plot3(get(lineB, 'xdata'), get(lineB, 'ydata'), get(lineB, 'zdata')), grid on % reusing the part b plot data
subplot(1, 2, 2), plot3(P(:, 1), P(:, 2), P(:, 3)), grid on
a = axes; t1 = title('Problem 4c - Global Affect of Modifying Bezier Control Points'); set(a, 'Visible', 'off'); set(t1, 'Visible', 'on'); % creating unified subplot title

fprintf('Problem 4c - The curve is affected globally\n')