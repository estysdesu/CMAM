clc, clear all, close all
%% Step I
[F, V] = stlread('part.stl');

figure()
patch('vertices', V, 'faces', F, 'facevertexcdata', jet(length(F)), 'facecolor', 'flat')
grid on
title('Original Object')

%% Rotation
thetaX = 30;
thetaY = 45;
Rx = [1, 0, 0;
    0, cosd(thetaX), -sind(thetaX);
    0, sind(thetaX),cosd(thetaX)]; % rotation x
Ry = [cosd(thetaY), 0 , sind(thetaY);
    0, 1, 0;
    -sind(thetaY) 0, cosd(thetaY)]; % rotation y
R = Rx*Ry;

for i = 1:size(V, 1)
    V(i, :) = V(i, :)*R;
end

figure()
patch('vertices', V, 'faces', F, 'facevertexcdata', jet(length(F)), 'facecolor', 'flat')
grid on
title('Rotated Object')
%% Step II - Step VII
P5 = slicer(F, V, 5, 1);
P30 = slicer(F, V, 30, 1);
%% Slice Area
SA5 = 0;
m = length(P5);
for i = 1:m-2
    a = [P5(i, :) - P5(m, :);
        P5(i+1, :) - P5(m, :)];
    n = cross(a(1, :), a(2, :));
    if n == zeros(1, 3)                                     % some normal vectors resulted in [0, 0, 0] because the vectors to P_m were in the same direction -- they shouldn't contribute to SA
        continue
    end
    u = n / norm(n);
    SA5 = SA5 + (0.5 * dot(u, cross(a(1, :), a(2, :))));
end
SA5 = abs(SA5)

SA30 = 0;
m = length(P30);
for i = 1:m-2
    a = [P30(i, :) - P30(m, :);
        P30(i+1, :) - P30(m, :)];
    n = cross(a(1, :), a(2, :));
    if n == zeros(1, 3)                                     % some normal vectors resulted in [0, 0, 0] because the vectors to P_m were in the same direction -- they shouldn't contribute to SA
        continue
    end
    u = n / norm(n);
    SA30 = SA30 + (0.5 * dot(u, cross(a(1, :), a(2, :))));
end
SA30 = abs(SA30)
