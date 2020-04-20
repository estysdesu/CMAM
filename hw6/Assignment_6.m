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
P = slicer(F, V, 5, 1);