clc, clear, close all
%% step -I
% Read the information from the stl file i.e Faces, Normals and Vertices
[F,V,N] = stlread('part.stl');

% Plotting of the original object using patch command
figure()
patch('vertices', V, 'faces', F, 'facevertexcdata', jet(length(F)), 'facecolor', 'flat')
grid on
title('Original Object')

% Rotation about X-axis
theta = 30; % 30 degress about X-axis
Rx = [1, 0, 0; 0, cosd(theta),-sind(theta); 0, sind(theta),cosd(theta)]; % Rotation Matrix along X-axis

% for loop to update the vertices after rotation along X axis by 30 degrees
for i= 1:length(V)
    V1(i,:) = V(i,:)*Rx;
end

% Rotation about Y-axis
theta= 45; % 45 degrees about Y-axis
Ry= [cosd(theta), 0 , sind(theta); 0,1,0; -sind(theta) 0, cosd(theta)]; % Rotation Matrix along Y-axis
% for loop to update the vertices after rotation along Y axis by 45 degrees
for i= 1:length(V1)
    V2(i,:) = V1(i,:)*Ry;
end

% Plotting of the rotated object using patch command
figure()
patch('vertices', V2, 'faces', F, 'facevertexcdata', jet(length(F)), 'facecolor', 'flat')
grid on
title('Rotated Object')

%% step-II
% Maximum of Z amongst all the vertices 
Z_max = max(V2(:,3));

% Minimum of Z amongst all the vertices
Z_min= min(V2(:,3));

% Number of Layers
Distance=1; % Given in the problem, uniform layer thickness of 1 mm
N_o_L = (Z_max-Z_min)/Distance; % Number of Layers
%% Step-III
% For the rotated part, find intersection of a slicing planes at Z=5mm and Z=30mm
Z_target = 5;
target_faces_5 = [];
for i = 1:length(F)
    v_i = F(i, :); % F(i) = some array like 80, 83, 97
    for j = 1:length(v_i)
        z(j) = V2(v_i(j), 3); % v_i(j) = some number like 104
    end
    if min(z) <= Z_target && max(z) >= Z_target
        target_faces_5 = [target_faces_5; F(i)];
    end
end

Z_target = 30;
target_faces_30 = [];
for i = 1:length(F)
    v_i = F(i, :); % F(i) = some array like 80, 83, 97
    for j = 1:length(v_i)
        z(j) = V2(v_i(j), 3); % v_i(j) = some number like 104
    end
    if min(z) <= Z_target && max(z) >= Z_target
        target_faces_30 = [target_faces_30; F(i)];
    end
end
