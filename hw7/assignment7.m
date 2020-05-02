%% Assignment 7
clc, clear all, close all

%% Read the STL file and rotate the part
[F, V, N] = stlread('part.stl');

figure()
patch('vertices', V, 'faces', F, 'facevertexcdata', jet(size(F, 1)), 'facecolor', 'flat')
grid on
title('Original Object')

thX = 30; % [deg]
Rx = [1, 0, 0; 0, cosd(thX), -sind(thX); 0, sind(thX), cosd(thX)]; % rotation matrix abt x-axis
thY = 45; % [deg]
Ry = [cosd(thY), 0 , sind(thY); 0,1,0; -sind(thY) 0, cosd(thY)]; % rotation matrix abt y-axis
R = Rx*Ry;
for i = 1:size(V, 1)
    V(i, :) = V(i, :)*R;
end
for i = 1:size(N, 1)
    N(i, :) = N(i, :)*R;
end

figure()
patch('vertices', V, 'faces', F, 'facevertexcdata', jet(length(F)), 'facecolor', 'flat')
grid on
title('Rotated Object')

%% Create substrate/base plate
Z_LOWER = 10; % [mm]
boundingBox = [
    min(V(:, 1)), max(V(:, 1)); % xLow, xHigh
    min(V(:, 2)), max(V(:, 2)); % yLow, yHigh
    min(V(:, 3))-Z_LOWER, max(V(:, 3)); % zLow - zDrop, zHigh
]';

figure()
ff = [1, 2, 3; 1, 5, 3; 4, 2, 3; 4, 5, 3];
patch('xdata', [boundingBox(1, 1), boundingBox(1, 1), boundingBox(2, 1), boundingBox(2, 1)], 'ydata', [boundingBox(1, 2), boundingBox(2, 2), boundingBox(2, 2), boundingBox(1, 2)], 'zdata', boundingBox(1, 3)*ones(1, 4))
grid off
title('Substrate Layer')

MESH_STEP = 2; % [mm]
[meshX, meshY] = meshgrid(boundingBox(1, 1):MESH_STEP:boundingBox(2, 1), boundingBox(1, 2):MESH_STEP:boundingBox(2, 2));

%% Find the faces that intersect each ray
SUPPORT_ANGLE = 135; % [deg]
BUILD_DIR = [0, 0, 1]; % [0, 0, 1] * [i, j, k]'
if exist('hitFaces.mat', 'file') == 2
    load('hitFaces.mat')
else
    hitFaces = zeros(size(meshX, 1), size(meshX, 2), size(F, 1), 2); # dims: [x, y, face, face_indx, dist]
    for i = 1:size(meshX, 1) % meshX and meshY are same size
        for j = 1:size(meshX, 2)
            for k = 1:size(F, 1)
                vv = V(F(k, :), :); % 3x3
                xBounds = [min(vv(:, 1)), max(vv(:, 1))];
                yBounds = [min(vv(:, 2)), max(vv(:, 2))];

                if meshX(i, j) < xBounds(1) || meshX(i, j) > xBounds(2) || meshY(i, j) < yBounds(1) || meshY(i, j) > yBounds(2) % ray not inside square bounding box
                    continue
                end

                A = vv(1, :);
                B = vv(2, :);
                C = vv(3, :);
                xx = [A(1), B(1), C(1)]; % x values
                yy = [A(2), B(2), C(2)]; % y values
                triangle_area = polyarea(xx, yy); % compute area

                A = vv(1, :);
                B = vv(2, :);
                C = [meshX(i, j), meshY(i, j)];
                xx = [A(1), B(1), C(1)]; % x values
                yy = [A(2), B(2), C(2)]; % y values
                triangle_area1 = polyarea(xx, yy); % compute area

                A = vv(2, :);
                B = vv(3, :);
                C = [meshX(i, j), meshY(i, j)];
                xx = [A(1), B(1), C(1)]; % x values
                yy = [A(2), B(2), C(2)]; % y values
                triangle_area2 = polyarea(xx, yy); % compute area

                A = vv(3, :);
                B = vv(1, :);
                C = [meshX(i, j), meshY(i, j)];
                xx = [A(1), B(1), C(1)]; % x values
                yy = [A(2), B(2), C(2)]; % y values
                triangle_area3 = polyarea(xx, yy); % compute area

                multiTriangleArea = triangle_area1 + triangle_area2 + triangle_area3;
                if abs(multiTriangleArea - triangle_area) > 1e-5 % areas not equal; ray not in triangle
                    continue
                end
                % hitFaces(i, j, k, 1) = 1;

                % A = vv(1, :);
                % B = vv(2, :);
                % C = vv(3, :);
                % AB = B - A;
                % BC = C - B;
                % n = cross(AB, BC);
                n = N(k, :);
                alpha = dot(BUILD_DIR, n);
                while alpha < 0
                    alpha = alpha + 360;
                end

                buildPlatePt = [meshX(i, j), meshY(i, j), boundingBox(1, 3)];
                d = dot( (A - buildPlatePt), n) / dot(BUILD_DIR, n);
                intersectPt = buildPlatePt + BUILD_DIR * d;
                hitFaces(i, j, k, 2) = intersectPt(3);
                if alpha < 135
                    hitFaces(i, j, k, 1) = 1;
                else
                    hitFaces(i, j, k, 1) = 2;
                end
            end
        end
    end
    save('hitFaces.mat', 'hitFaces')
end

figure(), grid on, hold on
title('Rotated Object with Generated Supports')
patch('xdata', [boundingBox(1, 1), boundingBox(1, 1), boundingBox(2, 1), boundingBox(2, 1)], 'ydata', [boundingBox(1, 2), boundingBox(2, 2), boundingBox(2, 2), boundingBox(1, 2)], 'zdata', boundingBox(1, 3)*ones(1, 4))
patch('vertices', V, 'faces', F, 'facevertexcdata', jet(length(F)), 'facecolor', 'flat')

vol = 0; % [mm]
for i = 1:size(hitFaces, 1)
    for j = 1:size(hitFaces, 2)
        intersectPtIndx = find(hitFaces(i, j, :, 1));
        overhangPtIndx = find(hitFaces(i, j, :, 1) == 2); % 1: face intersection with ray; 2: alpha > 135
        if size(overhangPtIndx, 1) > 0
            for k = 1:size(overhangPtIndx)
                if k == 1
                    bottomZ = boundingBox(1, 3);
                else
                    bottomZ = hitFaces(i, j, intersectPtIndx(k), 2);
                end
            topZ = hitFaces(i, j, overhangPtIndx(k), 2);
            vol = vol + MESH_STEP^2*(topZ - bottomZ);
            plot3(meshX(i,j)*ones(1, 2), meshY(i, j)*ones(1, 2), [bottomZ, topZ], 'color', [.8, .8, .8], 'linestyle', '--')
            end
        end
    end
end
hold off

fprintf('The volume of the generated supports is %.2f mm\n', vol)