% #####
% Homework 2 Problem 1
% Tyler Estes
% Utkarsh Wani
% #####

%% p1
verts = [4, 6, -5;
    14, 6, -5;
    4, 16, -5;
    4, 6, 5;
    14, 16, 5;
    4, 16, 5;
    4, 6, 15;
    14, 6, 15;
    14, 16, 15]; % vertices

faces = [1, 3, 2, NaN, NaN;
    1, 2, 4, NaN, NaN;
    2, 3, 5, NaN, NaN;
    7, 6, 3, 1, 4;
    2, 5, 4, NaN, NaN;
    3, 6, 9, 5, NaN;
    4, 5, 8, NaN, NaN;
    4, 8, 7, NaN, NaN;
    5, 9, 8, NaN, NaN;
    6, 7, 9, NaN, NaN;
    7, 8, 9, NaN, NaN];

normDirs = [+1, -1, +1, -1, +1, +1, +1, -1, +1, -1, -1];

% plotting
figure()
patch('vertices', verts, 'faces', faces, 'facevertexcdata', jet(length(faces)), 'facecolor', 'flat')
grid on
title('Problem 1 - Area of Polyhedron')
view(3)

% sum the volume of each of the faces then divide by 1/6
total_vol = 0;
for n = 1:length(faces)
    face = faces(n, not(isnan(faces(n, :))));
    vol = VolumeOfFace(verts( face, :));
    vol = vol*normDirs(n);
    total_vol = total_vol + vol;
end
total_vol = total_vol*1/6;
fprintf('Problem 1 - Total Volume = %.2f\n', total_vol)
% error('not implemented: direction from origin')