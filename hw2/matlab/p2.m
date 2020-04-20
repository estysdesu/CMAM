% #####
% Homework 2 Problem 2
% Tyler Estes
% Utkarsh Wani
% #####

%% p2

verts = [5, 1, 0, 1;
    5, 1, 2, 1;
    6, 1, 1, 1;
    4, 1, 1, 1;
    5, 0, 1, 1;
    5, 2, 1, 1]; % vertices

% verts = [5, 1, 0, 1;
%     6, 1, 1, 1;
%     5, 2, 1, 1]; % vertices

% translate V1 to origin
tx1 = -5; ty1 = -1; tz1 = 0;
T1 = [1, 0, 0, 0;
    0, 1, 0, 0;
    0, 0, 1, 0;
    tx1, ty1, tz1, 1];

% rotate V3 to x-axis (rotation is about y-axis)
th2 = 45; % degs
T2 = [cosd(th2), 0, -sind(th2), 0;
    0, 1, 0, 0;
    sind(th2), 0, cosd(th2), 0;
    0, 0, 0, 1];

% rotate V6 to x-z plane (rotation is about x-axis)
th3 = (90-atand(1/sqrt(2)/1)); % degs
T3 = [1, 0, 0, 0;
    0, cosd(th3), sind(th3), 0;
    0, -sind(th3), cosd(th3), 0;
    0, 0, 0, 1];

% calculate transformed vertices
vertsTrans = verts*T1*T2*T3;
% project onto x-z plane (drop y aka second column)
projXZ = vertsTrans; projXZ(:, 2) = zeros(length(projXZ), 1);
% calculate the inverse transforms back to original coordinate axes
projV1V3V6 = projXZ*inv(T3)*inv(T2)*inv(T1);

figure()
plot3(projV1V3V6(:, 1), projV1V3V6(:, 2), projV1V3V6(:, 3), 'r*'), hold on % projected points
patch(verts([1, 3, 6], 1), verts([1, 3, 6], 2), verts([1, 3, 6], 3), 'blue'), hold off % plane of V1, V3, V6
grid on
