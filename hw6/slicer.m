function [p] = slicer(f, v, z_slice, th, show)
%slicer 
%   Inputs:
%       f (faces) [3xn array] - indexes into the vertex array (ex.
%       [4, 6, 7] -> 4th, 6th, and 7th values of vertexes `v` 
%       v (vertexes) [3xn array] = [x, y, z] values used to compose
%       faces `f`
%       z_slice (z height) [uint] - the z height to slice at
%       th (thickness) [uint] - thickness of slice in mm
%       show [bool] - whether or not to plot
%   Outputs:
%       p (points) [3xn array] - the points on z slice plane
%% Initialize
if nargin < 4
    error('not enough arguments')
elseif nargin == 4
    show = 1;
end
%% Step II
z_max = max(v(:, 3));
z_min = min(v(:, 3));
layer_cnt = (z_max - z_min) / th;                           % number of layers
%% Step III
target_faces = [];
for i = 1:size(f, 1)                                        % for face
    v_i = f(i, :);                                          % get 3 points (ex. [104, 106, 107])
    z = v(v_i, 3);                                          % get z-val for 3 pts
    if min(z) <= z_slice && max(z) >= z_slice               % if some vertices of face lower than z_slice and some higher
        target_faces = [target_faces; f(i, :)];             % then that face is part of slice
    end
end
%% Step IV
t = zeros(size(target_faces, 1), 3);                        % preallocate
for i = 1:size(target_faces, 1)                             % for target_face
    v_i = target_faces(i, :);                               % get 3 points (ex. [104, 106, 107])
    z = v(v_i, 3);                                          % get z-val for 3 pts
    tt = zeros(1, 3);                                       % preallocate
    tt(1) = (z_slice - z(1)) / (z(2) - z(1));               % ratio along edge 1
    tt(2) = (z_slice - z(2)) / (z(3) - z(2));               % ratio along edge 2
    tt(3) = (z_slice - z(3)) / (z(1) - z(3));               % ratio along edge 3
    t(i, :) = tt;                                           % store ratios along 3 edges for each target_face
end
%% Step V
r = zeros(3, 3, size(target_faces, 1));                     % rows:vertexes, cols:directions, aisles:faces
for i = 1:size(target_faces, 1)
    v_i = target_faces(i, :);                               % get 3 points (ex. [104, 106, 107])
    x = v(v_i, 1);                                          % get x-val for 3 pts
    y = v(v_i, 2);                                          % get y-val for 3 pts
    z = v(v_i, 3);                                          % get z-val for 3 pts
    tt = t(i, :);
    
    if tt(1) >= 0 && tt(1) <= 1
        r(1, 1, i) = x(1) + tt(1) * (x(2) - x(1));
        r(1, 2, i) = y(1) + tt(1) * (y(2) - y(1));
        r(1, 3, i) = z(1) + tt(1) * (z(2) - z(1));
    else
        r(1, :, i) = zeros(1, 3);
    end
    if tt(2) >= 0 && tt(2) <= 1
        r(2, 1, i) = x(2) + tt(2) * (x(3) - x(2));
        r(2, 2, i) = y(2) + tt(2) * (y(3) - y(2));
        r(2, 3, i) = z(2) + tt(2) * (z(3) - z(2));
    else
        r(2, :, i) = zeros(1, 3);
    end
    if tt(3) >= 0 && tt(3) <= 1
        r(3, 1, i) = x(3) + tt(3) * (x(1) - x(3));
        r(3, 2, i) = y(3) + tt(3) * (y(1) - y(3));
        r(3, 3, i) = z(3) + tt(3) * (z(1) - z(3));
    else
        r(3, :, i) = zeros(1, 3);
    end
end
%% Step VI
for i = 1:size(r, 3)                                        % for face indexes
    vv = r(:, :, i);                                        % rows:vertexes, cols:directions
    if all(vv(1, :) == vv(2, :))                            % if two vertexes (in the same face) are on the slicing plane and equal
        r(1, :, i) = zeros(1, 3);                           % then drop one of the vertexes (zeros will drop out later)
    elseif all(vv(1, :) == vv(3, :))
        r(1, :, i) = zeros(1, 3);
    elseif all(vv(2, :) == vv(3, :))
        r(2, :, i) = zeros(1, 3);
    end
    if sum(~all(vv, 2)) > 1                                 % if the face only has one intersection with the slicing plane
        r(:, :, i) = zeros(3, 3);                           % then drop the face (zeros will drop out later)
    end
end
%% Step VII
p = zeros(size(r, 1)*size(r, 3), size(r, 2));               % 3D array -> 2D array (couldn't get reshape to work how I needed)
for i = 1:size(r, 3)
    j = (i - 1) * 3 + 1;                                    % i = 1, j = 1 -> 3; i = 2, j = 4 -> 6; etc.
    p(j:j+2, :) = r(:, :, i);
end
% p_indx = find(~all(p, 2) == 0);                             % find nonzero rows
% p = p(p_indx, :); % mask out/drop zero-filled rows
% p = unique(p, 'stable', 'rows');                            % drop non-unique rows/vertexes ('stable' stops unique from sorting)

contours = [];
for i = 1:size(p, 1)
    pToFind = p(i, :);
    [~, row_indx] = ismember(pToFind, p, 'rows')
end
    
%% Post
if show == 1
    figure()
    
    subplot(1, 2, 1)
    plot3(p(:, 1), p(:, 2), p(:, 3), 'r*'), hold on
    patch('vertices', v, 'faces', f, 'facecolor', 'k'), hold off
    title('Vertexes on Slicing Plane on Object')
    grid on, axis square
    
    subplot(1, 2, 2)
    plot(p(:, 1), p(:, 2), 'r')
    title('2D Profile of Object on Slicing Plane')
    axis square
end

