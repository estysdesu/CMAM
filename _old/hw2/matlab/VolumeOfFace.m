function vol = VolumeOfFace(varargin)
    if nargin == 1 % input is of 2d matrix form
        m = cell2mat(varargin);
    else % input is series of 1d matrix form
        m = cell2mat(varargin');
    end

    success = 0;
    i = 1;
    while success == 0
        p1 = m(i, :);
        p2 = m(i+1, :);
        if or(length(m) < i+2, any(m(i+2, :) == NaN))
            error('no more points on face')
        end
        p3 = m(i+2, :);

        try % try finding the normal of the face from 3 points (2 vectors)
            n = faceUnitNormalVector(p1, p2, p3); % unit normal vector
            success = 1;
        catch % if there is an error, check the next 3 points (windowed)
            i = i + 1;
        end

        a = area(n, p1, p2, p3); % area of face
        vol = p1 * n' .* a; % volume
    end


function unitNormVec = faceUnitNormalVector(p1, p2, p3)
    v1 = p2 - p1; % vector 1
    v2 = p3 - p1; % vector 2
    normVec = cross(v1, v2); % cross product of vector 1 w/ vector 2 yields normal vector
    if normVec == 0
        error('vectors in same direction')
    end
    unitNormVec = normVec/norm(normVec); % normalized normal vector yields unit normal vector

function A = area(norm_vec, p1, p2, p3)
    m = [norm_vec; (p3 - p1); (p2 - p3)]; % build determinant matrix
    A = abs(det(m)); % positive determinant yields area