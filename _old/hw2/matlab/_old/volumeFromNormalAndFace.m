function [ volume ] = volumeFromNormalAndFace(vertexMatrix)
vertex1 = vertexMatrix(1, :);
vertex2 = vertexMatrix(2, :);
vertex3 = vertexMatrix(3, :);

normal = face_normal(vertex1, vertex2, vertex3)
area = area(normal, vertex1, vertex2, vertex3);
volume = vertex1 * normal' .* area;
end

function [ normal_vector ] = face_normal(vertex1, vertex2, vertex3)
inside1 = vertex2 - vertex1;
inside2 = vertex3 - vertex1;
normal_vector = cross(inside1, inside2)/norm(cross(inside1, inside2));
end

function [ area ] = area(normal_vector, vertex1, vertex2, vertex3 )
m = [normal_vector; (vertex3 - vertex1); (vertex2 - vertex3)];
area = abs(det(m));
end