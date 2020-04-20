clear all
close all

%% step -I
%% Read the information from the stl file i.e Faces, Normals and Vertices
[F,V,N] = stlread('part.stl');
%% Plotting of the original object using patch command
figure()
patch('vertices', V, 'faces', F, 'facevertexcdata', jet(length(F)), 'facecolor', 'flat')
grid on
title('Original Object')


%% Rotation about X-axis
theta= 30; % 30 degress about X-axis
Rx= [1, 0, 0; 0, cosd(theta),-sind(theta); 0, sind(theta),cosd(theta)]; % Rotation Matrix along X-axis
%% for loop to update the vertices after rotation along X axis by 30 degrees
 for i= 1:length(V)
V1(i,:)= V(i,:)*Rx;
 end

%% Rotation about Y-axis
theta= 45; % 45 degrees about Y-axis
Ry= [cosd(theta), 0 , sind(theta); 0,1,0; -sind(theta) 0, cosd(theta)]; % Rotation Matrix along Y-axis
%% for loop to update the vertices after rotation along Y axis by 45 degrees
 for i= 1:length(V1)
V2(i,:)= V1(i,:)*Ry;
 end
%% Plotting of the rotated object using patch command
figure()
patch('vertices', V2, 'faces', F, 'facevertexcdata', jet(length(F)), 'facecolor', 'flat')
grid on
title('Rotated Object')

%% step-II
%% Maximum of Z amongst all the vertices 
Z_max = max(V2(:,3));

%% Minimum of Z amongst all the vertices
Z_min= min(V2(:,3));

%% Number of Layers
Distance=1; % Given in the problem, uniform layer thickness of 1 mm
N_o_L = (Z_max-Z_min)/Distance; % Number of Layers

%% Step-III
%% For the rotated part,intersection of a slicing planes at Z=5mm and Z=30mm
%% Step-III

%% For z= 5mm 
% For the rotated part, find intersection of a slicing planes at Z=5mm and Z=30mm
Z_target = 5;
target_faces_5 = [];
for i = 1:length(F)
    v_i = F(i, :); % F(i) = some array like 80, 83, 97
    for j = 1:length(v_i)
        z(j) = V2(v_i(j), 3); % v_i(j) = some number like 104
    end
    if min(z) <= Z_target && max(z) >= Z_target
        target_faces_5 = [target_faces_5;F(i,:)];
    end
end

%% For z= 30mm 
Z_target = 30;
target_faces_30 = [];
for i = 1:length(F)
    v_i = F(i, :); % F(i) = some array like 80, 83, 97
    for j = 1:length(v_i)
        z(j) = V2(v_i(j), 3); % v_i(j) = some number like 104
    end
    if min(z) <= Z_target && max(z) >= Z_target
        target_faces_30 = [target_faces_30; F(i,:)];
    end
end

%% step IV 
%Finding the intersection points of the slicing plane with the STL facets 
%% Slicing height, Z=5mm
Z_slicing_plane= 5;
T_5 = [];
for i= 1:length(target_faces_5) % extracting all the faces
    v = V2(target_faces_5(i, :), :); %% extracting all the vertices for each face (x1 x2 x3, y1 y2 y3, z1 z2 z3)
    z = v(:, 3); %% Extracting the Z coordinate from all the vertices
    
    t(1) = (Z_slicing_plane - z(1)) / (z(2) - z(1));  % Intersection point on edge 1
    t(2) = (Z_slicing_plane - z(2)) / (z(3) - z(2));  % Intersection point on edge 2
    t(3) = (Z_slicing_plane - z(3)) / (z(1) - z(3));  % Intersection point on edge 3
    T_5 = [T_5; t];
end

%% Slicing height, Z=30mm
Z_slicing_plane= 30;
T_30 = [];
for i= 1:length(target_faces_30) % extracting all the faces
    v = V2(target_faces_30(i, :), :); %% extracting all the vertices for each face (x1 x2 x3, y1 y2 y3, z1 z2 z3)
    z = v(:, 3); %% Extracting the Z coordinate from all the vertices
    
    t(1) = (Z_slicing_plane - z(1)) / (z(2) - z(1));  % Intersection point on edge 1
    t(2) = (Z_slicing_plane - z(2)) / (z(3) - z(2));  % Intersection point on edge 2
    t(3) = (Z_slicing_plane - z(3)) / (z(1) - z(3));  % Intersection point on edge 3
    T_30 = [T_30; t];
end

%% Step-V
%% Slicing plane z = 5; Find the intersection point r for every edge per face 

for i = 1:length(target_faces_5) % extracting the faces
    v = V2(target_faces_5(i, :), :); % extracting all the vertices for each face
    x = v(:, 1); % x1 x2 x3 for each face
    y = v(:, 2); % y1 y2 y3 for each face
    z = v(:, 3); % z1 z2 z3 for each face
    t = T_5(i, :); % Extracting the point of intersection for each face
    
    if t(1) >= 0 && t(1) <= 1
        r1 = [x(1) + t(1) * (x(2) - x(1)), y(1) + t(1) * (y(2) - y(1)), z(1) + t(1) * (z(2) - z(1))];
    else
        r1 = [0, 0, 0];
    end
    if t(2) >= 0 && t(2) <= 1
        r2 = [x(2) + t(2) * (x(3) - x(2)), y(2) + t(2) * (y(3) - y(2)), z(2) + t(2) * (z(3) - z(2))];
    else
        r2 = [0, 0, 0];
    end
    if t(3) >= 0 && t(3) <= 1
        r3 = [x(3) + t(3) * (x(1) - x(3)), y(3) + t(3) * (y(1) - y(3)), z(3) + t(3) * (z(1) - z(3))];
    else
        r3 = [0, 0, 0];
    end
    r_target_5(:, :, i) = [r1; r2; r3]; % rows = vertexes, col = directions, i = faces
end

%% Slicing plane z = 30; Find the intersection point r for every edge per face 

for i = 1:length(target_faces_30) % extracting the faces
    v = V2(target_faces_30(i, :), :); % extracting all the vertices for each face
    x = v(:, 1); % x1 x2 x3 for each face
    y = v(:, 2); % y1 y2 y3 for each face
    z = v(:, 3); % z1 z2 z3 for each face
    t = T_30(i, :); % Extracting the point of intersection for each face
    
    if t(1) >= 0 && t(1) <= 1
        r1 = [x(1) + t(1) * (x(2) - x(1)), y(1) + t(1) * (y(2) - y(1)), z(1) + t(1) * (z(2) - z(1))];
    else
        r1 = [0, 0, 0];
    end
    if t(2) >= 0 && t(2) <= 1
        r2 = [x(2) + t(2) * (x(3) - x(2)), y(2) + t(2) * (y(3) - y(2)), z(2) + t(2) * (z(3) - z(2))];
    else
        r2 = [0, 0, 0];
    end
    if t(3) >= 0 && t(3) <= 1
        r3 = [x(3) + t(3) * (x(1) - x(3)), y(3) + t(3) * (y(1) - y(3)), z(3) + t(3) * (z(1) - z(3))];
    else
        r3 = [0, 0, 0];
    end
    r_target_30(:, :, i) = [r1; r2; r3]; % rows = vertexes, col = directions, i = faces
end

%% Step VI
%% special cases
%% For z= 5mm
for i = 1:size(r_target_5, 3) % for face indexes
    v = r_target_5(:, :, i); % each vertex is a row
    if all(v(1, :) == v(2, :))
        r_target_5(1, :, i) = [0, 0, 0];
    elseif all(v(1, :) == v(3, :))
        r_target_5(1, :, i) = [0, 0, 0];
    elseif all(v(2, :) == v(3, :))
        r_target_5(2, :, i) = [0, 0, 0];
    end
    if sum(~all(v, 2)) > 1
        r_target_5(:, :, i) = zeros(3, 3) % make zeros b/c zeros will drop out during contour tracing
    end
end

%% For Z=30mm
%% For z= 5mm
for i = 1:size(r_target_30, 3) % for face indexes
    v = r_target_30(:, :, i); % each vertex is a row
    if all(v(1, :) == v(2, :))
        r_target_30(1, :, i) = [0, 0, 0];
    elseif all(v(1, :) == v(3, :))
        r_target_30(1, :, i) = [0, 0, 0];
    elseif all(v(2, :) == v(3, :))
        r_target_30(2, :, i) = [0, 0, 0];
    end
    if sum(~all(v, 2)) > 1
        r_target_30(:, :, i) = zeros(3, 3) % make zeros b/c zeros will drop out during contour tracing
    end
end

%% Step VII
%% contour construction
%% Z=5mm
% creating zeros for the C matrix of size (3*3*83)
for i= length(r_target_5)
    C_5= zeros(3,3,i);
end

D_5 = zeros( size(C_5, 1)*size(C_5, 3), size(C_5, 2) ); % size(D) = (249, 3)
for i=1:size(C_5, 3)
    j = (i-1)*3+1; % i = 1, j = 1->3; i = 2; j = 4->6
    D_5(j:j+2, :) = r_target_5(:, :, i);
end

%% drop zeros
nonzeroes_index_5 = find(~all(D_5, 2) == 0);
non_zeroes_D_5 = D_5(nonzeroes_index_5, :) % Non zero rows of D are obtained 

%% drop non-unique rows
unique_D_5 = unique(non_zeroes_D_5, 'stable', 'rows') % Unique rows of D

%% Z= 30mm
% creating zeros for the C matrix of size (3*3*83)
for i= length(r_target_30)
    C_30= zeros(3,3,i);
end

D_30 = zeros( size(C_30, 1)*size(C_30, 3), size(C_30, 2) ); % size(D) = (186, 3)
for i=1:size(C_30, 3)
    j = (i-1)*3+1; % i = 1, j = 1->3; i = 2; j = 4->6
    D_30(j:j+2, :) = r_target_30(:, :, i);
end

%% drop zeros
nonzeroes_index_30 = find(~all(D_30, 2) == 0);
non_zeroes_D_30 = D_30(nonzeroes_index_30, :) % Non zero rows of D are obtained 

%% drop non-unique rows
unique_D_30 = unique(non_zeroes_D_30, 'stable', 'rows') % Unique rows of D

% figure(3)
% % plot(non_zeroes_D_5(:,1), non_zeroes_D_5(:,2))
% % plot(unique_D_5(:,1), unique_D_5(:,2))

figure(4)
% plot(non_zeroes_D_5(:,1), non_zeroes_D_5(:,2))
plot3(unique_D_5(:,1), unique_D_5(:,2),unique_D_5(:,3),'y' )
hold on
patch('vertices', V2, 'faces', F, 'facevertexcdata', jet(length(F)), 'facecolor', 'flat')




%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
    
    
    
    % %% Step IV
% 
% Z_slicing_plane= 5;
% for i= 1:length(target_faces_5) % for face
%     V_i = target_faces_5(i,:);
%     z = V2(V_i, 3); % Z coordinates can be extracted from each of the vertices
%         z1= z(1,1);
%         z2= z(2,1);
%         z3= z(3,1);
%      t(1)= (Z_slicing_plane-z1)/(z2-z1); 
%      t(2)= (Z_slicing_plane-z2)/(z3-z2);
%      t(3)= (Z_slicing_plane-z3)/(z1-z3);
%      T{i,1}=t; % Cell has been created with the to store the intersection points for each facet
% end
% 
% Z_slicing_plane= 30;
% for i= 1:length(target_faces_30) % for face
%     V_i = target_faces_30(i,:);
%     z = V2(V_i, 3); % Z coordinates can be extracted from each of the vertices
%         z1= z(1,1);
%         z2= z(2,1);
%         z3= z(3,1);
%      t(1)= (Z_slicing_plane-z1)/(z2-z1); 
%      t(2)= (Z_slicing_plane-z2)/(z3-z2);
%      t(3)= (Z_slicing_plane-z3)/(z1-z3);
%      T1{i,1}=t; % Cell has been created with the to store the intersection points for each facet
% end
% %% Step V
% for i= 1:length(target_faces_5) % for face
%     v = target_faces_5(i,:);
%     x = V2(v, 1);
%     y= V2(v,2);
%     z= V2(v,3);
%     for j= 1:length(T) % for {t1, t2, t3} set
%         if T{j,1}(1)>=0 &&  T{j,1}(1)<=1
%             r1= [x(1,1)+T{j,1}(1)*(x(2,1)-x(1,1)) 
%                 y(1,1)+T{j,1}(1)*(y(2,1)-y(1,1))
%                 z(1,1)+T{j,1}(1)*(z(2,1)-z(1,1))]
%         else
%             r1= [0;0;0]
%         end
%         if T{j,1}(2)>=0 &&  T{j,1}(2)<=1
%             r2= [x(2,1)+T{j,1}(2)*(x(3,1)-x(2,1)); 
%                 y(2,1)+T{j,1}(2)*(y(3,1)-y(2,1)); 
%                 z(2,1)+T{j,1}(2)*(z(3,1)-z(2,1))];
%         else
%             r2= [0;0;0];
%         end
%         if T{j,1}(3)>=0 &&  T{j,1}(3)<=1
%             r3= [x(3,1)+T{j,1}(3)*(x(1,1)-x(3,1)); 
%                 y(3,1)+T{j,1}(3)*(y(1,1)-y(3,1)); 
%                 z(3,1)+T{j,1}(3)*(z(1,1)-z(3,1))];
%         else
%             r3= [0;0;0];
%         end
%     end
% end
    
    
    
