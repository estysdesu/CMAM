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

%% Step IV

Z_slicing_plane= 5;
for i= 1:length(target_faces_5) % for face
    V_i = target_faces_5(i,:);
    z = V2(V_i, 3); % Z coordinates can be extracted from each of the vertices
        z1= z(1,1);
        z2= z(2,1);
        z3= z(3,1);
     t(1)= (Z_slicing_plane-z1)/(z2-z1); 
     t(2)= (Z_slicing_plane-z2)/(z3-z2);
     t(3)= (Z_slicing_plane-z3)/(z1-z3);
     T{i,1}=t; % Cell has been created with the to store the intersection points for each facet
end

Z_slicing_plane= 30;
for i= 1:length(target_faces_30) % for face
    V_i = target_faces_30(i,:);
    z = V2(V_i, 3); % Z coordinates can be extracted from each of the vertices
        z1= z(1,1);
        z2= z(2,1);
        z3= z(3,1);
     t(1)= (Z_slicing_plane-z1)/(z2-z1); 
     t(2)= (Z_slicing_plane-z2)/(z3-z2);
     t(3)= (Z_slicing_plane-z3)/(z1-z3);
     T1{i,1}=t; % Cell has been created with the to store the intersection points for each facet
end
%% Step V
for i= 1:length(target_faces_5) % for face
    v = target_faces_5(i,:);
    x = V2(v, 1);
    y= V2(v,2);
    z= V2(v,3);
    for j= 1:length(T) % for {t1, t2, t3} set
        if T{j,1}(1)>=0 &&  T{j,1}(1)<=1
            r1= [x(1,1)+T{j,1}(1)*(x(2,1)-x(1,1)); 
                y(1,1)+T{j,1}(1)*(y(2,1)-y(1,1)); 
                z(1,1)+T{j,1}(1)*(z(2,1)-z(1,1))];
        else
            r1= [0;0;0];
        end
        if T{j,1}(2)>=0 &&  T{j,1}(2)<=1
            r2= [x(2,1)+T{j,1}(2)*(x(3,1)-x(2,1)); 
                y(2,1)+T{j,1}(2)*(y(3,1)-y(2,1)); 
                z(2,1)+T{j,1}(2)*(z(3,1)-z(2,1))];
        else
            r2= [0;0;0];
        end
        if T{j,1}(3)>=0 &&  T{j,1}(3)<=1
            r3= [x(3,1)+T{j,1}(3)*(x(1,1)-x(3,1)); 
                y(3,1)+T{j,1}(3)*(y(1,1)-y(3,1)); 
                z(3,1)+T{j,1}(3)*(z(1,1)-z(3,1))];
        else
            r3= [0;0;0];
        end
    end
end




% %%
% Z_slicing_plane= 5;
% T = []
% for i= 1:length(target_faces_5) % for face
%     v = V2(target_faces_5(i, :), :);
%     z = v(:, 3);
%     
%     t(1) = (Z_slicing_plane - z(1)) / (z(2) - z(1));
%     t(2) = (Z_slicing_plane - z(2)) / (z(3) - z(2));
%     t(3) = (Z_slicing_plane - z(3)) / (z(1) - z(3));
%     T = [T; t];
% end
% 
% %%
% for i = 1:length(target_faces_5)
%     v = V2(target_faces_5(i, :), :);
%     x = v(:, 1);
%     y = v(:, 2);
%     z = v(:, 3);
%     t = T(i, :);
%     
%     if t(1) >= 0 && t(1) <= 1
%         r1 = [x(1) + t(1) * (x(2) - x(1)), y(1) + t(1) * (y(2) - y(1)), z(1) + t(1) * (z(2) - z(1))];
%     else
%         r1 = [0, 0, 0];
%     end
%     if t(2) >= 0 && t(2) <= 1
%         r2 = [x(2) + t(2) * (x(3) - x(2)), y(2) + t(2) * (y(3) - y(2)), z(2) + t(2) * (z(3) - z(2))];
%     else
%         r2 = [0, 0, 0];
%     end
%     if t(3) >= 0 && t(3) <= 1
%         r3 = [x(3) + t(3) * (x(1) - x(3)), y(3) + t(3) * (y(1) - y(3)), z(3) + t(3) * (z(1) - z(3))];
%     else
%         r3 = [0, 0, 0];
%     end
%     r = [r1; r2; r3]
% end
%     
%     
%     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
%       for j= 1:length(T)
%         v_t= T{j,:};
%         if v_t(1,1)>=0 && v_t(1,1)<=1
%     X_slicing_edge_1= X(1,1)+v_t(1,1)*(X(2,1)-X(1,1));
%     Y_slicing_edge_1= Y(1,1)+v_t(1,1)*(Y(2,1)-Y(1,1));
%     Z_slicing_edge_1= Z(1,1)+v_t(1,1)*(Z(2,1)-Z(1,1));
%         else
%     X_slicing_edge_1= 0;
%     Y_slicing_edge_1= 0;
%     Z_slicing_edge_1= 0;
%     
%     R_1= [ X_slicing_edge_1; Y_slicing_edge_1; Z_slicing_edge_1];
%     
%     if v_t(1,2)>=0 && v_t(1,2)<=1
%     X_slicing_edge_2= X(2,1)+v_t(2,1)*(X(3,1)-X(2,1))
%     Y_slicing_edge_2= Y(2,1)+v_t(2,1)*(Y(3,1)-Y(2,1))
%     Z_slicing_edge_2= Z(2,1)+v_t(2,1)*(Z(3,1)-Z(2,1))
%         else
%     X_slicing_edge_2= 0
%     Y_slicing_edge_2= 0
%     Z_slicing_edge_2= 0
%     
%    R_2= [ X_slicing_edge_2; Y_slicing_edge_2; Z_slicing_edge_2]
%     
%     
%     end
%         end
%        
%         end
%       end
% 
% 
%             
%     
%     
%     
%   
% 
% 
% % for i= 1:length(T) % Total span of the intersection
% %     v_t= T{i,:}
% %     x = V2(V_i, 3)
% %   for j = 1:length(v_t)  
% %       if v_t(j)>=0 && v_t(j)<=1
% %           
% % end
% % 
% % 
% 
%      
%         
%         
%         
%         
% %         
% % for j= 1:length(V_i) % for vertex
% %     Z(j)= V2(V_i(j),3)  
% %     for k =1:length(Z)
% %         
% % end
% % end
% 
