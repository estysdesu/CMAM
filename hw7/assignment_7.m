clear all
close all
clc

%% read the stl file and plot the rotated part.
%% Read the information from the stl file i.e Faces, Normals and Vertices
[F,V,N] = stlread('Part_hw.stl');
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

%% step-I, creation of substrate/base plate
%% for the given build orientation find the maximum X and Y and minimum Z.
%% Assume the support substrate is placed 10 mm below the lowermost Z-coordinate the transformed part. (Zmin-10)

%% Calculating xmin ymin zmin
xmin= min(V2(:,1)); % minimum x coordinate
ymin= min(V2(:,2)); % minimum y coordinate
zmin= min(V2(:,3)); % minimum z coordinate
zmin= zmin-10;      %support substrate placed 10mm below lowermost z-coordinate

%% calculating xmax ymax
xmax= max(V2(:,1)); % maximum of x coordinate
ymax= max(V2(:,2)); % maximum of y coordinate

%% Plotting the substrate/base plate
figure()
patch([xmin xmax xmax xmin xmin],[ymin ymin ymax ymax ymin],[zmin-10 zmin-10 zmin-10 zmin-10 zmin-10]);
grid off
title('Substrate')

%% shoot rays from the baseplate plane in +z direction

% x= xmin:2:xmax; % variation of x
% y= ymin:2:ymax;
% [X,Y] = meshgrid(x,y);
% 
% for i= 1:length(X(:,1))
%     for j= 1:length(Y(1,:))
%         plot(X(i,:),Y(i,:));
%         hold on
%         plot(X(:,j),Y(:,j));
%         hold on
%     end
% end
% 
% figure()
% patch('vertices', V2, 'faces', F, 'facevertexcdata', jet(length(F)), 'facecolor', 'flat')
% hold on
% for i= 1:length(X(:,1))
%     for j= 1:length(Y(1,:))
%         plot(X(i,:),Y(i,:),'k');
%         hold on
%         plot(X(:,j),Y(:,j),'k');
%         hold on
%     end
% end
% 
% 
% 
% FF=[]
% for i= 1:length(X)
%     for j=1:length(Y(:,1))
%         for k= 1:length(F(:,1))
%             v_i = F(k, :);
%             x = V2(v_i(1:length(v_i)), 1);                                          % get x-val for 3 pts
%             y = V2(v_i(1:length(v_i)), 2);
%             X1= X(j,i); % Extracting x coordinate from the mesh grid
%             Y1= Y(j,i);  % % Extracting y coordinate from the mesh grid
%             A = [x(1),y(1)];
%             B = [x(2),y(2)];
%             C = [x(3),y(3)];
%             % collect coordinates
%             xx = [A(1), B(1), C(1)]; % x values
%             yy = [A(2), B(2), C(2)]; % y values
%             % compute area
%             triangle_area = polyarea(xx,yy);
%             
%             A1= [x(1),y(1)];
%             B1= [x(2),y(2)];
%             D1= [X1, Y1];
%             % collect coordinates
%             xx1 = [A1(1), B1(1), D1(1)]; % x values
%             yy1 = [A1(2), B1(2), D1(2)]; % y values
%             % compute area
%             triangle_area_1 = polyarea(xx1,yy1);
%             
%             A2= [x(2),y(2)];
%             B2= [x(3),y(3)];
%             D2= [X1, Y1];
%             % collect coordinates
%             xx2 = [A2(1), B2(1), D2(1)]; % x values
%             yy2 = [A2(2), B2(2), D2(2)]; % y values
%             % compute area
%             triangle_area_2 = polyarea(xx2,yy2);
%             
%             A3= [x(1),y(1)];
%             B3= [x(3),y(3)];
%             D3= [X1, Y1];
%             % collect coordinates
%             xx3 = [A3(1), B3(1), D3(1)]; % x values
%             yy3 = [A3(2), B3(2), D3(2)]; % y values
%             % compute area
%             triangle_area_3 = polyarea(xx3,yy3);
%             summation= triangle_area_1+triangle_area_2+triangle_area_3;
%             
%             if abs(triangle_area-summation) < 1e-5
%                 FF=[FF;F(k,:)]
%             end
%         end
%     end
% end

%% 

[meshX, meshY] = meshgrid(xmin:2:xmax, ymin:2:ymax);

hitFaces = zeros(size(meshX, 1), size(meshX, 2), size(F, 1));
for i = 1:size(meshX, 1) % meshX and meshY are same size
    for j = 1:size(meshX, 2)
        for k = 1:size(F, 1)
            vv = V2(F(k, :), 1:2); % 3x2
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
            if abs(multiTriangleArea - triangle_area) < 1e-5
                hitFaces(i, j, k) = 1;
            end
        end
    end
end

hitFaces