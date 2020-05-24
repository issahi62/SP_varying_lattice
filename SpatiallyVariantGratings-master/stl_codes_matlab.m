
%% sphere (3D stl files)
[X, Y, Z] = sphere(41); 
[F, V] = surf2patch(X, Y, Z, 'triangles');
subplot(3, 1,1); 
h = patch('faces', F, 'vertices', V);             
set(h, 'FaceColor', [.5 .5 .8], 'EdgeColor', 'k');

%% ellipsoids 
[X, Y, Z] = meshgrid(linspace(-1, 1)); 
rx = .2; ry=rx; rz=rx; 
OBJ =(X/rx).^2 +(Y/rx).^2 +(Z/rx).^2 ;
OBJ = (OBJ)<1; 
[F1, V1] = isosurface(X, Y, Z, OBJ, .5); 
subplot(312); 
h1 = patch('faces', F1, 'vertices', V1); 
set(h1, 'FaceColor', [.5 .5 .8], 'EdgeColor', 'k');

%% Using Isocaps to get 2D surface images
[F2, V2] = isocaps(X, Y, Z, OBJ, .5);
F3 = [ F1; F2+length(V1(:, 1))];
V3 = [ V1; V2]; 
subplot(3,1,3); 
h2 = patch('faces', F3, 'vertices', V3);
set(h2, 'FaceColor', [.5 .5 .8], 'EdgeColor', 'k');



%% images to stl files in matlab
% Load image
B = imread('imag.jpg'); 
Sx =1; 
Sy = Sx; 
% RESIZE IMAGE 
B = imresize(B, .2); 
[Nx, Ny, Nc] =size(B); 
dx = Sx/Nx; 
dy = Sy/Ny; 
xa = [1:Nx]*dx; 
ya = [1:Ny]*dy; 
%FLATTEN COLOR IMAGE 
B = double(B); 
B = B(:, :, 1) + B(:, :, 2) + B(:, :, 3);
B = 1-B/756; %% 3*255; 

% STACK IMAGE 
B(:, :, 2) = B; 

% CREATE 2D MESH

[F4, V4] = isocaps(ya, xa, [0, 1], B, .5, 'zmin'); 

%svlcad('Letter.stl', F4, v4); 
figure(2); 
h4 = patch('faces', F4, 'vertices', V4);
set(h4, 'FaceColor', [.5 .5 .8], 'EdgeColor', 'k');




