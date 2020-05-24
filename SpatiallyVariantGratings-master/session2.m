% Session1.m 

%% INITIALIZE MATLAB 

close all; 
clc; 
clear all; 

% UNITS 
degrees = pi/180; 

% OPEN FIGURE WINDOW 
figure('Color' , 'w', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PARAMETERS FOR THE PLANAR GRATING

a     = 1;          % period of grating 


% GRID PARAMETERS 
Sx = 10*a; 
Sy = 10*a; 
NRESLO = 6; 
NRESHI = 10; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OUR GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOW RESOLUTION GID 
dx = a/NRESLO; 
dy = a/NRESLO; 


% NUMBER OF POINTS 
Nx = ceil(Sx/dx); 
dx = Sx/Nx; 

Ny = ceil(Sx/dy); 
dy = Sy/Ny; 

% GRID AXES
xa = [0:Nx-1]*dx; xa = xa-mean(xa); 
ya = [0:Ny-1]*dy; ya = ya-mean(ya); 

% MESHGRID
[Y, X] = meshgrid(ya, xa); 


% HIGH RESOLUTION GID 
dx2 = a/NRESHI; 
dy2 = a/NRESHI; 


% NUMBER OF POINTS 
Nx2 = ceil(Sx/dx2); 
dx2 = Sx/Nx2; 

Ny2 = ceil(Sx/dy2); 
dy2 = Sy/Ny; 

% GRID AXES
xa2 = linspace(xa(1), xa(Nx), Nx2);
ya2 = linspace(ya(1), ya(Ny), Ny2); 

dx2 = xa2(2)- xa2(1); 
dy2 = ya2(2) - ya2(1); 
% MESHGRID
[Y2, X2] = meshgrid(ya2, xa2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE INPUT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PERIOD 
PER = X-Y; 
PER = PER - min(PER(:)); 
PER = PER/max(PER(:)); 
PER = .5*a + a *PER; 

% ORIENTATION PHI
r = .25*Sx; 
RSQ = X.^2+Y.^2; 
THETA = (45*degrees)*(RSQ<r.^2);
THETA = svlblur(THETA, [2 2]); 

% FILL FACTOR 
FF = X2+Y2; 
FF = FF-min(FF(:));
FF = FF/max(FF(:)); 
FF = .2 + .6*FF; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE UNIFORM SPATIALLY-VARIANT GRATING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CONSTRUCT K-FUNCTIONS 

Kx = (2*pi./PER) .* cos(THETA); 
Ky = (2*pi./PER) .* sin(THETA);
% 
% CONSTRUCT DERIVATIVE OPERATION 
NS             = [Nx, Ny]; 
RES            = [dx, dy];
BC             = [1,1]; 
[DX, ~, DY, ~] = fdder(NS, RES, BC); 

% SOLVE FOR GRATING PHASE 
A = [ DX; DY]; 
b = [Kx(:) ; Ky(:)];

PHI = (A.'*A)\(A.'*b); 
PHI = reshape(PHI, Nx, Ny); 

% INTERPOLATE TO HIGH RESOLUTION GRID
PHI = interp2(ya,xa', PHI, ya2, xa2'); 

% CALCULATE ANALOG GRATING 
GA = cos(PHI);

% CALCULATE BINARY GRATING; 
 gth = cos(pi*FF); 
 GB = GA>gth; 
 GB = double(GB); 
 
%% MESHING GRATING TO STL FILE 
% STACK IMAGE
GH = GB; 
GH(:, :, 2) = GB; 
% CREATE 2D MESH
[F4, V4] = isocaps(ya2, xa2, [0, 1], GH, .5, 'zmin'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHOW GRATINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SHOW ANALOG GRATING 
subplot(231); 
pcolor(xa,ya, PER'); 
shading interp; 
colorbar; 
axis equal tight; 
title('PERIOD'); 
colormap(jet); 


subplot(232); 
pcolor(xa,ya, THETA'); 
shading interp; 
colorbar; 
axis equal tight; 
title('THETA'); 
colormap(jet); 

subplot(233); 
pcolor(xa2,ya2, FF'); 
shading interp; 
colorbar; 
axis equal tight; 
title('FF'); 
colormap(jet); 

subplot(234); 
pcolor(xa2, ya2, PHI'); 
shading interp; 
colorbar; 
axis equal tight; 
title('PHI'); 
colormap(jet); 

subplot(235); 
pcolor(xa2, ya2, GA'); 
shading interp; 
colorbar; 
axis equal tight; 
title('Analog Grating'); 
colormap(jet); 

subplot(236); 
h = imagesc(xa2, ya2, GB'); 
h2 = get(h, 'Parent'); 
set(h2, 'YDir', 'normal'); 
shading interp; 
colorbar; 
axis equal tight; 
title('Binary Grating'); 
colormap(jet); 


figure2 = figure('Color', 'w'); 
h4 = patch('faces', F4, 'vertices', V4);
h5 = get(h4, 'Parent'); 
set(h5, 'YDir', 'normal');
axis xy;
svlcad('GRATING.stl', F4, V4); 