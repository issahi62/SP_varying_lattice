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
theta = 0*degrees; % slant of the grating
ff    = .6; 


% GRID PARAMETERS 
Sx = 10*a; 
Sy = 10*a; 
NRES = 10; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OUR GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx = a/NRES; 
dy = a/NRES; 

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE UNIFORM PLANAR GRATING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE GRATING VECTORS FUNCTION 

Kx = (2*pi/a) * cos(theta); 
Ky = (2*pi/a) * sin(theta); 

% CALCULATE ANALOG GRATING

GA = cos(Kx*X + Ky*Y);

% CALCULATE BINARY GRATING; 
gth = cos(pi*ff); 
GB = GA>gth; 
GB = double(GB); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHOW GRATINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SHOW ANALOG GRATING 
subplot(121); 
pcolor(xa,ya, GA'); 
shading interp; 
colorbar; 
axis equal tight; 
title('ANALOG GRATING'); 


subplot(122); 
pcolor(xa,ya, GB'); 
shading interp; 
colorbar; 
axis equal tight; 
title('BINARY GRATING'); 


