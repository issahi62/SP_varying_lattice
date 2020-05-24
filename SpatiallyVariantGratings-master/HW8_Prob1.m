% HW8_Prob1.m
%
% Homework #8, Problem #1
% ECE 5322 – 21st Century Electromagnetics
% 21ST CENTURY ELECTROMAGNETICS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZE MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clc;
clear all;

% UNITS
degrees = pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEVICE PARAMETERS
a = 1;
Sx = 10*a;
Sy = 10*a;

% GRID PARAMETERS
dx = a/Sx;
Nx = Sx/dx;
dy = a/Sy;
Ny = Sy/dy;

% CREATE SOME AXES
xa      = [0 : Nx-1]*dx; xa = xa - mean(xa);
ya      = [0 : Ny-1]*dy; ya = ya - mean(ya);
[Y,X]   = meshgrid(ya,xa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROBLEM 1: Direct Construction from Uniform K-Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CREATE UNIFORM K-FUNCTION
A = (25*degrees) * ones(Nx,Ny);
R = (2*pi/a) * ones(Nx,Ny);
[KX,KY] = pol2cart(A,R);

% GENERATE GRATING
ER = cos(KX.*X + KY.*Y);

% PLOT K FUNCTION AND GRATING
figure('Color','w');
subplot(2,2,1)
imagesc(xa,ya,KX')
colorbar;
axis equal tight
caxis([0 7]);
title('$\textrm{KX}$','FontSize',15,'Interpreter','LaTex');
xlabel('$\textrm{x}$','FontSize',12,'Interpreter','LaTex')
ylabel('$\textrm{y}$','FontSize',12,'Interpreter','LaTex','Rotation',0);

subplot(2,2,3)
imagesc(xa,ya,KY')
colorbar;
axis equal tight
caxis([0 7]);
title('$\textrm{KY}$','FontSize',15,'Interpreter','LaTex');
xlabel('$\textrm{x}$','FontSize',12,'Interpreter','LaTex')
ylabel('$\textrm{y}$','FontSize',12,'Interpreter','LaTex','Rotation',0)

subplot(2,2,[2 4]);
imagesc(xa,ya,ER')
colorbar;
axis equal tight
title('$\textrm{ER}$','FontSize',15,'Interpreter','LaTex');
xlabel('$\textrm{x}$','FontSize',12,'Interpreter','LaTex')
ylabel('$\textrm{y}$','FontSize',12,'Interpreter','LaTex','Rotation',0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROBLEM 2: Direct Construction from Spatially Variant K-Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CREATE SPATIALLY VARIANT K-FUNCTION
A = exp(-(X.^2 + Y.^2)/(2*a)^2);
A = (45*degrees) * (A>0.5);
R = (2*pi/a) * ones(Nx,Ny);
[KX,KY] = pol2cart(A,R);

% GENERATE GRATING
ER = cos(KX.*X + KY.*Y);

% PLOT K FUNCTION AND GRATING
figure('Color','w');
subplot(2,2,1)
imagesc(xa,ya,KX')
colorbar;
axis equal tight
caxis([0 8]);
title('$\textrm{KX}$','FontSize',15,'Interpreter','LaTex');
xlabel('$\textrm{x}$','FontSize',12,'Interpreter','LaTex')
ylabel('$\textrm{y}$','FontSize',12,'Interpreter','LaTex','Rotation',0);

subplot(2,2,3)
imagesc(xa,ya,KY')
colorbar;
axis equal tight
caxis([0 8]);
title('$\textrm{KY}$','FontSize',15,'Interpreter','LaTex');
xlabel('$\textrm{x}$','FontSize',12,'Interpreter','LaTex')
ylabel('$\textrm{y}$','FontSize',12,'Interpreter','LaTex','Rotation',0)

subplot(2,2,[2 4]);
imagesc(xa,ya,ER')
colorbar;
axis equal tight
title('$\textrm{ER}$','FontSize',15,'Interpreter','LaTex');
xlabel('$\textrm{x}$','FontSize',12,'Interpreter','LaTex')
ylabel('$\textrm{y}$','FontSize',12,'Interpreter','LaTex','Rotation',0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROBLEM 3: Construction using grating phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALL fdder() TO OBTAIN DERIVATIVE MATRICES
NS  = [Nx Ny];
RES = [dx dy];
BC  = [1 1];
[DX,~,DY,~] = fdder(NS,RES,BC);

% GENERATE NECESSARY MATRICES
A = [DX ; DY];
b = [KX(:) ; KY(:)];

% COMPUTE GRATING PHASE
PHI = (A.'*A)\(A.'*b);
PHI = reshape(PHI,Nx,Ny);

% COMPUTE ANALOG GRATING
ERA = cos(PHI);

% PLOT RESULTS
figure('Color','w');
subplot(1,2,1);
imagesc(xa,ya,ER');
axis equal tight
title('$\textrm{Direct Method}$','FontSize',15,'Interpreter','LaTex');
xlabel('$\textrm{x}$','FontSize',12,'Interpreter','LaTex')
ylabel('$\textrm{y}$','FontSize',12,'Interpreter','LaTex','Rotation',0)

subplot(1,2,2);
imagesc(xa,ya,ERA');
axis equal tight
title('$\textrm{Grating Phase Method}$','FontSize',15,'Interpreter'...
    ,'LaTex');
xlabel('$\textrm{x}$','FontSize',12,'Interpreter','LaTex')
ylabel('$\textrm{y}$','FontSize',12,'Interpreter','LaTex','Rotation',0)

