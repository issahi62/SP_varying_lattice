%% COMPLETE SPATIALLY-VARYING LATTICE 


% INITIALIZE MATLAB 
clear; 
clc; 
close all; 

% OPEN FIGURE WINDOW 
figure('Color', 'w', 'Units', 'normalized', 'Outerposition', [0 0 1 1]); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UNIT CELL 
a   = 1;        %lattice period 
NSP = 10;       % Number of lattice

% GRID PARAMETERS 
Nxh = 255;      %unit cell grid
Nyh = Nxh; 

Sx = NSP * a ;
Sy = Sx; 

NRESLO = 10;   %number of points per period
NRESHI = 20; 

% FOURIER EXPANSION PARAMETERS
cth = 0.01;  %threshold of fourier co-efficients
NP  = 7; 
NQ  = 7; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE GRIDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOW RESOLUTION GRID

dx = a/NRESLO; 
dy = a/NRESLO; 

Nx = ceil(Sx/dx);  dx = Sx/Nx; 
Ny = ceil(Sy/dy);  dy = Sy/Ny; 

xa = (1:Nx)*dy; 
ya = (1:Ny)*dy; 

[Y, X] = meshgrid(ya, xa); 


% HIGH RESOLUTION GRID
Kmax = (2*pi/a) * [floor(NP/2); floor(NQ/2)]; 
amin = (2*pi)/norm(Kmax); %min lattice for high resolution 
dx2 = amin/NRESHI; 
dy2 = amin/NRESHI; 

Nx2 = ceil(Sx/dx2); 
Ny2 = ceil(Sy/dy2); 

xa2 = linspace(xa(1), xa(Nx), Nx2); 
ya2 = linspace(ya(1), ya(Ny), Ny2); 

dx2 = xa2(2) - xa2(1); 
dy2 = ya2(2) - ya2(1); 

[Y2, X2] = meshgrid(ya2, xa2); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD GRAYSCALE UNIT CELL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UNIT CELL GRID 
dxh = a/Nxh; 
dyh = a/Nyh; 
xah = (0:Nxh-1)*dxh; xah = xah - mean(xah); 
yah = (0:Nyh-1)*dyh; yah = yah - mean(yah);

[YH, XH] = meshgrid(yah, xah);


% DEFINE VERTICES OF THE TRIANGLE (USING VERTICES)

w = .9*a;   %width of triangle 
h = w*sqrt(3)/2; 

v1 = [0; h/2]; 
v2 = [-w/2; -h/2]; 
v3 = [+w/2; -h/2]; 

% SECTION 1 
p1 = v1; 
p2 = v2; 
D1 = (p2(2) - p1(2))*XH - (p2(1)-p1(1))*YH - p2(1)*p1(2) - p2(2)*p1(1); 
D1 = -D1./sqrt((p2(2) - p1(2))^2 + (p2(1) - p1(1))^2); 


% SECTION 2 
p1 = v1; 
p2 = v3; 
D2 = (p2(2) - p1(2))*XH - (p2(1)-p1(1))*YH - p2(1)*p1(2) - p2(2)*p1(1); 
D2 = +D2./sqrt((p2(2) - p1(2))^2 + (p2(1) - p1(1))^2);

% SECTION 3 
p1 = v2; 
p2 = v3; 
D3 = (p2(2) - p1(2))*XH - (p2(1)-p1(1))*YH - p2(1)*p1(2) - p2(2)*p1(1); 
D3 = -D3./sqrt((p2(2) - p1(2))^2 + (p2(1) - p1(1))^2);

% BUILD THE UNIT CELL ON GRAYSCALE
UC = min(D1, D2); 
UC = min(UC, D3);
UC = UC- min(UC(:)); 
UC = UC / max(UC(:)); 

% SHOW IMAGES 

subplot(241); 
h1 = imagesc(xah, yah,  UC'); 
h2 = get(h1, 'Parent'); 
set(h2, 'YDir', 'normal'); 
axis equal tight; 
colorbar; 
colormap('jet'); 
title('UNIT CELL', 'Fontsize', 14); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FILL FRACTION SWEEP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gth_dat = linspace(0, 1, 100); 

% CALCULATE FILL FRACTION DATA

ff_dat = 0*gth_dat; 

for n = 1: length(gth_dat)
    
    % GENERATE BINARY UNIT CELL
    UCB = UC > gth_dat(n); 
    
    %CALCULATE FILL FRACTION 
    ff_dat(n) = sum(UCB(:))/(Nxh*Nyh); 
    
    % SHOW STATUS OF BINARY GRATING
    
  
end 
subplot(242); 
h1 = imagesc(xah, yah,  (UC> gth_dat(n/2))'); 
h2 = get(h1, 'Parent'); 
set(h2, 'YDir', 'normal'); 
axis equal tight; 
colorbar; 
colormap('jet'); 
title('BINARY UNIT CELL', 'Fontsize', 14); 

% FILL FRACTION PLOT
subplot(243); 
plot(gth_dat(1:n), ff_dat(1:n), 'Linewidth', 2.5); 
xlim([gth_dat(1), max(gth_dat)]);
ylim([0, 1]); 
axis equal tight
xlabel('Thereshold Value'); 
ylabel('Fill Fraction'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE FOURIER EXPANSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% CALCULATE FFT
A = fftshift(fft2(UC))/(Nxh*Nyh); 

% TRUNCATION TO TAKE SECTION OF FFT AMPS 

pO = 1+floor(Nxh/2); 
qO = 1+floor(Nyh/2); 


np1 = pO - floor(NP/2); 
np2 = pO + floor(NP/2); 
nq1 = qO - floor(NQ/2); 
nq2 = qO + floor(NQ/2); 

AT = A(np1:np2, nq1:nq2); 

pa = (-floor(NP/2):floor(NP/2)); 
qa = (-floor(NQ/2):floor(NQ/2));

% GRATING VECTOR EXPANSION 
KX = 2*pi*pa/a; 
KY = 2*pi*qa/a; 

[KY, KX] = meshgrid(KY, KX); 

% SHOW TRUNCATED FFT
subplot(244); 
h1 = imagesc(pa, qa,  real(AT)'); 
h2 = get(h1, 'Parent'); 
set(h2, 'YDir', 'normal'); 
axis equal tight; 
colorbar; 
colormap('jet'); 
title('FFT', 'Fontsize', 14); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE INPUT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% PERIOD (LOW RESOLUTION) 
PER = a*ones(Nx, Ny); 
% ORIENTATION (LOW RESOLUTION GRID)
THETA = atan2(Y, X); 

% FILL FACTOR (HIGH RESOLUTION GRID)

FF = (X2- Sx/2).^2 + (Y2-Sy/2).^2; 
FF = FF-min(FF(:)); 
FF = FF/ max(FF(:));

FF = 1- FF; 
FF = 0.05 + .3*FF; 


% SHOW  PERIOD
subplot(245); 
h1 = imagesc(pa, qa,  PER'); 
h2 = get(h1, 'Parent'); 
set(h2, 'YDir', 'normal'); 
axis equal tight; 
colorbar; 
colormap('jet'); 
title('PER', 'Fontsize', 14); 


% SHOW  THETA 
subplot(246); 
h1 = imagesc(pa, qa,  THETA'); 
h2 = get(h1, 'Parent'); 
set(h2, 'YDir', 'normal'); 
axis equal tight; 
colorbar; 
colormap('jet'); 
title('THETA', 'Fontsize', 14); 


% SHOW  FF
subplot(247); 
h1 = imagesc(pa, qa,  FF'); 
h2 = get(h1, 'Parent'); 
set(h2, 'YDir', 'normal'); 
axis equal tight; 
colorbar; 
colormap('jet'); 
title('FF', 'Fontsize', 14); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE SPATIALLY-VARIANT LATTICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GENERATE LIST OF PLANAR GRATINGS
KLIST = [ KX(:)'; KY(:)'];
CLIST = AT(:); 

% TRUNCATE LIST 
cmax  = max(abs(CLIST)); 
ind   = find(abs(CLIST)>cth*cmax); 

CLIST = CLIST(ind); 
KLIST = KLIST(:, ind); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN LOOP OF ITERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NK = length(CLIST); 

SVL = zeros(Nx2, Ny2); 

for nk = 1: NK
    
    % CALCULATE K-FUNCTIONS
    Kx = KLIST(1, nk); 
    Ky = KLIST(2, nk); 
    
    [theta, r] = cart2pol(Kx,Ky); 
    
    [Kx, Ky]   = pol2cart(THETA+theta, r*a./PER); 
    
    % SOLVE FOR GRATING PHASE
    PHI = svlsolve(Kx, Ky, dx, dy); 
    
    % INTERPOLATE TO HIGH RESOLUTION GRID
    
    PHI = interp2(ya, xa', PHI, ya2, xa2'); 
    
    % Calculate Analog Planar GRating 
    
    GA = exp(1i*PHI); 
    
    % Add Planar GRating to Overall Lattice 
    
    SVL = SVL + CLIST(nk).*GA; 
    % REFRESH Graphics.
    %clf; 
     
        subplot(248); 
        h1 = imagesc(pa, qa,  real(GA)'); 
        h2 = get(h1, 'Parent'); 
        set(h2, 'YDir', 'normal'); 
        axis equal tight; 
        colorbar; 
        colormap('jet'); 
        title('PLANAR GRATING', 'Fontsize', 14); 
 figure(2);
        subplot(121); 
        h1 = imagesc(pa, qa,  real(SVL)'); 
        h2 = get(h1, 'Parent'); 
        set(h2, 'YDir', 'normal'); 
        axis equal tight; 
        colorbar; 
        colormap('jet'); 
        title('ANALOG LATTICE', 'Fontsize', 14); 
    
    drawnow(); 
    %end
end

% CLEAN UP NUMERICAL NOISE 

SVL = real(SVL); 
SVL = SVL - min(SVL(:)); 
SVL = SVL/max(SVL(:)); 

% GENERATE BINARY LATTICE

GTH = interp1(ff_dat, gth_dat, FF(:)); 

GTH = reshape(GTH, Nx2, Ny2); 

SVLS = SVL > GTH; 

% SHOW FINAL LATTICE
%clf; 
subplot(122); 
h1 = imagesc(pa, qa,  SVLS'); 
h2 = get(h1, 'Parent'); 
set(h2, 'YDir', 'normal'); 
axis equal tight; 
colorbar; 
colormap('jet'); 
title('ANALOG LATTICE', 'Fontsize', 14); 


% SAVE LATTICE TO FILE 
save SVLATTICE2D SVLS dx2 dy2
