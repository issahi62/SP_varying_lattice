%% SEssion 3 

% INITIALIZE 
clc
close all 
clear 

% FIGURE 
figure('Color', 'w', 'Units', 'normalized', 'OuterPosition',  [ 0 0 1 1])

%%%%%%%%%%%%%%%%%%%%%%% 
%%   DASHBOARD 
%%%%%%%%%%%%%%%%%%%%%%%

%%% GRATING PARAMETER
a    = 1;         %% period of grating 
NRES = 100;        %% resolution

% SPATIAL HARMONICS PARAMETERS 
NP   = 3;  % p in the reciprocal space
NQ   = 3;  % q in the reciprocal space
NumS = 1;  % number of cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   BUILD GRID FOR CALCULATION  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%GRID SIZE

Sx = NumS*a;  % horizontal scale 
Sy = Sx;   % vertical '
dx = a/NRES; 
dy = a/NRES; 
Nx = ceil(Sx/dx); 
Ny = ceil(Sy/dy); 
dx = Sx/Nx; 
dy = Sy/Ny; 

% Generate Axis 
xa = [0:Nx-1]*dx; xa = xa-mean(xa); 
ya = [0:Ny-1]*dy; ya = ya-mean(ya); 

% Meshgrid Calculation

[Y, X] = meshgrid(ya, xa); 

%%%%%%%%%%%%%%%%%%%%%%% 
%%   BUILD UNIT CELL 
%%%%%%%%%%%%%%%%%%%%%%%

% UNIT CELL PARAMETERS 
w = .9*Sx; 
h = w*sqrt(3)/2;  

ny = round(h/dy); 
ny1 = 1+ floor((Ny-ny)/2); 
ny2 = ny1 + ny -1;

% GENERATING UNIT CELL
UC = zeros(Nx, Ny);

for ny = ny1:ny2
    f  = 1- (ny-ny1 + 1)/(ny2-ny1 + 1); 
    nx = round(f*w/dx); 
    nx1 = 1+ floor((Nx-nx)/2); 
    nx2 = nx1 + nx -1;
    
    UC(nx1:nx2, ny) = 1; 
end 

%% SHOW UNIT CELL 
subplot(231);
h = imagesc(xa, ya, UC'); 
h2 = get(h, 'Parent'); 
set(h2, 'YDir', 'normal'); 
colormap('jet'); 
axis equal tight; 
title('UNIT CELL', 'FontSize', 12); 


% GENERATE FFT OF UNIT CELL

UCF = fftshift(fftn(UC))/(Nx*Ny); 

    
%% SHOW  
subplot(232);
h = imagesc(xa, ya, real(UCF)'); 
h2 = get(h, 'Parent'); 
set(h2, 'YDir', 'normal'); 
colormap('jet'); 
axis equal tight; 
title('FFT', 'FontSize', 12); 


p0 = 1+ceil(Nx/2);
q0 = 1+ceil(Ny/2); 
np1 = p0 - floor(NP/2); 
np2 = p0 + floor(NP/2); 
nq1 = q0 - floor(NQ/2); 
nq2 = q0 + floor(NQ/2); 

% TRUNCATION SECTION 
UCT = UCF(np1:np2, nq1:nq2); 

pa = [-floor(NP/2):floor(NP/2)];
qa = [-floor(NQ/2):floor(NQ/2)]; 

KX = 2*pi*pa/a; 
KY = 2*pi*qa/a; 

[KY, KX] = meshgrid(KY, KX); 

% Quick way for IFFT
UCA = zeros(Nx, Ny); 
UCA(np1:np2, nq1:nq2) = UCT; 
UCI = ifft2(ifftshift(UCA))*(Nx*Ny); 

%% SHOW  
subplot(233);
h = imagesc(pa, qa, real(UCI)'); 
h2 = get(h, 'Parent'); 
set(h2, 'YDir', 'normal'); 
colormap('jet'); 
axis equal tight; 
title('Truncated UNIT CELL', 'FontSize', 12);


%% BRUTE FORCE METHOD
UCK = zeros(Nx, Ny); 
for nq = 1:NQ
    for np = 1:NP
        % CALCULATE PHASOR GRATING
        G = exp(1i *(KX(np, nq)*X + KY(np, nq)*Y)); 
        
        % ADD GRATING TO OVERALL UNIT CELL
        UCK = UCK + UCT(np, nq)*G;  
        
      
        % show Truncated FFT 
        %% SHOW  
        subplot(234);
        h = imagesc(pa, qa, real(UCT')); 
        h2 = get(h, 'Parent'); 
        set(h2, 'YDir', 'normal'); 
        colormap('jet'); 
        axis equal tight; 
        title('Truncated FFT', 'FontSize', 12);
        hold on 
        x = pa(np) -0.5 + [ 0 1 1 0 0]; 
        y = qa(nq) -0.5 + [0, 0, 1, 1 0]; 
        line(x,y, 'Color', 'r', 'LineWidth', 2); 
        hold off 
        %% SHOW 
        
        subplot(235);
        h = imagesc(xa, ya, real(G)'); 
        h2 = get(h, 'Parent'); 
        set(h2, 'YDir', 'normal'); 
        colormap('jet'); 
        axis equal tight; 
        title([ 'P = ' num2str(pa(np)) ', Q = ' num2str(qa(nq))]); 
        
        
          subplot(236);
        h = imagesc(pa, qa, ifftshift(real((UCK)'))); 
        h2 = get(h, 'Parent'); 
        set(h2, 'YDir', 'normal'); 
        colormap('jet'); 
        axis equal tight; 
        title(['RECONSTRUCTED UNIT CELL: No: Periods = ' num2str(NumS)], 'FontSize', 12);
        
        % FORCE MATLAB TO DRAW 
        drawnow(); 
    end 
    
    
end 

