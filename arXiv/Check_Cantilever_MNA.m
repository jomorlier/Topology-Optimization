clear all; close all;
%check first test case with paper Identifying boundaries of topology optimization results using basic parametric features
%here we put x0 close to the solution, the optimizer still oscillate
 nelx = 2*40;          % Number of finite elements along x
 nely = 2*25;          % Number of finite elements along y
 volfrac = 0.5;     % Targeted Volume fraction
% %% Starting configuration
 nX = 10;             % Number of deformable elements along x
 nY = 5;             % Number of deformable elements along y
 aspect= 9;         % maximal difference on aspect ratio (Ly/Lx)
 ratio=1.2;         % tolerence on the node centroide and its actual dimension Lx,Ly
 tolchange=1e-2;
 %because it's difficult to satisfy the mass constraints as the problem
 %size is changing: underestimate of the total weight
 tune=1.2;
%Creation of the initial x0 (regular grid of moving node)
d = sqrt(volfrac*nelx*nely/(0.5*nX*nY));
[xElts,yElts] = meshgrid(linspace(1/(nX+1)*nelx,nX/(nX+1)*nelx,nX),linspace(1/(nY+1)*nely,nY/(nY+1)*nely,nY));
x0 = [reshape(xElts(~isnan(xElts)),1,numel(xElts(~isnan(xElts))));
      reshape(yElts(~isnan(yElts)),1,numel(yElts(~isnan(yElts))));
      zeros(1,numel(xElts(~isnan(xElts))));
      d*ones(2,numel(xElts(~isnan(xElts))))];
x0 = reshape(x0,numel(x0),1); % Vector of size 5: xposition, yposition, beam orientation, beam length, beam width IN the local coordinate of the node
%% case paper87 : initial nodal positions
 x0=[
    25;45;0;60;7;...
    20;35;-pi/4;45;7;... % 0 in angle for perturbation
    20;15;pi/4;45;7;...
    25;5;0;60;7;...
    65;35;-pi/4;45;7;...
    50;35;pi/4;35;7;...
    50;15;-pi/4;45;7;
    65;15;pi/4;45;7];
%call to topmna.m
disp('MNA')
x = topmna(x0,nelx,nely,volfrac*tune,3,[ratio;aspect],tolchange);
%call to top88.m
figure;
disp('SIMP')
top88(nelx,nely,volfrac,3,2,1)
