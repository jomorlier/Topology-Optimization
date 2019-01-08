clear all; close all;
%Second test case
 nelx = 2*20;          % Number of finite elements along x
 nely = 2*20;          % Number of finite elements along y
 volfrac = 0.25;     % Targeted Volume fraction
% %% Starting configuration
 nX = 7;             % Number of deformable elements along x
 nY = 5;             % Number of deformable elements along y
 aspect= 2;         % maximal difference on aspect ratio (Ly/Lx)
 ratio=1.1;         % tolerence on the node centroide and its actual dimension Lx,Ly
 tolchange=1e-2;
 %because it's difficult to satisfy the mass constraints as the problem
 %size is changing: underestimate of the total weight
 tune=1.0;
%Creation of the initial x0 (regular grid of moving node)
d = sqrt(volfrac*nelx*nely/(0.5*nX*nY));
[xElts,yElts] = meshgrid(linspace(1/(nX+1)*nelx,nX/(nX+1)*nelx,nX),linspace(1/(nY+1)*nely,nY/(nY+1)*nely,nY));
rm = and(xElts > 0.5*nelx,yElts > 0.5*nely);
xElts(rm) = nan;
yElts(rm) = nan;
x0 = [reshape(xElts(~isnan(xElts)),1,numel(xElts(~isnan(xElts))));
      reshape(yElts(~isnan(yElts)),1,numel(yElts(~isnan(yElts))));
      zeros(1,numel(xElts(~isnan(xElts))));
      d*ones(2,numel(xElts(~isnan(xElts))))];
x0 = reshape(x0,numel(x0),1); % Vector of size 5: xposition, yposition, beam orientation, beam length, beam width IN the local coordinate of the node
%call to topmnaL.m
disp('MNA')
x = topmnaL(x0,nelx,nely,volfrac*tune,3,[ratio;aspect],tolchange);
%call to top88.m
figure;
disp('SIMP')
top88L(nelx,nely,volfrac,3,2,1)
