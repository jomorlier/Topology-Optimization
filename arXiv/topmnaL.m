%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
% Modified by G. Raze, under the supervision of Prof. J. Morlier
% April 2017

function x = topmnaL(x,nelx,nely,volfrac,penal0,tol,tolchange)
%tol, criteria for merging + ratio lenght/width beam
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
rhoMax = 1.05;
mMax = nelx*nely*volfrac;
penal = 1;
dp = (penal0-penal)/5;
f1 = figure('units','normalized','position',[0.1,0.25,0.3,0.5]);
f2 = figure('units','normalized','position',[0.6,0.25,0.3,0.5]);
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (Cantilever)
% F = sparse(2*(nely+1)*nelx+nely,1,-1,2*(nely+1)*(nelx+1),1);
% U = zeros(2*(nely+1)*(nelx+1),1);
% fixeddofs = 1:2*(nely+1);
% emptyelts = [];
% fullelts = [];
% alldofs = 1:2*(nely+1)*(nelx+1);
% freedofs = setdiff(alldofs,fixeddofs);
% % DEFINE LOADS AND SUPPORTS (L-shape)
F = sparse(2*(nely+1)*nelx+nely,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = kron(2*nely:2*(nely+1):(nelx+2)*(nely+1),[1,1]) + repmat([1,2],1,nelx/2+1);
emptyelts = repmat(0.5*nely*(nelx+1)+1:1:nely*(0.5*nelx+1),1,0.5*nelx)+kron(0:nely:(0.5*nelx-1)*nely,ones(1,0.5*nely));
fullelts = [];
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
%% INITIALIZE ITERATION
loop = 0;
nodesMerged = 1;
suppress = [];
while nodesMerged
  %% START ITERATION
  dtmin = 0;
  dtmax = min([x(4:5:end);x(5:5:end)]);
  x = constraint(x,mMax);
  [cprev,dc,d] = objective(x,nelx,nely,rhoMax,E0,Emin,KE,U,freedofs,iK,jK,F,edofMat,penal,emptyelts,fullelts);
  dc = rescale(dc);
  dcprev = dc;
  change = 1;
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f\n',loop,cprev,mean(full(d(:))));
  while change > tolchange
    loop = loop + 1;
    correctStep = false;
    while ~correctStep

      dt = 0.5*(dtmax+dtmin);
      [xnew,dcdx] = constraint(x - dt*dc,mMax);
      [c,dc,d] = objective(xnew,nelx,nely,rhoMax,E0,Emin,KE,U,freedofs,iK,jK,F,edofMat,penal,emptyelts,fullelts);
      dc = rescale(dc);
      if c > cprev + 0.0001*(xnew-x)'*dcprev
        dtmax = dt;
        dtmin = 0.8*dtmin;
      elseif (xnew-x)'*dc < 0.5*(xnew-x)'*dcprev
        dtmin = dt;
        dtmax = 1.2*dtmax;
      else
        correctStep = true;
      end
      if ~isempty(suppress)
        correctStep = true;
      end
      if norm(dcdx) > 0
        dc = dc - (dc'*dcdx)/norm(dcdx)^2*dcdx;
      end
    end
    change = abs(cprev-c)/c;
    %change = max(abs(xnew(:)-x(:)));
    suppress = suppressNodes(x,nelx,nely);
    if ~isempty(suppress)
      xnew(suppress) = []; dcdx(suppress) = []; dc(suppress) = [];
    end
    x = xnew;
    cprev = c;
    dcprev = dc;
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,cprev, ...
      mean(full(d(:))),change);
    %% PLOT MEMBERS
     figure(f1)
    plotConfig(x); axis equal; axis([0 nelx 0 nely]);  drawnow;
    %saveas(h1,sprintf('MNAa%d.eps',loop)); % will create FIG1, FIG2,...
    figure(f2)
    colormap(flipud(gray)); imagesc(flipud(reshape(d,nely,nelx))); caxis([0 1]); axis equal; axis off; drawnow;
  %saveas(h2,sprintf('MNAb%d.eps',loop)); % will create FIG1, FIG2,...
    
  end
  %to comment or not, change compliance to compare ! positive effect?
  if penal >= penal0 - 1e-6
    [x,nodesMerged] = mergeNodes(x,tol);
    fprintf('%5i nodes merged\n',nodesMerged);
  else
    penal = penal+dp;
    fprintf('Penalty factor changed to %5i\n',penal);
  end
end

function [c,dc,d] = objective(x,nelx,nely,rhoMax,E0,Emin,KE,U,freedofs,iK,jK,F,edofMat,penal,emptyelts,fullelts)
  %% FE-ANALYSIS
    [d,dd] = density(x,nelx,nely,rhoMax);
    d(emptyelts) = 0; dd(emptyelts,:) = 0;
    d(fullelts) = 1; dd(fullelts,:) = 0;
    sK = reshape(KE(:)*(Emin+d'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = sum((U(edofMat)*KE).*U(edofMat),2);
    c = sum(sum((Emin+d.^penal*(E0-Emin)).*ce));
    dc = -penal*(E0-Emin)*dd'*(d.^(penal-1).*ce);
end

function [d,dd] = density(x,nelx,nely,rhoMax)
  d = zeros(nelx*nely,1);
  dd = zeros(nelx*nely,length(x));
  for in=1:length(x)/5
    n = 5*(in-1)+1;
    ln = sqrt(x(n+3)^2+x(n+4)^2);
    sn = sin(x(n+2));
    cn = cos(x(n+2));
    %% FIND NEIGHBOURS
    limx = [max(min(floor(x(n)-ln),nelx),1) ; max(min(floor(x(n)+ln),nelx),1)];
    limy = [max(min(floor(x(n+1)-ln),nely),1) ; max(min(floor(x(n+1)+ln),nely),1)];
    [xN,yN] = meshgrid((limx(1):limx(2))-0.5,(limy(1):limy(2))-0.5);
    xN = reshape(xN,1,numel(xN));
    yN = reshape(yN,1,numel(yN));
    neigh = [ ((xN-x(n))*cn + (yN-x(n+1))*sn)/x(n+3);
              (-(xN-x(n))*sn + (yN-x(n+1))*cn)/x(n+4)];
    nc = and(abs(neigh(1,:))<1,abs(neigh(2,:))<1);
    ind = reshape((ones(limy(2)-limy(1)+1,1)*(limx(1):limx(2))-1)*nely + ...
      (limy(1):limy(2))'*ones(1,limx(2)-limx(1)+1),(limx(2)-limx(1)+1)*(limy(2)-limy(1)+1),1);
    ind = ind(nc);                                              % Global index of neighbours
    %% COMPUTE DENSITIES
    if(~isempty(ind))
      [f,dfdx,dfdy] = shapeFunction(neigh(:,nc));
      dfdx = dfdx/x(n+3); dfdy = dfdy/x(n+4);
      d(ind) = d(ind) + f;
      dd(ind,n) = - cn*dfdx + sn*dfdy;
      dd(ind,n+1) = - sn*dfdx - cn*dfdy;
      dd(ind,n+2) = (-(xN(nc)-x(n))*sn + (yN(nc)-x(n+1))*cn)'.*dfdx - ...
                    ((xN(nc)-x(n))*cn + (yN(nc)-x(n+1))*sn)'.*dfdy;
      dd(ind,n+3) = - ((xN(nc)-x(n))*cn + (yN(nc)-x(n+1))*sn)'.*dfdx/x(n+3);
      dd(ind,n+4) = ((xN(nc)-x(n))*sn - (yN(nc)-x(n+1))*cn)'.*dfdy/x(n+4);
    end

  end
  %% ASYMPTOTIC DENSITY
%   b = 1/(rhoMax-1)+1;
%   a = rhoMax^b/(rhoMax-1);
%   dd = sparse(diag((a*(1-b)*d.^b + a^2)./((d.^b+a).^2))*dd);
%   d = sparse(a*d./(d.^b+a));
  dd(d<1,:) = diag((1+12*d(d<1).^2-28*d(d<1).^3+15*d(d<1).^4))*dd(d<1,:);
  dd(d>=1,:) = 0;
  dd = sparse(dd);
  d = d+4*d.^3-7*d.^4+3*d.^5;
  d = sparse(min(d,1));
end

function [f,dfdx,dfdy] = shapeFunction(x)
  sx = sign(x(1,:)');
  sy = sign(x(2,:)');
  x = abs(x);

  p1 = x(1,:)<0.5;
  p2 = x(2,:)<0.5;

  wx = zeros(size(x,2),1);
  wy = zeros(size(x,2),1);
  dfdx = zeros(size(x,2),1);
  dfdy = zeros(size(x,2),1);

  wx(p1) = 1-6*x(1,p1).^2+6*x(1,p1).^3;
  wx(~p1) = 2-6*x(1,~p1)+6*x(1,~p1).^2-2*x(1,~p1).^3;
  wy(p2) = 1-6*x(2,p2).^2+6*x(2,p2).^3;
  wy(~p2) = 2-6*x(2,~p2)+6*x(2,~p2).^2-2*x(2,~p2).^3;

  f = abs(wx.*wy);

  dfdx(p1) = (-12*x(1,p1)+18*x(1,p1).^2);
  dfdx(~p1) = (-6+12*x(1,~p1)-6*x(1,~p1).^2);
  dfdx = dfdx.*sx.*wy;
  dfdy(p2) = (-12*x(2,p2)+18*x(2,p2).^2);
  dfdy(~p2) = (-6+12*x(2,~p2)-6*x(2,~p2).^2);
  dfdy = dfdy.*sy.*wx;
end

function dcn = rescale(dc)
  dcn = zeros(size(dc));
  for nn=1:5:length(dc)
    if norm(dc(nn+(0:4)))
      dcn(nn+(0:4)) = dc(nn+(0:4))/norm(dc(nn+(0:4)));
    end
  end
end

function [x,dcdx,suppress] = constraint(x,mMax)
  nm = length(x)/5;
  ind = zeros(2*nm,1);
  suppress = [];
  for i = 1 : nm
    ind(2*i-1) = 5*i-1;
    ind(2*i) = 5*i;
  end
  H = 0.373592899303953*sparse(kron(eye(nm),[0 1; 1 0]));
  m = x(ind)'*H*x(ind);
  cm = mMax - m;
  dcdx = zeros(length(x),1);
  if cm < 0
    dmdx = 2*H*x(ind);
    a = dmdx'*H*dmdx;
    b = dmdx'*H*x(ind)+x(ind)'*H*dmdx;
    x1 = x(ind) + (-b+sqrt(b^2+4*a*cm))/(2*a)*dmdx;
    x2 = x(ind) + (-b-sqrt(b^2+4*a*cm))/(2*a)*dmdx;
    if sum(x1<1) < sum(x2<1)
      x(ind) = x1;
    else
      x(ind) = x2;
    end
    dcdx(ind) = -dmdx;
  end
end

function suppress = suppressNodes(x,nelx,nely)
  suppress =[];
  for n = 1 : 5: length(x)
    ln = sqrt(x(n+3)^2 + x(n+4)^2);
    if x(n)+ln < 0 || x(n)-ln > nelx || x(n+1)+ln < 0 || x(n+1)-ln > nely ...
        || x(n+3) < 1 || x(n+4) < 1
      suppress = [suppress ; n + (0:4)];
    end
  end
end

function [x,nodesMerged] = mergeNodes(x,tol)
  nodesMerged = 0;
  i=1;
  while i<length(x)
    j=i+5;
    while j<length(x)
      theta = atan2(x(j+1)-x(i+1),x(j)-x(i));
      [lli,lti] = directionalComponents(x(i:i+4),theta);
      [llj,ltj] = directionalComponents(x(j:j+4),theta);
      if (x(j+1)-x(i+1))^2+(x(j)-x(i))^2 <= (tol(1)*0.37*(lli+llj))^2 && abs(lli/lti-llj/ltj) <= tol(2)
        xe = newElement(x(i:i+4),x(j:j+4),theta,lli+llj,0.5*(lti+ltj));
        x([i:i+4,j:j+4]) = [];
        x = [x;xe];
        nodesMerged = nodesMerged+1;
        i=1;j=1;
      end
      j=j+5;
    end
    i=i+5;
  end
end

function [ll,lt] = directionalComponents(x,theta)
  p = tan(theta-x(3));
  ll = min(abs(x(4)*sqrt(1+p^2)),abs(x(5)*sqrt(1+(1/p)^2)));
  lt = min(abs(x(4)*sqrt(1+(1/p)^2)),abs(x(5)*sqrt(1+p^2)));
end

function xe = newElement(xi,xj,theta,lx,ly)
  xe = zeros(5,1);
  mi = xi(4)*xi(5);
  mj = xj(4)*xj(5);
  me = mi+mj;
  xe(1:2) = mi/me*xi(1:2) + mj/me*xj(1:2);
  xe(3) = theta;
  xe(4:5) = sqrt(me/(lx*ly))*[lx;ly];
end

function plotConfig(x)
  fact = 0.5;%.37
  clf
  hold on
  for i=1:length(x)/5
    n = 5*(i-1)+1;
    sn = sin(x(n+2));
    cn = cos(x(n+2));
    fill([x(n)+fact*x(n+3)*cn-fact*x(n+4)*sn,...
        x(n)-fact*x(n+3)*cn-fact*x(n+4)*sn,...
        x(n)-fact*x(n+3)*cn+fact*x(n+4)*sn,...
        x(n)+fact*x(n+3)*cn+fact*x(n+4)*sn,...
        x(n)+fact*x(n+3)*cn-fact*x(n+4)*sn],...
        [x(n+1)+fact*x(n+3)*sn+fact*x(n+4)*cn,...
        x(n+1)-fact*x(n+3)*sn+fact*x(n+4)*cn,...
        x(n+1)-fact*x(n+3)*sn-fact*x(n+4)*cn,...
        x(n+1)+fact*x(n+3)*sn-fact*x(n+4)*cn,...
        x(n+1)+fact*x(n+3)*sn+fact*x(n+4)*cn],[0 0.2 0.4])
  end
  hold off
end


end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,%
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
%  Technical University of Denmark,                                        %
%  DK-2800 Lyngby, Denmark.                                                %
% Please sent your comments to: sigmund@fam.dtu.dk                         %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


