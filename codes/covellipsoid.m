function [] = covellipsoid(U,L,M,speed,Dir,time)
%% C: the covariance matrix.
%% Dir: direction of the estimate, to be plotted together with the DT ellipsoid.
%% M: the mean vector, usually 0 in case of DT.
%% speed: time to pause between plotting, lowervalue = faster.
%% Dir: 1 or -1 for clockwise, anticlockwise.
%% time: number of iterations for which rotation is needed, higher = longer.
%% example: visualizeDTrot(diag([17 2 2]),[0 0 0],0.4,1,100)
% For N standard deviations spread of data, the radii of the eliipsoid will
% be given by N*SQRT(eigenvalues).
N = 1; % choose your own N
radii = N*sqrt(diag(L));
% generate data for "unrotated" ellipsoid
[xc,yc,zc] = ellipsoid(0,0,0,radii(1),radii(2),radii(3));
% rotate data with orientation matrix U and center M
a = kron(U(:,1),xc); b = kron(U(:,2),yc); c = kron(U(:,3),zc);
data = a+b+c; n = size(data,2);
x = data(1:n,:)+M(1); y = data(n+1:2*n,:)+M(2); z = data(2*n+1:end,:)+M(3);
% now plot the rotated ellipse
if Dir == 1
	i = 1:time;
else 
	i = time:-1:1;
end
for j = i
sc = surf(x,y,z);
colormap copper
shading interp
zdir = [1 0 0];
axis equal
% rotate(sc,zdir,j)
% pause(speed)
title('actual ellipsoid represented by data: C and M')
axis equal
xlabel('x axis ---->')
ylabel('y axis ---->')
zlabel('z axis ---->')
alpha(0.5)
end
%% This should do the job!
