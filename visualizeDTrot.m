function [] = visualizeDTrot(C, M, speed, Dir, time)

%% C: the covariance matrix.

%% Dir: direction of the estimate, to be plotted together with the DT ellipsoid.

%% M: the mean vector, usually 0 in case of DT.

%% speed: time to pause between plotting, lowervalue = faster.

%% Dir: 1 or -1 for clockwise, anticlockwise.

%% time: number of iterations for which rotation is needed, higher = longer.

%% example: visualizeDTrot(diag([17 2 2]),[0 0 0],0.4,1,100)
C = [6, 1, 0; 1, 2, 0; 0, 0, 1];
M = [0, 0, 0];
speed = 0.4;
Dir = 1;
time = 100;
[U, L] = eig(C);
% For N standard deviations spread of data, the radii of the eliipsoid will
% be given by N*SQRT(eigenvalues).
N = 1; % choose your own N
radii = N * sqrt(diag(L));
% generate data for "unrotated" ellipsoid
[xc, yc, zc] = ellipsoid(0, 0, 0, radii(1), radii(2), radii(3), 50);
% rotate data with orientation matrix U and center M
a = kron(U(:, 1), xc);
b = kron(U(:, 2), yc);
c = kron(U(:, 3), zc);
data = a + b + c;
n = size(data, 2);
x = data(1:n, :) + M(1);
y = data(n+1:2*n, :) + M(2);
z = data(2*n+1:end, :) + M(3);
% now plot the rotated ellipse
if Dir == 1
    i = 1:time;
else
    i = time:-1:1;
end
% for j = i
mc = meshc(x, y, z);
% set(mc,'linewidth',2);
colormap hsv
shading interp
zdir = [0, 0, 1];
axis equal
% rotate(sc,zdir,j)
% pause(speed)
title('Correlation Ellipsoid Represented by $[a \ h \ b \ f \ g]^T/c$');
axis equal
xlabel('$x_{ij}$');
ylabel('$y_{ij}$');
zlabel('$\tau$');
alpha(0.5);
tightfig;
set(mc(2), 'visible', 'off');
saveas(gcf, '../ellipsoid.pdf');
view([0, 90]);
set(mc(2), 'visible', 'on');
set(mc(1), 'visible', 'off');
% set(mc(2),'showtext','on');
tightfig;
title('Correlation Ellipse at $\tau$ = 0 Represented by $[a \ h \ b]^T/c$');
saveas(gcf, '../ellipse0.pdf');
close;
% end

%% This should do the job!
