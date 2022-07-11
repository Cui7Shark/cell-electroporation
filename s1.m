clear;clc;close all;
A=[-50,-50;-50,50;50,-50;50,50];
theata = 0:2*pi/(360*3) :2*pi;
R_theata = (abs(cos(3*theata/2)) +abs(sin(3*theata/2))).^(1/2);
cm_n(:,1) = 16.8* R_theata .* cos(theata);%%x
cm_n(:,2) = 16.8* R_theata .* sin(theata);%%y

cm_w(:,1) = (16.8 * R_theata + 0.005) .* cos(theata);
cm_w(:,2) = (16.8 * R_theata + 0.005) .* sin(theata);
pfix = [A;cm_n;cm_w];




%%
% figure(1);
% hold on;
% plot(cm_n(:,1),cm_n(:,2),cm_w(:,1),cm_w(:,2));
% hold off;

%%
figure(1)
hold on;
fd = @(p) drectangle(p,-50,50,-50,50);

fh = @(p) (0.3 + 0.3*abs( sqrt( p(:,1).^2 + p(:,2).^2 ) - 16.8*(abs(cos(3/2*atan(p(:,2) ./ p(:,1)))) + abs(sin(3/2*atan(p(:,2)./p(:,1))))).^(1/2) ));

[p,t] =  distmesh2d(fd,fh,0.3,[-50,-50;50,50],pfix);

patch('vertices', p, 'faces', t, 'facecolor', [.9, .9, .9] );
%癌样细胞的直角坐标方程
% d = sqrt( p(:,1).^2 + p(:,2).^2 ) - 16.8*(abs(cos(3/2*atan(p(:,2) ./ p(:,1)))) + abs(sin(3/2*atan(p(:,2)./p(:,1))))).^(1/2);
% save(p.mat);save(cm_n.mat);save(cm_w.mat);