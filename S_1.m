clc;clear;close all;tic;

A=[-50,50;-50,-50;50,50;50,-50;0,0];%%BOX:剖分区域
%%
r_c = 10;theta=0:pi/250:2*pi;
cm_nei(:,1) = r_c*cos(theta);
cm_nei(:,2) = r_c*sin(theta);
cm_wai(:,1) = (r_c+0.005)*cos(theta);                              
cm_wai(:,2) = (r_c+0.005)*sin(theta);


%%
r_l = 3;theta2=0:pi/100:2*pi;
Lm_nei(:,1) = -2+r_l*cos(theta2);
Lm_nei(:,2) = 2+r_l*sin(theta2);
Lm_wai(:,1) = -2+(r_l+0.005)*cos(theta2);
Lm_wai(:,2) = 2+(r_l+0.005)*sin(theta2);



%%
r_s= 0.5;theta1=0:pi/30:2*pi;
Sm_nei(:,1) = 4+r_s*cos(theta1);
Sm_nei(:,2) = -4+r_s*sin(theta1);
Sm_wai(:,1) = 4+(r_s+0.005)*cos(theta1);
Sm_wai(:,2) = -4+(r_s+0.005)*sin(theta1);
points= [A;cm_nei;cm_wai;Lm_nei;Lm_wai;Sm_nei;Sm_wai];
%%
%%拟合图形
% figure(1)
% plot(cm_nei(:,1),cm_nei(:,2),'r',cm_wai(:,1),cm_wai(:,2),'g',Lm_nei(:,1),Lm_nei(:,2),'r',Lm_wai(:,1),Lm_wai(:,2),'c',Sm_nei(:,1),Sm_nei(:,2),'m',Sm_wai(:,1),Sm_wai(:,2),'y');
% hold on
%%
%%delaunay三角形剖分


figure(1)
fd=@(p) drectangle(p,-50,50,-50,50);%设定计算区域

fh=@(p) min(min(min(0.003+0.006*abs(dcircle(p,0,0,10)),...
          0.003+0.006*abs(dcircle(p,-2,2,3))),...
          0.0022+0.0044*abs(dcircle(p,4,-4,0.5))),0.3);
[p,t]=distmesh2d(fd,fh,0.1,[-50,-50;50,50],points);


