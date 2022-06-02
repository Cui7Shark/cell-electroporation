clc;clear;tic;
global a a1 a2 a3 K  Cap Res Rus;
load(                           'p.mat');
load('Boun.mat');
load('cm_suoyin.mat');
load('lm_suoyin.mat');
load('sm_suoyin.mat');
load('Cap.mat');
load('Res.mat');
load('Rus.mat');

p=roundn(p,-15);
k=find(roundn(abs(p(:,2)),-4)~=50);

%%对于细胞膜计算极角并排序
a=size(k,1);
a1=size(cm_suoyin,1);
a2=size(lm_suoyin,1);
a3=size(sm_suoyin,1);

k1=zeros(length(cm_suoyin),4);
for i=1:size(cm_suoyin,1)
    for j = 1:size(cm_suoyin,2)
        k1(i,j) = find(k == cm_suoyin(i,j));
    end
    k1(i,3) = Boun(cm_suoyin(i,1),cm_suoyin(i,2));
    if abs(atand(p(k(k1(i,1)),2)/p(k(k1(i,1)),1))) == 90
        if p(k(k1(i,1)),2)>0
            k1(i,4) = 0;
        elseif p(k(k1(i,1)),2)<0
            k1(i,4) = 180;
        end
    elseif p(k(k1(i,1)),1)>0
        k1(i,4)=-atand(p(k(k1(i,1)),2)/p(k(k1(i,1)),1))+90;
    elseif p(k(k1(i,1)),1)<0
        k1(i,4)=-atand(p(k(k1(i,1)),2)/p(k(k1(i,1)),1))+270;
    end    
end
k1=sortrows(k1,4);

%%对于细胞核膜
p(:,1)=p(:,1)+2;p(:,2)=p(:,2)-2;
k2=zeros(length(sm_suoyin),4);
for i=1:size(lm_suoyin,1)
   for j = 1:size(lm_suoyin,2)
       k2(i,j) = find(k == lm_suoyin(i,j)); 
    end
    k2(i,3) = Boun(lm_suoyin(i,1),lm_suoyin(i,2));
    if abs(atand(p(k(k2(i,1)),2)/p(k(k2(i,1)),1))) == 90
        if p(k(k2(i,1)),2)>0
            k2(i,4) = 0;
        elseif p(k(k2(i,1)),2)<0
            k2(i,4) = 180;
        end
    elseif p(k(k2(i,1)),1)>0
        k2(i,4)=-atand(p(k(k2(i,1)),2)/p(k(k2(i,1)),1))+90;
    elseif p(k(k2(i,1)),1)<0
        k2(i,4)=-atand(p(k(k2(i,1)),2)/p(k(k2(i,1)),1))+270;
    end    
end
k2=sortrows(k2,4);

%%
%%移动坐标
p(:,1)=p(:,1)-2;p(:,2)=p(:,2)+2;
p(:,1)=p(:,1)-4;p(:,2)=p(:,2)+4;

%%
%%对于小细胞器膜

k3=zeros(length(sm_suoyin),4);
for i=1:size(sm_suoyin,1)
    for j = 1 :size(sm_suoyin,2)
        k3(i,j) = find(k == sm_suoyin(i,j));
    end
    k3(i,3) = Boun(sm_suoyin(i,1),sm_suoyin(i,2));
    if abs(atand(p(k(k3(i,1)),2)/p(k(k3(i,1)),1))) == 90
        if p(k(k3(i,1)),2)>0
            k3(i,4) = 0;
        elseif p(k(k3(i,1)),2)<0
            k3(i,4) = 180;
        end
    elseif p(k(k3(i,1)),1)>0
        k3(i,4)=-atand(p(k(k3(i,1)),2)/p(k(k3(i,1)),1))+90;
    elseif p(k(k3(i,1)),1)<0
        k3(i,4)=-atand(p(k(k3(i,1)),2)/p(k(k3(i,1)),1))+270;
    end    
end
k3=sortrows(k3,4);
p(:,1)=p(:,1)+4;p(:,2)=p(:,2)-4;
%%
K=[k1;k2;k3]; 

%%
M = eye(a+a1+a2+a3);      %%0-a 计算电流 a--a1 a1--a2 a2--a3分别计算膜
M(1:a,1:a)=6.00e13* Cap;  %%初始孔密度
opts = odeset('Mass',M);  %%创建结构体

%%
%%%%初始条件
u0(1:a)=0;              
u0(k1(:,1))=-86e-3;    %细胞膜内膜电位
u0(k2(:,1))=-1e-15;    %%细胞核膜内膜电位 原为0 但因计算机计算精度问题 设一极小数字
u0(k3(:,1))=-174e-3;   %细胞器膜内膜电位
u0(a+1:a+a1+a2+a3)=1.5e9; %%初始孔隙密度N0

%%
%%解微分方程
tspan = [0 101e-6]; %%0-101微秒
[t,u] = ode15s(@transmem,tspan,u0,opts);
%
disp('完成S5');
toc
 



