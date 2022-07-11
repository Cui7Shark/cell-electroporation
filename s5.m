clc;clear;tic;
global a a1  K  Cap Res Rus;
load('p.mat');
load('Boun.mat');
load('cm_suoyin.mat');
load('Cap.mat');
load('Res.mat');
load('Rus.mat');

p=roundn(p,-15);
k=find(roundn(abs(p(:,2)),-4)~=50);

%%对于细胞膜计算极角并排序
a=size(k,1);
a1=size(cm_suoyin,1);

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


%%
K=k1; 

%%
M = eye(a+a1);            %%0-a 计算电流 a--a1 a1--a2 a2--a3分别计算膜
M(1:a,1:a)=3.90e13* Cap;  %%初始孔密度  （算一算Cap的秩或者行列式 确定乘的系数的大小）--如果M的行列式为0或inf 后面计算微分方程出错   
opts = odeset('Mass',M);  %%创建结构体

%%
%%%%初始条件
u0(1:a)=0;              
u0(k1(:,1))=-86e-3;       %%细胞膜内膜电位
% u0(k2(:,1))=-1e-15;     %%细胞核膜内膜电位 原为0 但因计算机计算精度问题 设一极小数字
% u0(k3(:,1))=-174e-3;    %%细胞器膜内膜电位
u0(a+1:a+a1)=1.5e9;       %%初始孔隙密度N0

%%
%%解微分方程
tspan = [0 1000e-09]; %%0-1000ns
[t,u] = ode15s(@transmem,tspan,u0,opts);
%
disp('完成S5');
toc
 



