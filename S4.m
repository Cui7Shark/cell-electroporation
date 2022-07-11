clc;clear;close all;tic;
load('p.mat');
load('Adj.mat');
load('Boun.mat');
load('cm_n.mat');
load('cm_w.mat');
p=roundn(p,-100);
cm_n=roundn(cm_n,-100);
cm_w=roundn(cm_w,-100);
%%Adj.mat存储邻接矩阵，Boun.mat存储点间距离
Th=16.8* 2^(1/4);          %%模型深度 取胞体中心至边缘最远处的距离 即R(pi/6)微米
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%初始条件%%%%%%%%%%%%%%%%%%%%%%%%
%
% sigma_ele  1.2  电介质电导率
% sigma_mem  9.5*10^(-9)     膜电导率
% epse epsi  电介质介电常数
% epsm       膜介电常数
%%
Mem_d = 0.005;  %%膜厚5nm
sigma_ele=1.2;
sigma_mem=9.5*1e-9;
epslong_ele=7.08*1e-10; %电介质介电常数
epslong_mem=4.43*1e-11; %细胞质介电常数

%%
%%cm 在膜上的P点的索引列向量
cm_nei_suoyin = zeros(length(cm_n),1);
for i = 1:length(p)
   for j=1:length(cm_n)
    if p(i,:) == cm_n(j,:)
        cm_nei_suoyin(i,1) = i;
    end
   end
end
cm_nei_suoyin(cm_nei_suoyin == 0) =[];

cm_wai_suoyin = zeros(length(cm_w),1);
for i = 1:length(p)
   for j=1:length(cm_w)
    if p(i,:) == cm_w(j,:)
        cm_wai_suoyin(i,1) = i;
    end
   end
end
cm_wai_suoyin(cm_wai_suoyin == 0) =[];
%%计算膜上点的索引 --- 便于下面给区域分配介电参数和电导率
%%%CM
cm_suoyin =zeros(length(cm_nei_suoyin),2);
for i=1:size(cm_nei_suoyin,1)
    for j=1:size(cm_wai_suoyin,1)
        D=sqrt((p(cm_nei_suoyin(i,1),1)-p(cm_wai_suoyin(j,1),1))^2+(p(cm_nei_suoyin(i,1),2)-p(cm_wai_suoyin(j,1),2))^2);%% 求D膜内外两点间的距离即等于膜的厚度5nm
        if roundn(D,-4) == Mem_d
            cm_suoyin(i,1)=cm_nei_suoyin(i);
            cm_suoyin(i,2)=cm_wai_suoyin(j);%Mem n*2的数组 第一列为内膜点索引，第二列为外膜点索引 //(i,j)是一对在膜上的节点
            break
        end
    end
end

%%
%%%%%%%%%%%%%%%%%分配导电系数和介电常数%%%%%%%%%%%%%%%%%%%%%
sig=zeros(size(p,1),1);
eps=zeros(size(p,1),1);

for i=1:size(p,1)
    %if ((size(find(cm_suoyin == i),1) == 1) || (size(find(lm_suoyin == i),1) == 1) || (size(find(sm_suoyin == i),1) == 1))
    if (size(find(cm_suoyin == i),1) == 1)
        sig(i,1) = sigma_mem;
        eps(i,1) = epslong_mem;
    else
        sig(i,1) = sigma_ele;
        eps(i,1) = epslong_ele;
    end
end


%%
%%%%%%%%%%%%%为VC单元分配电阻值电容值%%%%%%%%%%%%%%%%%%%%%
%%Y---->导纳
Res_VC=zeros(size(Boun));Cap_VC=zeros(size(Boun));
for i = 1:size(Adj,1)
    for j = 1:size(Adj{i},2)
        sigma_av = (sig(i) + sig(Adj{i}(j)))/2;
        epslong_av = (eps(i) + eps(Adj{i}(j)))/2;
        l_ij=sqrt((p(i,1)-p(Adj{i}(j),1))^2+(p(i,2)-p(Adj{i}(j),2))^2);
        Cap_VC(i,Adj{i}(j)) =(epslong_av * (Boun(i,Adj{i}(j)) * 1e-6) * (Th * 1e-6))/(l_ij * 1e-6);
        Res_VC(i,Adj{i}(j)) = (sigma_av * (Boun(i,Adj{i}(j)) * 1e-6) * (Th * 1e-6))/(l_ij * 1e-6) ;

    end
end

%%
%%%%%%%%%%%%%节点导纳矩阵%%%%%%%%%%%%%%%%%%%%%%%

k2=find(roundn(abs(p(:,2)),-4)~=50); %%K1是上下边界 K2是内部
Res=zeros(size(k2,1));
Cap=zeros(size(k2,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%《推导正负号》%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(k2,1)
    for j=1:size(k2,1)
        if i==j    %%i=j 主对角线上，求自导纳
            Res(i,j)=sum(Res_VC(k2(i),:));
            Cap(i,j)=-sum(Cap_VC(k2(i),:)); %%+-号的意义？节点电压法推导
        else     %%求互导纳 加-号
            Res(i,j)=-Res_VC(k2(i),k2(j));
            Cap(i,j)=Cap_VC(k2(i),k2(j));
        end
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%电源导纳矩阵%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1=find(roundn(p(:,2),-4)==50);k2=find(roundn(abs(p(:,2)),-4)~=50);
Rus=zeros(size(k2,1),1);
for i=1:size(k2,1)           %判断Res1中第k2行第k1列元素是否为0,
    for j=1:size(k1,1)
        if Res_VC(k2(i),k1(j))~=0
            Rus(i)=Rus(i) - Res_VC(k2(i),k1(j));  %%排除膜上临近节点电流的影响
        end
    end
end
%%

disp('完成！step_4');
toc