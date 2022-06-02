clc;clear;close all;tic;
load('p.mat');
load('Adj.mat');
load('Boun.mat');
p=roundn(p,-8);
%%Adj.mat存储邻接矩阵，Boun.mat存储点间距离
Th=10;          %%模型深度
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%初始条件%%%%%%%%%%%%%%%%%%%%%%%%
%
%  sige sigi 1.2  电介质电导率
% sigm  9.5*10^(-9)     膜电导率
% epse epsi  电介质介电常数
% epsm       膜介电常数
Mem_d = 0.005;  %%膜厚5nm
sigma_ele=1.2;         
sigma_mem=9.5*1e-9;
epslong_ele=7.08*1e-10;
epslong_mem=4.43*1e-11;

%%
%%cm 在膜上的P点的索引列向量
cm_nei_suoyin=find(roundn(sqrt(p(:,1).^2+p(:,2).^2),-4)==10);
cm_wai_suoyin=find(roundn(sqrt(p(:,1).^2+p(:,2).^2),-4)==10.005);
%%lom
lm_nei_suoyin=find(roundn(sqrt((p(:,1)+2).^2+(p(:,2)-2).^2),-4)==3);
lm_wai_suoyin=find(roundn(sqrt((p(:,1)+2).^2+(p(:,2)-2).^2),-4)==3.005);
%%som
sm_nei_suoyin=find(roundn(sqrt((p(:,1)-4).^2+(p(:,2)+4).^2),-4)==0.5);
sm_wai_suoyin=find(roundn(sqrt((p(:,1)-4).^2+(p(:,2)+4).^2),-4)==0.505);

%%%CM 
 cm_suoyin =zeros(length(cm_nei_suoyin),2);
for i=1:size(cm_nei_suoyin,1)
    for j=1:size(cm_wai_suoyin,1)
         D=sqrt((p(cm_nei_suoyin(i,1),1)-p(cm_wai_suoyin(j,1),1))^2+(p(cm_nei_suoyin(i,1),2)-p(cm_wai_suoyin(j,1),2))^2);%% 求D膜内外两点间的距离即等于膜的厚度5nm
        if roundn(D,-4) == Mem_d
            cm_suoyin(i,1)=cm_nei_suoyin(i);
            cm_suoyin(i,2)=cm_wai_suoyin(j);%Mem 400*2的数组 第一列为内膜点索引，第二列为外膜点索引 //(i,j)是一对在膜上的节点
            break
        end
    end
end


%%%%LOM
lm_suoyin =zeros(length(lm_nei_suoyin),2);
for i=1:size(lm_nei_suoyin,1)
    for j=1:size(lm_wai_suoyin,1)
         D=sqrt((p(lm_nei_suoyin(i,1),1)-p(lm_wai_suoyin(j,1),1))^2+(p(lm_nei_suoyin(i,1),2)-p(lm_wai_suoyin(j,1),2))^2);
        if roundn(D,-4) == Mem_d
            lm_suoyin(i,1)=lm_nei_suoyin(i);
            lm_suoyin(i,2)=lm_wai_suoyin(j);
            break
        end
    end      
end


%%%SOM
sm_suoyin =zeros(length(sm_nei_suoyin),2);
for i=1:size(sm_nei_suoyin,1)
    for j=1:size(sm_wai_suoyin,1)
         D=sqrt((p(sm_nei_suoyin(i,1),1)-p(sm_wai_suoyin(j,1),1))^2+(p(sm_nei_suoyin(i,1),2)-p(sm_wai_suoyin(j,1),2))^2);
        if roundn(D,-4)==Mem_d
            sm_suoyin(i,1)=sm_nei_suoyin(i);
            sm_suoyin(i,2)=sm_wai_suoyin(j);
            break
        end
    end
end

%%
%%%%%%%%%%%%%%%%%分配导电系数和介电常数%%%%%%%%%%%%%%%%%%%%%
sig=zeros(size(p,1),1);
eps=zeros(size(p,1),1);

for i=1:size(p,1)
    if ((size(find(cm_suoyin == i),1) == 1) || (size(find(lm_suoyin == i),1) == 1) || (size(find(sm_suoyin == i),1) == 1))
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