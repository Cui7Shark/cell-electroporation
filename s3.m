clc;clear; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%显示完成修改后凸包%%%%%%%%%%%%%%%%%%%
load('p.mat');p=roundn(p,-10);  %%四舍五入到小数点左侧第10位数
load('R.mat')
load('V.mat');V=roundn(V,-10);
%%绘制维诺图
[VX,VY] = voronoi(p(:,1),p(:,2));
h = plot(VX,VY,'-b',p(:,1),p(:,2),'.r');
xlim([-50,50]);
ylim([-50,50]);

%%
%%%%%%%%%%%邻接矩阵%%%%%%%%%%%
Adj=cell(length(p),1);
Boun=zeros(length(p));
for i=1:size(p,1)  %x 遍历
    k=1;           %%索引编号
    for j=1:size(p,1) %y 遍历
        if size(intersect(R{i},R{j}),2)==2   %%intersect 设置两个数组R{1},R{2}的交集 找出i=j的点，（对角线） 返回第二列的长度
            A=intersect(R{i},R{j});
            Adj{i,1}(k)=j;
            Boun(i,j)=sqrt((V(A(1),1)-V(A(2),1))^2+(V(A(1),2)-V(A(2),2))^2);%%计算两点间的距离
            if Boun(i,j)==0    %%如果邻接点间的距离为0即这是同一个点  两个单元的边连接一点点就忽略掉
                Adj{i,1}(k)=[];%%则把这个点的邻接矩阵中所表示的值置为0
            elseif Boun(i,j) == Inf
                Boun(i,j) = 50;
            elseif Boun(i,j) == -Inf
                Boun(i,j) = -50;
            else
                k=k+1;
            end
        end
    end
end

toc;
        