clc;clear;tic;
load('p.mat');
[VX,VY] = voronoi(p(:,1),p(:,2));%生成Voronoi图,返回 Voronoi边的二维顶点。
h = plot(VX,VY,'-b',p(:,1),p(:,2),'.r');%绘图

boundary=[-50,50];
xlim(boundary)
ylim(boundary)

%%
%%%%%%%%%%%%%%%%%显示凸包序号%%%%%%%%%
%% MATLAB官方文档 https://ww2.mathworks.cn/help/matlab/math/voronoi-diagrams.html?searchHighlight=voronoi&s_tid=srchtitle_voronoi_2
nump = size(p,1);
plabels = arrayfun(@(n) {sprintf('p%d', n)}, (1:nump)');
hold on
Hpl = text(p(:,1), p(:,2), plabels, 'color', 'r', ...
      'FontWeight', 'bold', 'HorizontalAlignment',...
      'center', 'BackgroundColor', 'none');
%%
%%%%%%%%%%%%%%%%% 计算voronoi单元  %%%%%%%%%%%%%%

dt = delaunayTriangulation(p); %%基于P中的点创建二维Delaunay三角剖分     %以delaunay为基础计算voronoi参数
[V,R] = voronoiDiagram(dt);    %%返回Delaunay三角剖分中点的Voronoi顶点V和Voronoi区域r,
                               %%r 中的每个区域表示围绕某个三角剖分顶点的点，它们比三角剖分中的其他顶点更靠近该顶点

[row,col]=find(abs(V)>50);     %%find 查找非零元素的索引和值  ---->为什么要大于50？答：绝对值>50就超出了边界
 

for i=1:size(R,1)              %%循环长度  size(R,1):返回第一列的长度                               
   if sum(R{i}==row)~=0        %% R{i}给出与点位i相关联的Voronoi顶点的索引判断有没有超出边界的点，如果和不等于0就表明有超过边界的点；    
      del=find(sum(R{i}==row));
        R{i}(del)=[];          %% 然后把这些点置为空数组。
  end
end

a1=size(V,1);                  %% a1=数组V的长度，即维诺图中单元的顶点个数
%%
hori=find(round(abs(p(:,1)),-4)==50 & p(:,2)<50 & p(:,2)>=-50);  
%%p的第一列取绝对值四舍五入，找到-50<=p第二列<50,p第一列=50 的点 （x=+-50,-50<=y<50）
%%hori表示边界框的左右边界

for i=1:size(hori,1)
    m=max(V(R{hori(i)},2));   %%取V的第二列的最大值 
    V(size(V,1)+1,:)=[p(hori(i),1),m];
end  %%这段作用？？？？
%%
vert=find(round(abs(p(:,2)),-4)==50 & p(:,1)<50 & p(:,1)>=-50);
%%vert表示边界框的上下边界(-50<=x<50,y=+-50)

for i=1:size(vert,1)
    m=max(V(R{vert(i)},1));
    V(size(V,1)+1,:)=[m,p(vert(i),2)];
end
%%这段作用？？？？VC单元排序
a2=size(V,1);
Adara(:,1)=(a1:a2)';Adara(:,2:3)=V(a1:a2,:);
Bound=find(round(abs(p(:,2)),-4)==50 | round(abs(p(:,1)),-4)==50);

for i=1:size(Bound,1)
    Adara(:,4)=sqrt((Adara(:,2)-p(Bound(i),1)).^2+(Adara(:,3)-p(Bound(i),2)).^2);
    Adara=sortrows(Adara,4);
    R{Bound(i)}=[R{Bound(i)},Adara(1,1),Adara(2,1)];
end

%%
%%%%%%%%%%%%%显示voronoi单元顶点序号%%%%%%%%%%%%
numv = size(V,1);
vlabels = arrayfun(@(n) {sprintf('V%d', n)}, (2:numv)');
hold on
Hpl = text(V(2:end,1), V(2:end,2), vlabels, ...
      'FontWeight', 'bold', 'HorizontalAlignment',...
      'center', 'BackgroundColor', 'none');
hold off
%%
for i=1:length(V)
    if V(i,1) > 50 && V(i,1) ~= Inf
       V(i,1) = 50;
    end
    if V(i,1) < -50 && V(i,1) ~= Inf
       V(i,1) = -50;
    end
    if V(i,2) > 50 && V(i,2) ~= Inf
       V(i,2) = 50;
    end
    if V(i,2) < -50 && V(i,2) ~= Inf
       V(i,2) = -50;
    end

end
toc;


