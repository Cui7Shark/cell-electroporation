function dudt=transmem(t,u)
global a a1 K Res Rus  ;

T=300;
k_Boltzmann=1.38e-23;
q =2.46;
e=1.60e-19;
X=q*e/(k_Boltzmann*T);%%95.0725
rm=0.80e-9;
sigma= 1.2;
dm=5e-9;
w0=2.65;
yita=0.15;
N0=1.5e9;
Vep=0.258;
Th=16.8* 2^(1/4);
alpha = 1e9;
dudt=zeros(a+a1,1);
%%10ns   2MV/m  plus
Vapp=  0.*(t<0)+((400/3)*1e9.*t).*(t<1.5e-9 & t>=0) + 200.*(t>=1.5e-9 & t<8.5e-9) + (-(400/3)*1e9.*t+4000/3).*(t>=8.5e-9 & t<10e-9)+ 0.*(t>=10e-9);
t
%%
for i = 1:(a+a1)
    [row,col]=find(K(:,1:2)==i);                                           %%判断点i是否为膜上节点
    if size(col,1)==0 && i<=a                                              %%非膜上节点的电流
        
        
        dudt(i) = 3.90e13*Res(i,:)*u(1:a) + 3.90e13 *Rus(i)*Vapp;                           %求电流
    elseif  size(col,1)==1 && i<=a
        Vm=@(u) u(K(row,2))-u(K(row,1));                                   %外电势-内电势
        vm=@(u) Vm(u)*X;                                                   %X=q*e/(k_Boltzmann*T);%%95.0725

        Gp=@(u) (sigma*pi*rm*rm/dm) * (exp(vm(u))-1)/((w0*exp(w0-yita*vm(u))-yita*vm(u))*exp(vm(u))/(w0-yita*vm(u))- ...
            ((w0*exp(w0+yita*vm(u))+yita*vm(u))/(w0+yita*vm(u))));
     
%%计算跨膜电流
        %细胞膜内节点
        if col==1 
            dudt(i) = 3.90e13*Res(i,:)*u(1:a) - 3.90e13*Gp(u)*u(row+a)*Vm(u)*K(row,3)*Th*1e-12;      %1e-6*1e-6*;
        %细胞膜外-节点
        elseif col==2 
            dudt(i) = 3.90e13* Res(i,:)*u(1:a) + 3.90e13*Gp(u)*u(row+a)*Vm(u)*K(row,3)*Th*1e-12;      %1e-6*1e-6;
        end
    end
%%计算孔隙密度
    if  i>a
        Vm1=@(u) (u(K(i-a,2))-u(K(i-a,1)));
        dudt(i) = alpha * (exp((Vm1(u)/Vep)^2) * (1-((u(i)/N0)*exp(-q*(Vm1(u)/Vep)^2))));
    end
end
