clear all
close all

% condizioni iniziali
% U_rest = -54.387;                 %  HH
U_rest = -85.23;                    %  TT

%Tf=30*10;      %  HH
Tf=500;         %  TT
nt=5000;
dt=Tf/nt;
T=linspace(0,Tf,nt+1);

Xf=1;
nex=100;
x=linspace(0,Xf,nex+1)';
hx=diff(x);
x_B=linspace(Xf/2,Xf/2,1)';
x_2L=linspace(0,Xf/2-hx(1),nex/2)';
hx_2L=diff(x_2L);
x_2U=linspace(Xf/2+hx(1),Xf,nex/2)';
hx_2U=diff(x_2U);

Yf=2;
ney=100;
y=linspace(0,Yf,ney+1)';
hy=diff(y);
y_B=linspace(Yf/2,Yf/2,1)';
y_2L=linspace(0,Yf/2-hy(1),ney/2)';
hy_2L=diff(y_2L);
y_2U=linspace(Yf/2+hy(1),Yf,ney/2)';
hy_2U=diff(y_2U);

[X,Y]=ndgrid(x,y);

[M,L]=assembleMatrices(hx,hy);
[M1,L1]=assembleMatrices(hx_2L,hy_2L);
[M2,L2]=assembleMatrices(hx_2U,hy_2L);
[M3,L3]=assembleMatrices(hx_2L,hy_2U);
[M4,L4]=assembleMatrices(hx_2U,hy_2U);

return; 

Mat=M+dt*1e-3*L;
H=chol(Mat);

V=zeros(nex+1,ney+1,nt+1);

V(:,:,1)=U_rest;

%ionicModel=HodgkinHuxley(V,dt);
ionicModel=TenTusscher(V,dt);

Iapp=zeros(size(V(:,:,1)));
which=find(X<0.2 & Y<0.2);
Iapp(which)=1;
%Iapp(:,:)=1;

for i=2:nt+1
    t=T(i);
    if mod(i,100)==0
        disp(['t=', num2str(t),'Vmax=',num2str(max(max(V(:,:,i-1))))])
    end

    Vold=V(:,:,i-1);

    Iion=ionicModel.getCurr(Vold,i);

    if (5<t && t<5.3)
        Iapp2=280*Iapp;
    else
        Iapp2=0*Iapp;
    end

    Itot=Iapp2-Iion;

    rhs=Vold+dt*Itot;
    rhs=rhs(:);
    rhs=M*rhs;

    y=H'\rhs;
    Vn=H\y;

    % PROTOCOL FOR SPIRAL WAVES
    % if (i==3334)
    %     which=find(Y<0.5 | X<0.5);
    %     Vn(which)=U_rest;
    % end

    V(:,:,i)=reshape(Vn,[nex+1 ney+1]);

end
