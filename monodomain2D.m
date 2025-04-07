clear all
close all

% condizioni iniziali
% U_rest = -54.387;                 %  HH
U_rest = -85.23;                    %  TT

%Tf=30;           %  HH
Tf=500;         %  TT
nt=5000;
dt=Tf/nt;
T=linspace(0,Tf,nt+1);

Xf=1;
nex=200;
x=linspace(0,Xf,nex+1)';
hx=diff(x);

Yf=1;
ney=200;
y=linspace(0,Yf,ney+1)';
hy=diff(y);

[X,Y]=ndgrid(x,y);

[M,L]=assembleMatrices(hx,hy);

Mat=M+dt*1e-3*L;
%H=chol(Mat);

V=zeros(nex+1,ney+1,nt+1);

V(:,:,1)=U_rest;

%ionicModel=HodgkinHuxley(V,dt);
ionicModel=TenTusscher(V,dt);

Iapp=zeros(size(V(:,:,1)));
which=find(X<0.2 & Y<0.2);
Iapp(which)=1;
%Iapp(:,:)=1;   % stimolo su tutto il dominio

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

    %% FATTORIZZAZIONE
     
    % tic
    % y=H'\rhs;
    % Vn_fact=H\y;
    % toc
    % norm(Mat*Vn_fact-rhs)

    %% PCG

    [Vn,~,~,iter2] = pcg(Mat,rhs,1e-7,1000,[],[],Vold(:));
    % norm(Mat*Vn-rhs)

    % PROTOCOL FOR SPIRAL WAVES
    % if (i==3334)
    %     which=find(Y<0.5 | X<0.5);
    %     Vn(which)=U_rest;
    % end

    V(:,:,i)=reshape(Vn,[nex+1 ney+1]);


end
