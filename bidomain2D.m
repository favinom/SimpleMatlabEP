clear all
close all

% conductivities
Di=3.15e-4;
De=1.35e-3;

% condizioni iniziali
U_rest = -54.387;                 %  HH
%U_rest = -85.23;                 %  TT



Tf=30;           %  HH
%Tf=500;         %  TT
nt=5000;
dt=Tf/nt;
T=linspace(0,Tf,nt+1);

Xf=1;
nex=100;
x=linspace(0,Xf,nex+1)';
hx=diff(x);

Yf=1;
ney=100;
y=linspace(0,Yf,ney+1)';
hy=diff(y);

[X,Y]=ndgrid(x,y);

[M,L]=assembleMatrices_old(hx,hy);

Li=Di*L;
Lie=(Di+De)*L;

Mat=M+dt*Li;
Mat2=Lie; 

% PerPG1 con chol FATTORIZZAZIONE
H=chol(Mat);

% Per PG2 con ilu
% Mat2=[Lie sum(M,2); sum(M,1) 0];
% [L,U]=ilu(Mat2);


V=zeros(nex+1,ney+1,nt+1);
u=zeros(nex+1,ney+1,nt+1);

V(:,:,1)=U_rest;
u(:,:,1)=0;

ionicModel=HodgkinHuxley(V,dt);
%ionicModel=TenTusscher(V,dt);

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
    uold=u(:,:,i-1);

    Iion=ionicModel.getCurr(Vold,i);

    if (5<t && t<5.3)
        Iapp2=280*Iapp;
    else
        Iapp2=0*Iapp;
    end

    Itot=Iapp2-Iion;

    rhs=Vold+dt*Itot;
    rhs=rhs(:);
    rhs=M*rhs-dt*Li*uold(:);


    %% First step: solve the first PDE for V
    
    % FATTORIZZAZIONE
    % y=H'\rhs;
    % Vn=H\y;
    
    %PCG
    [Vn,~,~,iter1] = pcg(Mat,rhs,1e-7,1000,H,H',Vold(:));
    clear rhs
    
    %if (i==3334)
    %    which=find(Y<0.5 | X<0.5);
    %    Vn(which)=U_rest;
    %end

    V(:,:,i)=reshape(Vn,[nex+1 ney+1]);


    %% Second step: solve the first PDE for u
    
    % FATTORIZZAZIONE 
    % rhs=[-Li*Vn; 0];
    % y=L\rhs;
    % un=U\y;
    % un=un(1:end-1);

    % PCG 
    rhs=-Li*Vn;
    [un,~,~,iter2] = pcg(Mat2,rhs,1e-7,1000,[],[],uold(:));

    % PCG con iLU         NON FUNZIONA
    % rhs=[-Li*Vn; 0];
    % [un,~,~,iter2] = pcg(Mat2,rhs,1e-7,1000,L,L',[uold(:);0]);
    % un=un(1:end-1);

    clear rhs

    u(:,:,i)=reshape(un,[nex+1 ney+1]);
    

end
