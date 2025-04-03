clear all
close all

% parametri
g_Na = 120.0;  % maximum conductances, in mS/cm^2
g_K = 36.0;
g_L = 0.3;
U_Na = 120;  % reversal potentials, in mV
U_K = -77.0; 
U_L = -54.387;

% conductivities
Di=3.15e-4;
De=1.35e-3;

% condizioni iniziali
m0 = 0.0000;
h0 = 0.9998;
n0 = 0.0040;
U_rest = -54.387;


% am = 0.1 * (25 - V_m) / (exp((25 - V_m) / 10) - 1);
% bm = 4.0 * exp(-V_m / 18);
% ah = 0.07 * exp(-V_m / 20);
% bh = 1 / (exp((30 - V_m) / 10) + 1);
% an = 0.01 * (10 - V_m) / (exp((10 - V_m) / 10) - 1);
% bn = 0.125 * exp(-V_m / 80);

am =@(V)( 0.1 * (25 - V) ./ (exp((25 - V) / 10) - 1) );
bm =@(V)( 4.0 * exp(-V / 18) );
ah =@(V)( 0.07 * exp(-V / 20) );
bh =@(V)( 1 ./ (exp((30 - V) / 10) + 1) );
an =@(V)( 0.01 * (10 - V) ./ (exp((10 - V) / 10) - 1) );
bn =@(V)( 0.125 * exp(-V / 80) );

rhs1=@(m,V)(am(V).*(1-m)-bm(V).*m);
rhs2=@(h,V)(ah(V).*(1-h)-bh(V).*h);
rhs3=@(n,V)(an(V).*(1-n)-bn(V).*n);

% notazione marco ricorda il -
I_Na=@(V,m,h)( g_Na * m.^3 .* h .* (V - U_Na));
I_K =@(V,n)  ( g_K * n.^4 .* (V - U_K)    );
I_L =@(V)    ( g_L * (V - U_L)          );


Tf=30;%*10;
nt=5000;%*10;
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

[M,L]=assembleMatrices(hx,hy);

Li=Di*L;
%Le=De*L;
Lie=(Di+De)*L;

Mat=M+dt*Li;
Mat2=Lie+1e-6*speye(size(Lie));

%tic
H=chol(Mat);
%toc
%tic
H2=chol(Mat2);
%toc

V=zeros(nex+1,ney+1,nt+1);
u=zeros(nex+1,ney+1,nt+1);
m=zeros(nex+1,ney+1,nt+1);
h=zeros(nex+1,ney+1,nt+1);
n=zeros(nex+1,ney+1,nt+1);

V(:,:,1)=U_rest;
u(:,:,1)=0;
m(:,:,1)=m0;
h(:,:,1)=h0;
n(:,:,1)=n0;

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
    uold=u(:,:,i-1);
    mold=m(:,:,i-1);
    hold=h(:,:,i-1);
    nold=n(:,:,i-1);
    rhs1eval=rhs1(mold,Vold);
    rhs2eval=rhs2(hold,Vold);
    rhs3eval=rhs3(nold,Vold);
    m(:,:,i)=mold+dt*rhs1eval;
    h(:,:,i)=hold+dt*rhs2eval;
    n(:,:,i)=nold+dt*rhs3eval;

    I_Na_eval=I_Na(Vold,m(:,:,i),h(:,:,i));
    I_K_eval=I_K(Vold,n(:,:,i));
    I_L_eval=I_L(Vold);

    if (5<t && t<5.3)
        Iapp2=280*Iapp;
    else
        Iapp2=0*Iapp;
    end

    Itot=Iapp2-I_Na_eval-I_K_eval-I_L_eval;

    rhs=Vold+dt*Itot;
    rhs=rhs(:);
    rhs=M*rhs-dt*Li*uold(:);

    %disp('1')
    %tic
    y=H'\rhs;
    Vn=H\y;
    %toc

    clear rhs
    
    %if (i==3334)
    %    which=find(Y<0.5 | X<0.5);
    %    Vn(which)=U_rest;
    %end

    V(:,:,i)=reshape(Vn,[nex+1 ney+1]);

    rhs=-Li*Vn;

    y=H2'\rhs;
    un=H2\y;
    %toc
    clear rhs

    media=sum(M,1)*un;
    un=un-media;

    u(:,:,i)=reshape(un,[nex+1 ney+1]);
    

end
