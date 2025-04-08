clear all
close all

% parametri
g_Na = 120.0;  % maximum conductances, in mS/cm^2
g_K = 36.0;
g_L = 0.3;
U_Na = 120;  % reversal potentials, in mV
U_K = -77.0; 
U_L = -54.387;

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

am =@(V)( 0.1 * (25 - V) / (exp((25 - V) / 10) - 1) );
bm =@(V)( 4.0 * exp(-V / 18) );
ah =@(V)( 0.07 * exp(-V / 20) );
bh =@(V)( 1 / (exp((30 - V) / 10) + 1) );
an =@(V)( 0.01 * (10 - V) / (exp((10 - V) / 10) - 1) );
bn =@(V)( 0.125 * exp(-V / 80) );

rhs1=@(m,V)(am(V)*(1-m)-bm(V)*m);
rhs2=@(h,V)(ah(V)*(1-h)-bh(V)*h);
rhs3=@(n,V)(an(V)*(1-n)-bn(V)*n);

% notazione marco ricorda il -
I_Na=@(V,m,h)( g_Na * m^3 * h * (V - U_Na));
I_K =@(V,n)  ( g_K * n^4 * (V - U_K)    );
I_L =@(V)    ( g_L * (V - U_L)          );


Tf=30;
nt=5000;
dt=Tf/nt;
T=linspace(0,Tf,nt+1);

V=zeros(1,nt+1);
m=zeros(1,nt+1);
h=zeros(1,nt+1);
n=zeros(1,nt+1);


V(1)=U_rest;
m(1)=m0;
h(1)=h0;
n(1)=n0;

for i=2:nt+1
    t=T(i);

    if mod(i,1000)
        te=num2str(t);
        Vma=num2str(V(i-1));
        disp(['t=',te,' Vmax=',Vma])
    end

    Vold=V(i-1);
    mold=m(i-1);
    hold=h(i-1);
    nold=n(i-1);
    rhs1eval=rhs1(mold,Vold);
    rhs2eval=rhs2(hold,Vold);
    rhs3eval=rhs3(nold,Vold);
    m(i)=mold+dt*rhs1eval;
    h(i)=hold+dt*rhs2eval;
    n(i)=nold+dt*rhs3eval;

    I_Na_eval=I_Na(Vold,m(i),h(i));
    I_K_eval=I_K(Vold,n(i));
    I_L_eval=I_L(Vold);
    

    if (5<t && t<5.3)
        Iapp=280;
    else
        Iapp=0;
    end

    Itot=Iapp-I_Na_eval-I_K_eval-I_L_eval;

    V(i)=Vold+dt*Itot;

end

plot(T,V)
